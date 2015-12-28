use strict;
use Tkx;
use GD::Graph::bars;
use GD::Graph::hbars;
Tkx::package_require("tile");
Tkx::package_require("BWidget");
#Tkx::package_require('tooltip');
#Tkx::namespace_import("::tooltip::tooltip");
use Cwd;
use DBI;
use Graph;
use Graph::Undirected;
use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request::Common;
use File::Copy qw(copy);
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; 		# number of significant digits to be returned

our $last_update="23 Dec 2015";


system ("cls");
print STDERR "
+-------------------------------------------------------+
|Exogeneous drug target identification tool (EDTI)	|
+-------------------------------------------------------+

 EDTI: Exogeneous Drug Target Identification tool
 Copyright (C) 2015  Aditya Narayan Sarangi, Kanhu Charan Moharana, 
 Shikha Agnihotry, Mohtashim Lohani, Rakesh Aggarwal

 This program is free software: you can redistribute it and/or modify
 it under the same terms as Perl itself. See (http://dev.perl.org/licenses/). 
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
";

############################
#  GLOBAL VARIABLES       ##
############################
our $Tproteome_file;	## Target proteome aa sequence; FASTA format
our $Hproteome_file;	## Host proteome aa sequence; FASTA format
our $Eproteome_file;	## Essential proteome aa sequence; FASTA format
our $PPI_id_map_file; 	## Sequence id mapped to STRING id
our $interactome_file;	## Interaction file , atleast three columns,ProteinA-ProteinB-Score
our ($skip_cdHit,$skip_host_BLAST)=(0,0);	## Trun on if want to skip steps
our $taxon_id;			## NCBI taxon id; trash
our $setup_error;		## Hold error messages if run for the first time;
our $cmd_hide=1;	##System prefre, batch query run cmd; update by reading sys_conf_file
our $wlc_msg_show=1;
our $blast_version=(check_executable_on_PATH("blastp.exe -help")?'blast+':'old_blastall'); ##Auto detect

our $project_name="New_Project";	
our $root_path=	get_my_document_path();		#getcwd();##update later
our $L_root_path = get_my_document_path('L');	## keeping an extra variable; storing path in Unix format
our $installation_path=getcwd();			##dont update; essential data files and folders;


our $front_page_status="File > Create a new project. ";
our $filter_param_settings_file;						#$root_path."/param.txt";
print STDERR "Intallation path:$installation_path\n";

our %sequence_summary=(
	-total_seq=>0,
	-very_short_seq=>0,
	-paralogslogs=>0,
	-host_orthologs=>0,
	-drug_target_homologs=>0,
	-putative_drug_targets=>0,	
);
##Progress bars
our $read_seq_prg=0;
our $cdhit_prg=0;
our $remove_short_seq_prg=0;
our $string_srch_prg=0;

		


## MAKING SCRIPT MORE INDEPENDANT
print STDERR "File setting:\n";

if(! -d "executables"){ print STDERR "\t".$installation_path."/executables/ : directory Not found. Creating one.\n"; mkdir  "./executables", 0755}
if(! -d "local_dat/"){ print STDERR "\t".$installation_path."/local_dat/ : directory Not found. Creating one.\n"; mkdir  "./local_dat", 0755}
if(!(-e "local_dat/sys_conf")){ print STDERR "\tlocal_dat/sys_conf file not found!!!.Creating..\n";
	open (O, "> local_dat/sys_conf") or die "$! local_dat/sys_conf\n";
	print O "SYS_CONF_HIDE_CMD_PROMPT= 1\n";
	print O "SYS_CONF_SHOW_WELCOME_MSG= 1\n";
	print O "SYS_CONF_BLAST_VER= $blast_version\n";		##Change later static dafault value
	close O;
}

##Making it on path
my $EXEC_PATH = win_path($installation_path."/executables");
$EXEC_PATH =~s/\"//g;
$ENV{'PATH'}.=';'.$EXEC_PATH;

## one may copy the directory to a new computer ; so rerun every time
open (O, ">HideCmd.vbs") or die "$! HideCmd.vbs\n";
print O '
Set objShell = CreateObject("WScript.Shell") 
Set objEnv = objShell.Environment("process")
Set sysEnv = objShell.Environment("system")
 
\'What we want to add
PathToAdd = "'.$EXEC_PATH.'" 
\'Set the new Path
objEnv("PATH") = objEnv("PATH")&";"&PathToAdd
\'Now get the bat file and run it without shownfg Command window
objShell.Run """" & WScript.Arguments(0) & """", 0, False	
';
close O;


print STDERR "Reading System preference:";
open (S, "<local_dat/sys_conf") or die "$! local_dat/sys_conf\n";
while(<S>)
{
	chomp; 
	if(/^SYS_CONF_HIDE_CMD_PROMPT=\s+([01]+)/){$cmd_hide=$1;}
	elsif(/^SYS_CONF_SHOW_WELCOME_MSG=\s+([01]+)/){$wlc_msg_show=$1}
	elsif(/^SYS_CONF_BLAST_VER=\s+(\S+)/){ $blast_version =$1}
}
close S;

if($blast_version eq 'old_blastall'){
	if(!check_executable_on_PATH("blastall --help")){ print STDERR "\nERROR: blastall not found on PATH.\nPlease download ncbi blastall (v.2.2.18) and put it in ENVIORNMENT PATH or 'executables' directory.EXIT\n"; exit(1)}
	if(!check_executable_on_PATH("formatdb --help")){ print STDERR "\nERROR: formatdb.exe not found on PATH.\nPlease download ncbi formatdb.exe and put it in ENVIORNMENT PATH or 'executables' directory.EXIT"; exit(1)}
}
elsif($blast_version eq 'blast+'){
	if(!check_executable_on_PATH("blastp.exe -help")){ print STDERR "\nERROR: blastp.exe not found on PATH.\nPlease download ncbi blast+ (v.2.2.20+) and put it in ENVIORNMENT PATH or 'executables' directory.EXIT\n"; exit(1)}
	if(!check_executable_on_PATH("makeblastdb.exe -help")){ print STDERR "\nERROR: makeblastdb.exe not found on PATH.\nPlease download ncbi blast+ (v.2.2.20+) and put it in ENVIORNMENT PATH or 'executables' directory.EXIT"; exit(1)}

}
else{die "ERROR!!! Unknown blast version: $blast_version\nUse old_blastall/blast+\n";}


if(!check_executable_on_PATH("cd-hit -help")){ print STDERR "\nERROR: ch-hit.exe (or cygwin1.dll) not found on PATH.\nPlease download cd-hit.exe and cygwin1.dll; and put it in ENVIORNMENT PATH or 'executables' directory.EXIT"; exit(1)}
if(!check_executable_on_PATH("sqlite3 -version")){ print STDERR "\nERROR: sqlite3.exe not found on PATH\n.Please download sqlite3.exe and put it in ENVIORNMENT PATH or 'executables' directory.EXIT"; exit(1)}


if(! -d "./local_dat/PATHOGENS"){ 
print STDERR "\t/local_dat/PATHOGENS : directory Not found. Creating one.\n\n"; mkdir  "local_dat/PATHOGENS", 0755; 
system ('sqlite3.exe local_dat/PATHOGENS/pathogen_taxonomy.db  "CREATE TABLE taxonomy (	fasta_file VARCHAR(40) Primary Key,	ORG_CODE VARCHAR (10),	species VARCHAR (100),	TaxID INT(10),	SuperKingdom VARCHAR(50),	phylum VARCHAR(50),	class VARCHAR(50),	order_ VARCHAR(50),	family VARCHAR(50),	genus VARCHAR(50),	Description TEXT NULL);"');
$setup_error.="*** Add Pathogen complete proteomes for broad-spectrum \nanalysis using Utility menu\n\n\n";
}
if(! -d "./local_dat/KNOWN_DRUG_TARGETS"){ print STDERR "\t".$installation_path."/local_dat/KNOWN_DRUG_TARGETS : directory Not found. Creating one.\n\n"; mkdir  "local_dat/KNOWN_DRUG_TARGETS", 0755;
open(D,"> local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt" )or die "$! local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt"; close D;
}
if(! -d "./local_dat/GO"){ print STDERR "\t".$installation_path."/local_dat/GO : directory Not found. Creating one.\n\n"; mkdir  "local_dat/GO", 0755;
}
if(!-e "local_dat/GO/GO.db"){
print STDERR "\t".$installation_path."/local_dat/GO/GO.db : file Not found. Creating one.\n\n";
system ('sqlite3 local_dat/GO/GO.db "CREATE TABLE go_term (	GO_id VARCHAR(10) PRIMARY KEY,	ontology_category VARCHAR(20),	term VARCHAR(150));"');
system ('sqlite3 local_dat/GO/GO.db "CREATE TABLE ecoli_go (GO_id VARCHAR(10) ,	ecoli_ac VARCHAR(8));"'); 
$setup_error.="*** Add E.coli GO annotations using Utility menu. Refer manual for details\n\n\n";
}



print STDERR "DONE\n";


our $min_aa_len=50; 
our $cd_hit_identity=60;
our $chk_cdhit=1;
print STDERR "Reading drug Target data:";
our $drug_db_names=read_drugTarget_db("./local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt");#' {All} {drugBank} {PTTD} ';	##read files to update it;
our $ref_drug_db_array=[];
open(G,"./local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt") or die"$! ./local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt"; while(<G>){chomp; push @$ref_drug_db_array,$_;};close G;
$setup_error.="*** Add Drug-target database and annotations using Utility menu. Refer manual for details\n\n\n" if(scalar @$ref_drug_db_array <1);
our $drug_blast_db_names=create_drugTarget_blast_db($ref_drug_db_array,"./local_dat/KNOWN_DRUG_TARGETS");
our $drug_target_annot=read_drugTarget_annot("./local_dat/KNOWN_DRUG_TARGETS");
print STDERR "DONE";

print STDERR "\nReading broad-spectrum data:";
our $broad_spectrum_pathogen_db_sq = win_path($installation_path.'\local_dat\PATHOGENS\pathogen_taxonomy.db');
##'"\"PATHOGENS\" "PATHOGENS\ACIBC" "PATHOGENS\ACIBS" "PATHOGENS\ACIBT" "PATHOGENS\BURM1\""';##read files to update it;
our $broad_spe_species_per_query=0;
our $tax_level="Family";
our $brd_sp_db_levels;
our $ref_brd_sel_db_array=[]; ## array ref ; array saves  tax_levels sel by user			
($brd_sp_db_levels,$ref_brd_sel_db_array) = fetch_tax_names($tax_level);		
our $broad_spectrum_pathogen_db_list=create_broad_spe_db_array($ref_brd_sel_db_array,"./local_dat/PATHOGENS");
print STDERR "DONE\n";

our $PPI_score_cutoff=700;		##STRING score cutoff
our $top_hub_perc=20;
#
our ($e_val_1,$out_fmt_1,$sub_matrix_1,$gap_score_1, $extra_params_BLAST1,$word_size_1,$threshold_1,$perc_identity_1,$blast_prg1)=(0.01,8,"BLOSUM62","11,1","-b 1",3,11,20,0);	
our ($e_val_2,$out_fmt_2,$sub_matrix_2,$gap_score_2, $extra_params_BLAST2,$word_size_2,$threshold_2,$perc_identity_2,$blast_prg2)=(0.0000000001,8,"BLOSUM62","11,1","-b 1",3,11,0,0);
our ($e_val_3,$out_fmt_3,$sub_matrix_3,$gap_score_3, $extra_params_BLAST3,$word_size_3,$threshold_3,$perc_identity_3,$blast_prg3)=(0.01,8,"BLOSUM62","11,1","",3,11,30,0);
our ($e_val_4,$out_fmt_4,$sub_matrix_4,$gap_score_4, $extra_params_BLAST4,$word_size_4,$threshold_4,$perc_identity_4,$blast_prg4)=(0.1,8,"BLOSUM62","11,1","-b 1",3,11,0,0);
our ($e_val_5,$out_fmt_5,$sub_matrix_5,$gap_score_5, $extra_params_BLAST5,$word_size_5,$threshold_5,$perc_identity_5,$blast_prg5)=(5,8,"BLOSUM62","11,1","-b 1",3,11,0,0);

if ($blast_version eq 'blast+'){
	($out_fmt_1,$out_fmt_2,$out_fmt_3,$out_fmt_4,$out_fmt_5)=(6,6,6,6,6);
	($extra_params_BLAST1,$extra_params_BLAST2,$extra_params_BLAST4,$extra_params_BLAST5)=("-num_alignments 1 ","-num_alignments 1 ","-num_alignments 1 ","-num_alignments 1 ");
}



###Variable for Ontology analysis
#params for finding sismilarity with E.coli
our $model_org_ref_proteome = win_path($installation_path."\\local_dat\\GO\\UP000000625_83333.fasta");		##E.coli is used as refernce for all ontology analaysis
our $GO_db = "./local_dat/GO/GO.db";		##contains two tables ecoli_go, go_terms
#our $model_org_ref_proteome_GO = process_GO_db($installation_path."/local_dat/GO/UP000000625_83333.GO.txt");		##E.coli is used as refernce for all ontology analaysis ; stores ref hash

chk_GO_db();		##Checking GO data in tables;


###Variables for subcellular localization
print STDERR "Checking Internet";
if(check_internet()){ print STDERR ": Internet working\n";}
else{print STDERR ": Internet NOT working, Subcellular localization cannot determined\n"; $setup_error.="*** Internet NOT working, Subcellular \nlocalization cannot determined\n\n\n";}



my $about_text="The availability of complete genome sequences of pathogenic bacteria and their protein complements in public domain has made it possible to determine potential drug targets in these pathogens using computer-based in-silico techniques. Intersection of two datasets, namely \n\t(i)a pathogen's subtractive proteome dataset with the host proteome, and \n\t(ii) the pathogen's minimal essential protein dataset, should represent a set of proteins whose manipulation may reasonably be expected to interfere with the pathogen's survival without adversely affecting the host. These proteins could thus act as potential targets for drugs acting against the particular pathogen.\n\nThis program comes with ABSOLUTELY NO WARRANTY; for details http://www.gnu.org/licenses/.
This is free software, and you are welcome to redistribute it under certain conditions.\n\n";
my $citation_text="Sarangi AN et al., Exogeneous drug target identifiation tool.\nPMID:000000";





####Detecting number of processors#####
my $avl_CPU_count=`wmic computersystem get numberofprocessors`;
$avl_CPU_count=~s/NumberOfProcessors//g;
$avl_CPU_count=~s/\n//g;
if($avl_CPU_count==1){
	$avl_CPU_count=	`wmic computersystem get NumberOfLogicalProcessors`;								#
	$avl_CPU_count=~s/NumberOfLogicalProcessors//g;	
	$avl_CPU_count=~s/\n//g;
}
$avl_CPU_count=~s/\D//g;
my $CPU_list;   
foreach(1..$avl_CPU_count){$CPU_list.=$_." ";   } 
my $use_cores=$avl_CPU_count;		##using max no of CPUs by default

if($filter_param_settings_file){
	print STDERR "Loading Parameter file....";
	
	open (P, "$filter_param_settings_file") or die "$!";
	while(<P>)
	{
		if(/^CPU\s+=\s+(\S+)/){  $use_cores=$1}
		if(/^MIN_AA_LEN\s+=\s+(\S+)/){ $min_aa_len=$1}
		if(/^CHK_CD_HIT\s+=\s+(\S+)/){ $chk_cdhit =$1}
		if(/^CD_HIT_IDN\s+=\s+(\S+)/){ $cd_hit_identity=$1}
		if(/^BRD_SPEC_SP_PER_Q\s+=\s+(\S+)/){ $broad_spe_species_per_query=$1}
       
		if(/^E_VAL_1\s+=\s+(\S+)/){ $e_val_1=$1}
		if(/^OUT_FMT_1\s+=\s+(\S+)/){ $out_fmt_1=$1}
		if(/^SUB_MAT_1\s+=\s+(\S+)/){ $sub_matrix_1=$1}
		if(/^GAP_SCOR_1\s+=\s+(\S+)/){ $gap_score_1=$1}
		if(/^EXTRA_PARAM_1\s+=\s+(\S+)/){ $extra_params_BLAST1=$1}
		if(/^WORD_SIZE_1\s+=\s+(\S+)/){ $word_size_1=$1}
		if(/^THRHOLD_1\s+=\s+(\S+)/){ $threshold_1=$1}
		if(/^PERC_IDENTY_1\s+=\s+(\S+)/){ $perc_identity_1=$1}
	
		if(/^E_VAL_2\s+=\s+(\S+)/){ $e_val_2=$1}
		if(/^OUT_FMT_2\s+=\s+(\S+)/){ $out_fmt_2=$1}
		if(/^SUB_MAT_2\s+=\s+(\S+)/){ $sub_matrix_2=$1}
		if(/^GAP_SCOR_2\s+=\s+(\S+)/){ $gap_score_2=$1}
		if(/^EXTRA_PARAM_2\s+=\s+(\S+)/){ $extra_params_BLAST2=$1}
		if(/^WORD_SIZE_2\s+=\s+(\S+)/){ $word_size_2=$1}
		if(/^THRHOLD_2\s+=\s+(\S+)/){ $threshold_2=$1}
		if(/^PERC_IDENTY_2\s+=\s+(\S+)/){ $perc_identity_2=$1}
      
		if(/^E_VAL_3\s+=\s+(\S+)/){ $e_val_3=$1}
		if(/^OUT_FMT_3\s+=\s+(\S+)/){ $out_fmt_3=$1}
		if(/^SUB_MAT_3\s+=\s+(\S+)/){ $sub_matrix_3=$1}
		if(/^GAP_SCOR_3\s+=\s+(\S+)/){ $gap_score_3=$1}
		if(/^EXTRA_PARAM_3\s+=\s+(\S+)/){ $extra_params_BLAST3=$1}
		if(/^WORD_SIZE_3\s+=\s+(\S+)/){ $word_size_3=$1}
		if(/^THRHOLD_3\s+=\s+(\S+)/){ $threshold_3=$1}
		if(/^PERC_IDENTY_3\s+=\s+(\S+)/){ $perc_identity_3=$1}
       
		if(/^E_VAL_4\s+=\s+(\S+)/){ $e_val_4=$1}
		if(/^OUT_FMT_4\s+=\s+(\S+)/){ $out_fmt_4=$1}
		if(/^SUB_MAT_4\s+=\s+(\S+)/){ $sub_matrix_4=$1}
		if(/^GAP_SCOR_4\s+=\s+(\S+)/){ $gap_score_4=$1}
		if(/^EXTRA_PARAM_4\s+=\s+(\S+)/){ $extra_params_BLAST4=$1}
		if(/^WORD_SIZE_4\s+=\s+(\S+)/){ $word_size_4=$1}
		if(/^THRHOLD_4\s+=\s+(\S+)/){ $threshold_4  =$1}
		if(/^PERC_IDENTY_4\s+=\s+(\S+)/){ $perc_identity_4=$1}
		
		if(/^E_VAL_5\s+=\s+(\S+)/){ $e_val_5=$1}
		if(/^OUT_FMT_5\s+=\s+(\S+)/){ $out_fmt_5=$1}
		if(/^SUB_MAT_5\s+=\s+(\S+)/){ $sub_matrix_5=$1}
		if(/^GAP_SCOR_5\s+=\s+(\S+)/){ $gap_score_5=$1}
		if(/^EXTRA_PARAM_5\s+=\s+(\S+)/){ $extra_params_BLAST5=$1}
		if(/^WORD_SIZE_5\s+=\s+(\S+)/){ $word_size_5=$1}
		if(/^THRHOLD_5\s+=\s+(\S+)/){ $threshold_5  =$1}
		if(/^PERC_IDENTY_5\s+=\s+(\S+)/){ $perc_identity_5=$1}
 
		if(/^PPI_THR\s+=\s+(\S+)/){ $PPI_score_cutoff=$1}
		if(/^TOP_HUB_PERC\s+=\s+(\S+)/){ $top_hub_perc=$1}
	     
	}
	
	close P;	
}
###################################################################################

print STDERR "Start GUI...\n";

my $mw = Tkx::widget->new(".");
$mw->g_wm_title("Exogeneous Drug Target Identification Tool");	
$mw->g_wm_resizable(0,0);
$mw->g_wm_geometry( "800x600" );
$mw->g_grid_columnconfigure(0, -weight => 1);
$mw->g_grid_rowconfigure(0, -weight => 1);

my $OS_platform=Tkx::tk_windowingsystem(); ## will return x11, win32 or aqua
# creating menu bar and menu
##It's important to put the following line in your application somewhere before you start creating menus. Without it, each of your menus (on Windows and X11) will start with what looks like a dashed line, and allows you to "tear off" the menu so it appears in its own window.
Tkx::option_add("*tearOff", 0);	
my $menu = $mw->new_menu();
$mw->configure(-menu => $menu);

#create menu object instances
my $file = $menu->new_menu;
my $settings = $menu->new_menu;
my $utils = $menu->new_menu;
my $dwn_str_anal=$menu->new_menu;
my $help = $menu->new_menu;

$menu->add_cascade(-menu => $file, -label => "File",-underline=>0);
$menu->add_cascade(-menu => $settings, -label => "Settings",-underline=>0);
$menu->add_cascade(-menu => $dwn_str_anal, -label => "Downstream analysis",-underline=>0);
$menu->add_cascade(-menu => $utils, -label => "Utilities",-underline=>0,);
$menu->add_cascade(-menu => $help, -label => "Help",-underline=>0,);


my $system = Tkx::widget->new(Tkx::menu($menu->_mpath . ".system"));
$menu->add_cascade(-menu => $system);
if (Tkx::tk_windowingsystem() eq "aqua") {
		$mw->g_bind("<2>", [sub {my($x,$y) = @_; $menu->g_tk___popup($x,$y)}, Tkx::Ev("%X", "%Y")] );
} 
	
#add menu items 
$file->add_command(-label => "Create new Project",-underline=>1, -command =>\&create_project, -accelerator=>"Ctrl+n");
$mw->g_bind("<Control-n>", sub{create_project()} );
#$file->add_command(-label => "Open Project",-underline=>0, -command => sub {});	
#$file->add_command(-label => "Save Project",-underline=>0, -command => sub {});	
$file->add_command(-label => "Quit",-underline=>0, -accelerator=>"Alt+F4",-command => sub {system("del $root_path\\*.tmp"); exit(1);});

$settings->add_command(-label => "Pipeline settings",-underline=>1, -accelerator=>"Ctrl+s",-command =>\&settings);
$mw->g_bind("<Control-s>", sub{settings()} );
$settings->add_command(-label => "Down-stream analysis",-underline=>1, -accelerator=>"Ctrl+d",-command =>\&down_str_anal);
$mw->g_bind("<Control-d>", sub{down_str_anal()} );
$settings->add_command(-label => "Reset setting",-underline=>1, -accelerator=>"Ctrl+r",-command =>\&reset_params);
$mw->g_bind("<Control-r>", sub{reset_params()} );
$settings->add_command(-label => "System preferences",-underline=>1, -accelerator=>"Ctrl+e",-command =>\&sys_preferences);
$mw->g_bind("<Control-e>", sub{sys_preferences()} );



$dwn_str_anal->add_command(-label =>"Broadspectrum analysis", -underline=>1, -accelerator=>"Ctrl+b",-command =>sub {});	## functions updated later
$dwn_str_anal->add_command(-label =>"Compare with known targets", -underline=>1, -accelerator=>"Ctrl+t",-command =>sub {});
$dwn_str_anal->add_command(-label =>"GO analysis",-underline=>1, -accelerator=>"Ctrl+g",-command =>sub {});
$dwn_str_anal->add_command(-label =>"Sub-cellular localization",-underline=>1, -accelerator=>"Ctrl+l",-command =>sub{});
    

$utils->add_command(-label =>"Add pathogens for broadspectrum analysis", -underline=>5, -command =>\&add_to_broad_spectrum_db);
$utils->add_command(-label =>"Add a drug taerget  database", -underline=>6, -command =>\&add_a_drug_target_db);
$utils->add_command(-label =>"Add Ontology database", -underline=>6, -command =>\&add_ecoli_go_db);
$utils->add_command(-label =>"Create PPI mapping file", -underline=>6, -command =>\&util_PPI);
$utils->add_command(-label =>"Fetch sequences by IDs", -underline=>6, -command =>\&tool_fetch_id_seq);



$help->add_command(-label => "Manual",-underline=>0, -command =>\&manual);	
$help->add_command(-label => "About",-underline=>0, -command =>	\&about);
$help->add_command(-label => "Citation",-underline=>0, -command =>	\&citation);

$mw->g_bind("<F1>", sub{manual()} );

my $frnt=$mw->new_ttk__frame(-borderwidth=>12,-relief => "sunken",-width => 600, -height => 800,-padding => "0 5 5 5" );
$frnt->g_grid(-column=>0,-row=>0,-sticky => "nwes");

my $frnt_top=$frnt->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 500);
$frnt_top->g_grid(-column=>0,-row=>0,-sticky=>"nswe");

 Tkx::image_create_photo( "BANER", -file => "banner.gif");
($frnt_top->new_ttk__label(-image=>'BANER'))->g_grid(-column=>0, -row=>0,-sticky=>"nwes",-columnspan=>2, -padx=>60);		

 my $heading = $frnt_top->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool",-justify=>"center",-foreground=>"blue",-font => "Helvetica 16 bold underline");
$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>50);
 
#my $message = $frnt_top->new_ttk__label(-textvariable =>\$front_page_status,-justify=>"left",-foreground=>"red",-font => "Helvetica 12 italic");
#$message->g_grid(-column=>0,-row=>2,-sticky=>"wn",-padx=>50,-pady=>10,-rowspan=>10);


##frnt_bottom frame 
my $frnt_bottom=$mw->new_ttk__frame(-borderwidth=>12,-relief => "sunken",-width => 600, -height => 150,-padding => "0 0 50 5"  );
$frnt_bottom->g_grid(-column=>0,-row=>1,-sticky=>"nwse");

my $run_but;
$run_but=$frnt_bottom->new_button(-text=>"Run program",-width=>18, -state=>"disabled",-command=>sub{ $frnt_top->g_destroy; $frnt->configure(-padding => "0 0 0 0"); main_script(\$frnt, \$run_but); });
$run_but->g_grid(-column=>0, -row=>0,-sticky=>"e");


my $save_options_but=$frnt_bottom->new_button(-text=>"Save options to file",-width=>18, -state=>"disabled",-command=>sub{
my $filter_param_save_file=$L_root_path."/Parameters.txt";
open (P, ">$filter_param_save_file") or die "$!";

print P "
CPU =  $use_cores
MIN_AA_LEN = $min_aa_len
CHK_CD_HIT = $chk_cdhit 
CD_HIT_IDN = $cd_hit_identity
BRD_SPEC_SP_PER_Q = $broad_spe_species_per_query

E_VAL_1 = $e_val_1
OUT_FMT_1 = $out_fmt_1
SUB_MAT_1 = $sub_matrix_1
GAP_SCOR_1 = $gap_score_1
EXTRA_PARAM_1 = $extra_params_BLAST1
WORD_SIZE_1 = $word_size_1
THRHOLD_1 = $threshold_1
PERC_IDENTY_1 = $perc_identity_1
 
E_VAL_2 = $e_val_2
OUT_FMT_2 = $out_fmt_2
SUB_MAT_2 = $sub_matrix_2
GAP_SCOR_2 = $gap_score_2
EXTRA_PARAM_2 = $extra_params_BLAST2
WORD_SIZE_2 = $word_size_2
THRHOLD_2 = $threshold_2
PERC_IDENTY_2 = $perc_identity_2

E_VAL_3 = $e_val_3
OUT_FMT_3 = $out_fmt_3
SUB_MAT_3 = $sub_matrix_3
GAP_SCOR_3 = $gap_score_3
EXTRA_PARAM_3 = $extra_params_BLAST3
WORD_SIZE_3 = $word_size_3
THRHOLD_3 = $threshold_3
PERC_IDENTY_3 = $perc_identity_3

E_VAL_4 = $e_val_4
OUT_FMT_4 = $out_fmt_4
SUB_MAT_4 = $sub_matrix_4
GAP_SCOR_4 = $gap_score_4
EXTRA_PARAM_4 = $extra_params_BLAST4
WORD_SIZE_4 = $word_size_4
THRHOLD_4 = $threshold_4  
PERC_IDENTY_4 = $perc_identity_4

E_VAL_5 = $e_val_5
OUT_FMT_5 = $out_fmt_5
SUB_MAT_5 = $sub_matrix_5
GAP_SCOR_5 = $gap_score_5
EXTRA_PARAM_5 = $extra_params_BLAST5
WORD_SIZE_5 = $word_size_5
THRHOLD_5 = $threshold_5  
PERC_IDENTY_5 = $perc_identity_5

PPI_THR = $PPI_score_cutoff
TOP_HUB_PERC = $top_hub_perc
 
 ";
Tkx::tk___messageBox(-message => "Parameters used in the current analysis \nwere saved to $L_root_path/Parameters.txt", -type=>"ok",-title=>"Success");
 close P;
 });
$save_options_but->g_grid(-column=>2, -row=>0,-sticky=>"e", -padx=>10);

$dwn_str_anal->entryconfigure("Broadspectrum analysis", -command =>sub{ $frnt_top->g_destroy; $frnt->configure(-padding => "0 0 0 0"); broad_spect_run(\$frnt,\$run_but); }); 
$dwn_str_anal->entryconfigure("Compare with known targets", -command =>sub{ $frnt_top->g_destroy; $frnt->configure(-padding => "0 0 0 0"); comp_known_DT(\$frnt,\$run_but); }); 
$dwn_str_anal->entryconfigure("GO analysis", -command =>sub{ $frnt_top->g_destroy; $frnt->configure(-padding => "0 0 0 0"); GO_analysis(\$frnt,\$run_but); }); 
$dwn_str_anal->entryconfigure("Sub-cellular localization", -command =>sub {$frnt_top->g_destroy; $frnt->configure(-padding => "0 0 0 0"); subCellLoc_analysis(\$frnt,\$run_but); });

$mw->g_bind("<Control-b>", sub{broad_spect_run(\$frnt,\$run_but)} );
$mw->g_bind("<Control-t>", sub{comp_known_DT(\$frnt,\$run_but)} );
$mw->g_bind("<Control-g>", sub{GO_analysis(\$frnt,\$run_but)} );
$mw->g_bind("<Control-l>", sub{subCellLoc_analysis(\$frnt,\$run_but)} );

Tkx::update();

welcome_message() if $wlc_msg_show;

if($setup_error){
Tkx::tk___messageBox(-title=>"Setting error: Metadata not found",-message => "Following settings are missing:\n\n\n".$setup_error, -type=>'ok', -icon=>'warning');
}


## Run project 
##Args: ref of front panel Frams, ref of Run button, ref of an global variable to store ref of all seq data;
##Returns: nothing
sub main_script
{
	my $frm=shift;
	my $run_button=shift;
	my $ppi_status="-";
		
	##Top panel destroyed; create a new at 0 0 
	my $frm_top=$$frm->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 500,-padding => "0 0 0 0");
	$frm_top->g_grid(-column=>0,-row=>0,-sticky=>"nswe");
	$file->entryconfigure("Create new Project",-state=>"disabled");
	
	Tkx::image_create_photo( "BANER", -file => "banner.gif");
	($frm_top->new_ttk__label(-image=>'BANER'))->g_grid(-column=>0, -row=>0,-sticky=>"nwes",-columnspan=>2, -padx=>60);	
	my $heading = $frm_top->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool",-justify=>"center",-foreground=>"blue",-font => "Helvetica 16 bold underline");
	$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>50);
	

	my $p = $frm_top->new_ttk__panedwindow(-orient => 'horizontal',);
	$p->g_grid(-column=>0,-row=>1,-pady=>0,-sticky=>"ne", );
	
	
	# first pane, which would get widgets gridded into it:
	my $frm_left = $p->new_ttk__labelframe(-text => "Run progress", -width => 780, -height => 150 );
	my $frm_right = $p->new_ttk__labelframe(-text => "No. of sequence at each filteration steps",-width => 280,); # second pane
	$p->add($frm_left);
	$p->add($frm_right);
	##########LEFT panel#######################
	$frm_left->new_ttk__label(-text=>"1. Process sequence files		:")->g_grid(-column=>0,-row=>1,-padx=>1,-pady=>1,-sticky=>"w");
	my $prg_read_seq=$frm_left->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$read_seq_prg);
	$prg_read_seq->g_grid(-column=>1,-row=>1,-padx=>5,-pady=>1,-sticky=>"w");
	$frm_left->new_ttk__label(-text=>"2. Removing very short sequences	:")->g_grid(-column=>0,-row=>2,-padx=>1,-pady=>1,-sticky=>"w");
	my $prg_rmv_short_seq=$frm_left->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$remove_short_seq_prg);
	$prg_rmv_short_seq->g_grid(-column=>1,-row=>2,-padx=>5,-pady=>1,-sticky=>"w");
	$frm_left->new_ttk__label(-text=>"3. Run CDhit			:")->g_grid(-column=>0,-row=>3,-padx=>1,-pady=>1,-sticky=>"w");
	my $prg_cdHit=$frm_left->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$cdhit_prg);
	$prg_cdHit->g_grid(-column=>1,-row=>3,-padx=>5,-pady=>1,-sticky=>"w");
	$frm_left->new_ttk__label(-text=>"4. Run BLASTp with host proteome	:")->g_grid(-column=>0,-row=>4,-padx=>1,-pady=>1,-sticky=>"w");
	my $prg_blast_host_homolog=$frm_left->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$blast_prg1);
	$prg_blast_host_homolog->g_grid(-column=>1,-row=>4,-padx=>5,-pady=>1,-sticky=>"w");
	
	##########RIGHT panel#################	
	$frm_right->new_ttk__label(-text=>"1. Total number of sequences	:")->g_grid(-column=>0,-row=>1,-padx=>1,-pady=>2,-sticky=>"w");
	my $entry_tot_seq=$frm_right->new_ttk__entry(-width=>5);
	$entry_tot_seq->g_grid(-column=>1,-row=>1,-padx=>8,-pady=>1,-sticky=>"w");
	$frm_right->new_ttk__label(-text=>"2. Number of short sequences(< $min_aa_len a.a) :")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>1,-sticky=>"w");
	my $entry_short_seq=$frm_right->new_ttk__entry(-width=>5);
	$entry_short_seq->g_grid(-column=>1,-row=>2,-padx=>8,-pady=>1,-sticky=>"w");
	$frm_right->new_ttk__label(-text=>"3. Number of paralogous sequences :")->g_grid(-column=>0,-row=>3,-padx=>4,-pady=>1,-sticky=>"w");
	my $entry_paolg_seq=$frm_right->new_ttk__entry(-width=>5);
	$entry_paolg_seq->g_grid(-column=>1,-row=>3,-padx=>8,-pady=>1,-sticky=>"w");
	$frm_right->new_ttk__label(-text=>"4. Number proteins homologous to host-proteome :")->g_grid(-column=>0,-row=>4,-padx=>4,-pady=>1,-sticky=>"w");
	my $entry_host_homolog_seq=$frm_right->new_ttk__entry(-width=>5);
	$entry_host_homolog_seq->g_grid(-column=>1,-row=>4,-padx=>8,-pady=>1,-sticky=>"w");
	$frm_right->new_ttk__label(-text=>"5. Number of proteins non-homologous to \nhost proteome :")->g_grid(-column=>0,-row=>5,-padx=>4,-pady=>1,-sticky=>"w");
	my $entry_non_host_proteome=$frm_right->new_ttk__entry(-width=>5);
	$entry_non_host_proteome->g_grid(-column=>1,-row=>5,-padx=>8,-pady=>1,-sticky=>"w");
	
	#########################################################################################
	##just givng a space between two sections
	$frm_top->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 80,-padding => "20 0 10 0")->g_grid(-column=>0,-row=>2,-sticky=>"nswe",-pady=>2 );
	
	
	################################################################################################
	my $p2 = $frm_top->new_ttk__panedwindow(-orient => 'horizontal');
	$p2->g_grid(-column=>0,-row=>3,-sticky=>"ne",-pady=>2 );
	my $frm_seq_app = $p2->new_ttk__labelframe(-text => "Sequence based approach", -width => 200, -height => 250);
	my $frm_ppi_app = $p2->new_ttk__labelframe(-text => "PPI network based approach", -width => 100, -height => 250); # second pane
	
	$p2->add($frm_seq_app);
	$p2->add($frm_ppi_app);
	
	my $entry_drug_cand_seq_app=$frm_seq_app->new_ttk__entry(-width=>5);
	my $prg_blast_ess_homolog=$frm_seq_app->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$blast_prg2);
	my $do_ess_pro_blast=$frm_seq_app->new_button(-text=>"Perform BlastP search",-width=>25, -state=>"disabled",-command=>sub{	});  ## command will be added later
	
	$do_ess_pro_blast->g_grid(-column=>0,-row=>0,-padx=>5,-pady=>1,-columnspan=>2);
	$frm_seq_app->new_ttk__label(-text=>"Run BLASTp with essential proteins	:")->g_grid(-column=>0,-row=>2,-padx=>1,-pady=>1,-sticky=>"w");
	$prg_blast_ess_homolog->g_grid(-column=>1,-row=>2,-padx=>5,-pady=>1,-sticky=>"w");
	$frm_seq_app->new_ttk__label(-text=>"Number of putative drug targets	:")->g_grid(-column=>0,-row=>3,-padx=>1,-pady=>1,-sticky=>"w");
	$entry_drug_cand_seq_app->g_grid(-column=>1,-row=>3,-padx=>5,-pady=>1,-sticky=>"w");
	
	###Approach -2 PPI	#################################
	my $do_PPI_search=$frm_ppi_app->new_button(-text=>" Identify hub proteins in PPI network",-width=>30, -state=>"disabled",-command=>sub{	}); ## command will be added later
	$do_PPI_search->g_grid(-column=>0,-row=>0,-padx=>5,-pady=>1,-columnspan=>2);
	$frm_ppi_app->new_ttk__label(-text=>"Searching hub proteins in PPI network	:")->g_grid(-column=>0,-row=>2,-padx=>1,-pady=>1,-sticky=>"w");
	my $prg_String_search=$frm_ppi_app->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$string_srch_prg);
	$prg_String_search->g_grid(-column=>1,-row=>2,-padx=>5,-pady=>1,-sticky=>"w");
	$frm_ppi_app->new_ttk__label(-text=>"Number of putative drug targets	:")->g_grid(-column=>0,-row=>3,-padx=>1,-pady=>1,-sticky=>"w");
	my $entry_drug_cand_ppi_app=$frm_ppi_app->new_ttk__entry(-width=>5);
	$entry_drug_cand_ppi_app->g_grid(-column=>1,-row=>3,-padx=>5,-pady=>1,-sticky=>"w");
	$frm_ppi_app->new_ttk__label(-textvariable=>\$ppi_status)->g_grid(-column=>0,-row=>6,-padx=>1,-pady=>5,-sticky=>"e",-columnspan=>2);
	###################panel widget END##########################	
	
	my $save_result = $frm_top->new_button(-text=>"Export putative target sequences",-width=>28, -state=>"disabled",-command=>sub{
			my $save= Tkx::tk___getSaveFile();
			#$crt_win->g_raise();
			#my $r=fetch_seq_by_id($all_t_seq,$novel_drug_targets);	##
			#write_fasta_seq($r,"$save_result") if $save_result ;;
			copy "$L_root_path/accepted_seq_step-4_1.fasta", $save;
				
	});
	$save_result->g_grid(-column=>0,-row=>6,-sticky=>"nw",-padx=>5, -pady=>15 );
	
	$$run_button->configure(-state=>"disabled");
	
	###########################################################
	###MAIN script goes here								  #
	###########################################################

	my $all_t_seq=read_fasta_sequence($Tproteome_file);
	open(R1,"> $L_root_path/excluded_seq_step-1.fasta") or die " $! $L_root_path/excluded_seq_step-1.fasta";
	open(A1,"> $L_root_path/accepted_seq_step-1.fasta") or die "$! $L_root_path/accepted_seq_step-1.fasta";
	$sequence_summary{-total_seq}=scalar keys %$all_t_seq;
	$entry_tot_seq->delete(0, "end"); $entry_tot_seq->insert(0, $sequence_summary{-total_seq}) ;
	$read_seq_prg=100;Tkx::update();
	foreach my $id(sort keys %$all_t_seq)
		{
			if ($all_t_seq->{$id}->{-len} <=$min_aa_len){ print R1 "$all_t_seq->{$id}->{-desc}\   |$all_t_seq->{$id}->{-len}\n$all_t_seq->{$id}->{-seq}\n" ;
			$sequence_summary{-very_short_seq}++;
			}
			else{print A1 "$all_t_seq->{$id}->{-desc}\   |$all_t_seq->{$id}->{-len}\n$all_t_seq->{$id}->{-seq}\n" ;}	##BUG fix
		}
	close R1;
	close A1;
	$entry_short_seq->delete(0, "end"); $entry_short_seq->insert(0, $sequence_summary{-very_short_seq}) ;
	$remove_short_seq_prg=100;Tkx::update();

	

	if(!$skip_cdHit)
		{
			my $cdHit_c=$cd_hit_identity/100;

			unlink "$L_root_path/cdHit_out.clstr.txt";		
			#`$cdhit_path -i $root_path/accepted_seq_step-1.fasta -o $root_path/cdHit_out -d 0 -c $cdHit_c -n 3`;
			`echo echo off > batch.bat`;
			`echo color B0 >> batch.bat`;
			`echo cls >> batch.bat`;
			`echo echo :: External program cd-hit.exe :: >> batch.bat`;
			`echo echo ---------------------------------------------- >> batch.bat`;
			
			`echo echo 	Input		: $root_path\\accepted_seq_step-1.fasta >> batch.bat`;
			`echo echo 	Sequence identity threshold	: $cd_hit_identity\% >> batch.bat`;
			`echo echo 	Word_length	: 3 >> batch.bat`;
			`echo echo 	Other parameters	: Default (refer to cd-hit website) >> batch.bat`;
			
			
			`echo echo Please wait.......... >> batch.bat`;
			`echo cd-hit.exe -i $root_path\\accepted_seq_step-1.fasta -o $root_path\\cdHit_out -d 0 -c $cdHit_c -n 3 >> batch.bat`;
			`echo rename $root_path\\cdHit_out.clstr cdHit_out.clstr.txt >> batch.bat`;	##mv works
			`echo exit >> batch.bat`;
			
			if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
			else{system("start batch.bat ");}
			
			while(!(-e "$L_root_path/cdHit_out.clstr.txt")){Tkx::update(); sleep(1); $cdhit_prg=25;}	##wait till cdHit_out.clstr.txt is available
			
			my ($exclude, $include)=process_cdHit_clstr("$L_root_path/cdHit_out.clstr.txt");
			#print "@$exclude\n$include\n";
			$cdhit_prg=75;Tkx::update();
			my $r=fetch_seq_by_id($all_t_seq,$exclude);
			write_fasta_seq($r,"$L_root_path/excluded_seq_step-2.fasta");
			$sequence_summary{-paralogslogs}=scalar @$exclude;
			my $r=fetch_seq_by_id($all_t_seq,$include);
			write_fasta_seq($r,"$L_root_path/accepted_seq_step-2.fasta");
			$cdhit_prg=95;Tkx::update();
		}
	else{ 
			$cdhit_prg=25;
			system "copy $root_path\\accepted_seq_step-1.fasta $root_path\\accepted_seq_step-2.fasta"; 
			$sequence_summary{-paralogslogs}=0;
		}	##cp works well
		$entry_paolg_seq->delete(0, "end"); $entry_paolg_seq->insert(0, $sequence_summary{-paralogslogs}) ;
		$cdhit_prg=100;	
		Tkx::update();
		
	$Hproteome_file=win_path($Hproteome_file);
	
	my ($host_like_proteins,$not_host_like_proteins);
	if(!$skip_host_BLAST){
	
			unlink "$L_root_path/host_orthologs_blast1.out.txt";
			my($gap_open, $gap_extns)=split /,/,$gap_score_1;
			my $blast1="";
			$blast1="blastall.exe -p blastp -d $Hproteome_file -i $root_path\\accepted_seq_step-2.fasta -e $e_val_1 -m $out_fmt_1 -W $word_size_1 -M $sub_matrix_1 -G $gap_open -E $gap_extns -o $root_path\\host_orthologs_blast1.out -a $use_cores -f $threshold_1"." $extra_params_BLAST1" if $blast_version eq 'old_blastall';
			$blast1="blastp.exe -db $Hproteome_file -query $root_path\\accepted_seq_step-2.fasta -evalue $e_val_1 -outfmt $out_fmt_1 -word_size $word_size_1 -matrix $sub_matrix_1 -gapopen $gap_open -gapextend $gap_extns -out $root_path\\host_orthologs_blast1.out -num_threads $use_cores -threshold $threshold_1"." $extra_params_BLAST1" if $blast_version eq 'blast+'; ##
			
					`echo echo off > batch.bat`;
					`echo color B0 >> batch.bat`;
					`echo cls >> batch.bat`;
					`echo echo :: External program BLAST :: >> batch.bat`;
					`echo echo ---------------------------------------------- >> batch.bat`;
					
					`echo echo Parameters >> batch.bat`;
					`echo echo 	Program: blastp >> batch.bat`;
					`echo echo 	Query		: $root_path\\accepted_seq_step-2.fasta >> batch.bat`;
					`echo echo 	Database	: $Hproteome_file >> batch.bat`;
					`echo echo 	E-value	: $e_val_1 >> batch.bat`;
					`echo echo 	Scoring matrix: $sub_matrix_1 >> batch.bat`;
					`echo echo 	Gap-penalty (Open,Extension): $gap_open,$gap_extns >> batch.bat`;
					`echo echo 	Word size: $word_size_1 >> batch.bat`;
					`echo echo 	Threshold for extending hits: $threshold_1 >> batch.bat`;
					`echo echo 	CPUs: $use_cores >> batch.bat`;
					
					`echo echo 	Identity percentage: $perc_identity_1 >> batch.bat`;
					
					`echo echo 	Please wait..... >> batch.bat`;
					`echo $blast1 >> batch.bat`;
					`echo rename $root_path\\host_orthologs_blast1.out host_orthologs_blast1.out.txt >> batch.bat`;	##mv works
					`echo exit >> batch.bat`;
					if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
					else{system("start batch.bat ");}

			while(!(-e "$L_root_path/host_orthologs_blast1.out.txt")){
					##wait till human_orthologs_blast1.out.txt is available
				$blast_prg1=blast_progress("$root_path\\host_orthologs_blast1.out","$L_root_path/host_orthologs_blast1.out",$sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-orthologous})); 
				$blast_prg1= 90 if $blast_prg1>90;
				sleep(3);Tkx::update(); 
			}	
			$blast_prg1=100;Tkx::update();
			if(count_blast_hits("$L_root_path/host_orthologs_blast1.out.txt")<1){
				Tkx::tk___messageBox(-message => "No BLAST hit found!!!Aborting analysis\nPlease change BLAST parameters and run analysis again.", -type=>"ok", -title=>"Alert",-icon=>'warning'); 
				return ();
				}
			
			($host_like_proteins,$not_host_like_proteins)=process_host_blast_out("$L_root_path/accepted_seq_step-2.fasta","$L_root_path/host_orthologs_blast1.out.txt", $perc_identity_1 );
			$sequence_summary{-host_orthologs}=scalar @$host_like_proteins;
					
			my $non_host_seq_id=fetch_seq_by_id($all_t_seq,$not_host_like_proteins);
			write_fasta_seq($non_host_seq_id,"$L_root_path/accepted_seq_step-3.fasta");
			my $r=fetch_seq_by_id($all_t_seq,$host_like_proteins);
			write_fasta_seq($r,"$L_root_path/excluded_seq_step-3.fasta");
			
	}
	else{
			$blast_prg1=50;
			system "copy $root_path\\accepted_seq_step-2.fasta $root_path\\accepted_seq_step-3.fasta"; 
			$sequence_summary{-host_orthologs}=0;
			my @ids = keys %{read_fasta_sequence("$L_root_path/accepted_seq_step-3.fasta")};
			($host_like_proteins,$not_host_like_proteins)=([],\@ids);		
			$blast_prg1=100;Tkx::update();	
	}
		
	$entry_host_homolog_seq->delete(0, "end"); $entry_host_homolog_seq->insert(0, $sequence_summary{-host_orthologs}) ;
	$entry_non_host_proteome->delete(0, "end"); $entry_non_host_proteome->insert(0, $sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-paralogslogs}+$sequence_summary{-host_orthologs})) ;	
	Tkx::update();	
	if (scalar @$not_host_like_proteins <2){Tkx::tk___messageBox(-message => "All input sequences have homologs in host proteome!!! ", -type=>"warning", -title=>"Warning"); exit(1);}
### Ask user if want to store result and terminate	
	
	Tkx::tk___messageBox(-message => "Select an target prediction approach\nPerform either 'Sequence-based approach'  or 'PPI network-based approach' ", -type=>"ok", -title=>"Alert");
	
###do_PPI_search	
	$do_PPI_search->configure(-command=>sub{
	
			if(!$interactome_file){
				Tkx::tk___messageBox(-type => "ok", 
					-message => "Probably you forgot to load PPI INTERACTION file.\nPress OK to import input.",
					-icon => "error", -title => "ERROR");
				my $interactome_file1=Tkx::tk___getOpenFile();	
				$interactome_file=$interactome_file1;
			}
			
			
			if(!$PPI_id_map_file){
				Tkx::tk___messageBox(-type => "ok",
					-message => "Probably you forgot to load ID MAPPING file.\nPress OK to import input.",
					-icon => "error", -title => "Input file missing");
				my $PPI_id_map_file1=Tkx::tk___getOpenFile();	
				$PPI_id_map_file=$PPI_id_map_file1;
			}
			
			$do_PPI_search->configure(-state=>"disabled");
			$do_ess_pro_blast->configure(-state=>"disabled");
			
			$string_srch_prg=0;
			#$interactome_file=~ s{/}{\\}g;
			#import PPI file to  sql database 25%
			#calculations; 40%
			#id mapping to $L_root_path/accepted_seq_step-3.fasta; 25%
			
			
			#sorting and mapping to accepted_seq_step-2.fasta ids; 10% $non_host_seq_id
			Tkx::tk___messageBox(-type => "ok",
					-message => "Shortest path calculation for all node pairs using Floyd-Warshall's algorithm may freez the GUI window. Please wait...\n",
					-title => "Alert", -icon=>'warning');
			Tkx::update();
			
			open (P,"$PPI_id_map_file") or die "$!";
			my %seq_id_to_PPI_id_map;##seq_id as key
			my %PPI_id_to_seq_id_map;##PPI_id as key
			my @no_PPI_id_found=0; my $ppi_found=0;
			while(<P>){
				chomp; next if (/^#/ or !$_);
				my ($k,$v)=split /\s+/,$_;
				#print N "K,V $k,$v\n";
				if(!$v){push @no_PPI_id_found,$k;}
				else{$PPI_id_to_seq_id_map{$v}=$k;	$seq_id_to_PPI_id_map{$k}=$v; $ppi_found++;}	
			}
			close P;
			#print STDERR "Warning: no ppi found:".$#no_PPI_id_found+1;
			#print STDERR "\nWarning:  ppi found:$ppi_found\n";
			
			my $database = $root_path."\\PPI_sqlite.db";
						
			open(S,">import.sql") or die "$!";
			
			my $bat="
DROP TABLE IF EXISTS PPI;
DROP TABLE IF EXISTS tmp;
CREATE TABLE tmp (id VARCHAR(20) NOT NULL);
CREATE TABLE PPI (proteinA VARCHAR(20) NOT NULL, proteinB VARCHAR(20) NOT NULL, score INT NOT NULL);
.separator ,
.import \"$interactome_file\" PPI";
			print S "$bat";
			close S;

			open(S,">batch.bat") or die "$!";
			print S "echo off \n";
			print S "color B0  \n";
			print S "cls  \n";
			print S "echo :: External program sqlite.exe ::  \n";
			print S "echo -------------------------------------  \n";
			print S "echo Input: $interactome_file  \n";
			print S "echo Database: $database  \n";
			print S "echo Separator: ','  \n";
			print S "echo ---------------------------\necho Please wait....  \n";
			print S "sqlite3.exe $database <import.sql\n";
			#print S "copy import.sql  mm.sql\n";
			print S "del import.sql  \n";
			print S "exit\n";
			close S;
			if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
			else{system("start batch.bat ");}
			
			while(-e "import.sql"){	$string_srch_prg=5;Tkx::update();sleep(1);	}
			$string_srch_prg=35;Tkx::update();

			my $g = Graph::Undirected->new(); # An undirected graph.
			my $driver   = "SQLite";
			my $dsn = "DBI:$driver:dbname=".unix_path($database);
			my $userid = "";
			my $password = "";
			my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 })
								  or die $DBI::errstr;
			$ppi_status="importing IDs"; Tkx::update();	
			my $c=0;
			foreach my $i(@$not_host_like_proteins)
			{				
				my $p=$seq_id_to_PPI_id_map{$i};
				#print "$i --> $p -\n";
				next if !$p;
				my $stmt = qq(INSERT INTO tmp (id) VALUES ('$p'););	
				my $rv = $dbh->do($stmt) or die $DBI::errstr;	
				$ppi_status="importing IDs ".sprintf("%d/%d",$c++,scalar@$not_host_like_proteins); Tkx::update();	
			}
			
			my $stmt = qq(SELECT PPI.proteinA, PPI.proteinB, PPI.score FROM PPI INNER JOIN tmp where (tmp.id=PPI.proteinA) AND PPI.score >=$PPI_score_cutoff);	##LIMIT 1500
			my ($row_c) = $dbh->selectrow_array(qq(SELECT COUNT(*) FROM PPI INNER JOIN tmp where (tmp.id=PPI.proteinA) AND PPI.score >=$PPI_score_cutoff));
					
			my $sth = $dbh->prepare( $stmt );
			my $rv = $sth->execute() or die $DBI::errstr;
			if($rv < 0){print $DBI::errstr;	}
			$c=0;
			while(my @row = $sth->fetchrow_array()) {
				$g->add_edge($row[0],$row[1]);
				$ppi_status="adding nodes ".sprintf("%d/%d",$c++,$row_c); Tkx::update();	
			}
			$ppi_status="adding nodes $c/$row_c\n"; Tkx::update();
			
			$sth->finish();
		##Calculating nodes		
			$string_srch_prg=40;Tkx::update();
			$ppi_status="vertices calculation"; Tkx::update();
			
			my @V = $g->vertices;		## nodes in array
			my $V = scalar @V;			##$g->vertices;		## network sizr (n)
			#$ppi_status="diameter calculation"; Tkx::update();
			#print STDERR "$ppi_status\n";
			#my $gd = $g->diameter;		#diamemter calculation taking a lot of time; so dropped
			$ppi_status="Exporting nodes";Tkx::update();
			open (K,">all_node.txt") or die "$! all_node.txt";
			$"="\n";
			print K "@V\n";
			close K;
			$"=" ";
			sleep(2);
			$ppi_status="shortest path calculation"; Tkx::update();
			print STDERR "$ppi_status in progress. Please wait....\t";
			my $apsp = $g->APSP_Floyd_Warshall();
			Tkx::update();print STDERR "DONE\n";
			#$prg_grph_cc->g_destroy();
			#$mw->configure(-cursor=>"arrow");
			#Return the all-pairs shortest path object computed from the graph using Floyd-Warshall's algorithm. The length of a path between two vertices is the sum of weight attribute of the edges along the shortest path between the two vertices. If no weight attribute name is specified explicitly the attribute weight is assumed.
			$ppi_status=""; Tkx::update();
			$string_srch_prg=50;Tkx::update();
			my %centrality_scores_per_node;
			my %total_normalized_score_per_node;		##in a separate one as can be sorted easily
			my($tot_radiality,$tot_Ec,$tot_closeness)=( 0, 0, 0);
			my($min_radiality,$min_Ec,$min_closeness,$max_radiality,$max_Ec,$max_closeness)=( 0, 0, 0,0,0,0);
			my %node_componenet_index_tograph;			## keeps componenet idex as key and created undirected graph object as value;later updated to diameter of the object;just to save time at the cost of memory
			
##Cytohubba
			$c=0;
			$ppi_status="Centrality calculation"; Tkx::update();
			foreach my $node(@V)			#(@V)
				{					
					$centrality_scores_per_node{$node}={
								radiality=>0,
								eccentricity=>0,
								closeness=>0,
								degree=>0,
								normalized_sum=>0,
								};
					
					my($radiality,$Ec,$closeness)=(0,0,0);		##centrality measires radiality, eccenticity, closeness
					my $d = $g->degree($node); 
					my $shortest_path_to_all=0; 
					my $longest_shortest_path=-1;		## longest distance node(v) to all other in componenet C(v)
					#my $longest_shortest_path_any_two_node=-1;	## longest distance any two nodes in componenet C(v)
					foreach my $w(@V)									
					{
						
						next  if $w eq $node;
						my $l = $apsp->path_length($node, $w);
						my $a;
						eval{$a=(1/$l)};
						if($@){$a=0};
						$shortest_path_to_all+=	$a; #;
					}
					$closeness=$shortest_path_to_all;
						my $i_cc = $g->connected_component_by_vertex($node);
						my @z = $g->connected_component_by_index($i_cc);
						my $factor= ($#z+1)/$V;
					foreach my $w(@z)									
					{
						next  if $w eq $node;
						my $l = $apsp->path_length($node, $w);
						$longest_shortest_path=($l>$longest_shortest_path?$l:$longest_shortest_path);
					}
					#print "$node compnenet index=$i_cc\n";
					
					foreach my $w(@z)									
					{
						foreach my $k(@z)									
						{
							next if ($k eq $w);
							my $l = $apsp->path_length($k, $w);
							#$longest_shortest_path_any_two_node=($l>$longest_shortest_path_any_two_node?$l:$longest_shortest_path_any_two_node);
							
						}	
					}
					if(!(exists $node_componenet_index_tograph{$i_cc})){
						open(S,">componenet.txt") or die "$!";
						foreach my $w(@z)									
						{
							foreach my $k(@z)									
							{
								last if $k eq $w;
								print S "$k,$w\n";
							}
						}	
						close S;
						open(S,">import.sql") or die "$!";	
						print S "
DROP TABLE IF EXISTS tmp2;
CREATE TABLE tmp2 (proteinA VARCHAR(20) NOT NULL, proteinB VARCHAR(20) NOT NULL);
.separator ,
.import \"componenet.txt\" tmp2";
						close S;
						system("sqlite3.exe $database < import.sql");
						$node_componenet_index_tograph{$i_cc} = Graph::Undirected->new();	##storing graph for that index
						
						my $stmt = qq(select PPI.proteinA, PPI.proteinB from PPI INNER JOIN tmp2 where ((tmp2.proteinA=PPI.proteinA) AND (tmp2.proteinB=PPI.proteinB)) AND PPI.score >=$PPI_score_cutoff);	##LIMIT 1500
						
						my $sth = $dbh->prepare( $stmt );
						my $rv = $sth->execute() or die $DBI::errstr;
						if($rv < 0){
								print $DBI::errstr;
							}
						while(my @row = $sth->fetchrow_array()) {
							$node_componenet_index_tograph{$i_cc}->add_edge($row[0],$row[1]);	
							#print "$row[0],$row[1])\n";	 $k++;	$s{$row[0]}=1;		$s{$row[1]}=1;		
						}
					#$node_componenet_index_tograph_apsp{$i_cc} = $node_componenet_index_tograph{$i_cc}->APSP_Floyd_Warshall();	
					$node_componenet_index_tograph{$i_cc}=$node_componenet_index_tograph{$i_cc}->diameter;
					
					}
					else{}
					my $no_of_nodes_in_comonenet= scalar @z;
					#print "Radiality:$node\t $no_of_nodes_in_comonenet\n";
					foreach my $w(@z)									
					{
						next  if $w eq $node;
						my $l = $g->path_length($node, $w);
						$radiality+=(($node_componenet_index_tograph{$i_cc} + 1)-$l);
						#print "\t$node_componenet_index_tograph{$i_cc} \t $l\n";
						
						#$longest_shortest_path_in_componenet=($l>$longest_shortest_path_in_componenet?$l:$longest_shortest_path_in_componenet);
					}
					
					
					$Ec=$factor*(1/$longest_shortest_path);
					#$radiality=$factor*($radiality/$longest_shortest_path);
					$radiality=$factor*($radiality/($no_of_nodes_in_comonenet-1));
										
					#my $ve = $g->vertex_eccentricity($node);
					$min_radiality=($min_radiality>$radiality?$radiality:$min_radiality);
					$min_Ec=($min_Ec>$Ec?$Ec:$min_Ec);
					$min_closeness=($min_closeness>$closeness?$closeness:$min_closeness);
					$max_radiality=($max_radiality<$radiality?$radiality:$max_radiality);
					$max_Ec=($max_Ec<$Ec?$Ec:$max_Ec);
					$max_closeness=($max_closeness<$closeness?$closeness:$max_closeness);
					
					
					$tot_radiality+=$radiality;
					$tot_Ec+=$Ec;
					$tot_closeness+=$closeness;
					
					$centrality_scores_per_node{$node}={
								radiality=>$radiality,
								eccentricity=>$Ec,
								closeness=>$closeness,
								degree=>$d,
								normalized_sum=>0,
								};
								
				#printf "%s\t%d\t%f\t%f\t%f\t%d\t%d\n", $node,$centrality_scores_per_node{$node}->{degree},$centrality_scores_per_node{$node}->{radiality},$centrality_scores_per_node{$node}->{eccentricity},$centrality_scores_per_node{$node}->{closeness},$no_of_nodes_in_comonenet,$longest_shortest_path;				
			Tkx::update();
			$ppi_status="Centrality measure ".sprintf("%d/%d",$c++,$V); Tkx::update();
			}
			$dbh->disconnect();
			$string_srch_prg=60;Tkx::update();
			#normalize hub data, sort and mapp
			foreach my $n(keys %centrality_scores_per_node)		##$n=PPI id
			{
			my $tot_relative_centrality_score=0;
		
				$tot_relative_centrality_score=(($centrality_scores_per_node{$n}->{radiality}-$min_radiality)/($max_radiality-$min_radiality))+(($centrality_scores_per_node{$n}->{eccentricity}-$min_Ec)/($max_Ec-$min_Ec))+(($centrality_scores_per_node{$n}->{closeness}-$min_closeness)/($max_closeness-$min_closeness));
				
			#my $tot_relative_centrality_score=(($centrality_scores_per_node{$n}->{radiality}/$tot_radiality)+($centrality_scores_per_node{$n}->{eccentricity}/$tot_Ec)+($centrality_scores_per_node{$n}->{closeness}/$tot_closeness));
			$centrality_scores_per_node{$n}->{normalized_sum}=$tot_relative_centrality_score;
			}
			#my $=read_fasta_sequence("$L_root_path/accepted_seq_step-3.fasta"); 
			#my $o=scalar keys %centrality_scores_per_node;
			#print STDERR "$#V nodes\ntotal in centrality_scores_per_node = $o\n";
			#
			
			$string_srch_prg=70;Tkx::update();
			my %non_host_seq_with_centrality_score;	##seq_id as key
			$ppi_status="Normalization"; Tkx::update();
			foreach my $i(@$not_host_like_proteins)
			{
				my $ppi_id=$seq_id_to_PPI_id_map{$i};	#$non_host_seq_with_centrality_score{$i}=0;
				if(exists $centrality_scores_per_node{$ppi_id}){ $non_host_seq_with_centrality_score{$i}=$centrality_scores_per_node{$ppi_id}->{normalized_sum}}
				else{$non_host_seq_with_centrality_score{$i}=0;}				
			}
			#print "non host prot:".scalar @$not_host_like_proteins;
			#print "\nnon_host_seq_with_centrality_score".scalar keys %non_host_seq_with_centrality_score;
			my @sorted_filterd_ids;	##stores non_host_ids after mapping with PPI ids
			my %total_summary;		##Holds all calculation in tabular hash
			foreach my $u(sort {$non_host_seq_with_centrality_score{$b} <=> $non_host_seq_with_centrality_score{$a}} keys %non_host_seq_with_centrality_score){ 
			push @sorted_filterd_ids, $u;
			}
			$ppi_status="Exporting results"; Tkx::update();
			#print "\nsorted_filterd_ids:$#sorted_filterd_ids\n";
			open (N,">$L_root_path/Centrality_measures.txt") or die "$!";			
			print N "#No of nodes: $#V +1\n#Filtered ids: $#sorted_filterd_ids +1\n";
			#print N "$filt\t$#sorted_filterd_ids\n";
			print N  "#Sequence_id\tPPI_database_id\tdegree\tradiality\tcloseness\teccentricity\tnormalized_total_score\n";
			my (@ids_with_CC,@filtered_id,@not_filtered_id);	##@ids_with_CC keeping those with Centrality measures>0; used to calc top#% hub genes	
			foreach my $i(0..scalar @sorted_filterd_ids-1)
			{	my $r=$seq_id_to_PPI_id_map{$sorted_filterd_ids[$i]};
				if($centrality_scores_per_node{$r}->{normalized_sum} ){
				push @ids_with_CC,$sorted_filterd_ids[$i]; 					
				#$total_summary{$sorted_filterd_ids[$i]}->{filt=>1}; 
				print N "$sorted_filterd_ids[$i]\t$r\t$centrality_scores_per_node{$r}->{degree}\t$centrality_scores_per_node{$r}->{radiality}\t$centrality_scores_per_node{$r}->{closeness}\t$centrality_scores_per_node{$r}->{eccentricity}\t$centrality_scores_per_node{$r}->{normalized_sum}\n";	#$min_radiality,$min_Ec,$min_closeness,$max_radiality,$max_Ec,$max_closeness\n"; 
				}
				else{push @not_filtered_id,$sorted_filterd_ids[$i]; $total_summary{$sorted_filterd_ids[$i]}->{filt=>0} }			
			}
			$string_srch_prg=80;Tkx::update();						
			close N;
			
#=cut		### IMPORTANT:  commnet the following one lines of code while DELETING cut	;just hacking to add shotcut @filtered_id;
			
			#my @sorted_filterd_ids = @{ids_in_fasta_seq("$L_root_path/accepted_seq_step-3.fasta")}; ## remove category
			#my (@ids_with_CC,@filtered_id,@not_filtered_id);		## remove category
			#@ids_with_CC = keys(%{process_centrality_measure_file("$L_root_path/Centrality_measures.txt")}); ## remove category
			
			
			my $filt=int (($top_hub_perc/100)*($#ids_with_CC+1));				
			foreach my $i(0..scalar @sorted_filterd_ids-1)
			{
				if($i<=$filt){push @filtered_id,$sorted_filterd_ids[$i]; }		##Keep only top 20% of ids; other discarded
				else{push @not_filtered_id,$sorted_filterd_ids[$i];}
			}
			
			my $drug_like_seq_ids=fetch_seq_by_id($all_t_seq,\@filtered_id);
				write_fasta_seq($drug_like_seq_ids,"$L_root_path/accepted_seq_step-4_1.fasta");
				my $r=fetch_seq_by_id($all_t_seq,\@not_filtered_id);	##
				write_fasta_seq($r,"$L_root_path/excluded_seq_step-4_1.fasta");
				$sequence_summary{-putative_drug_targets}=scalar @filtered_id;
			$string_srch_prg=90;Tkx::update();
			$entry_drug_cand_ppi_app->delete(0, "end"); $entry_drug_cand_ppi_app->insert(0, $sequence_summary{-putative_drug_targets}) ;

		##Export SIF file cytoscape visualization
			open (SIF, ">$L_root_path/filtered_hub_genes.SIF") or die "$!";
			foreach my $seq_id(@filtered_id){
				my $i = $g->connected_component_by_vertex($seq_id_to_PPI_id_map{$seq_id});	#return an index identifying the connected component the vertex belongs to
				my @Connected_v = $g->connected_component_by_index($i); #For an undirected graph, return the vertices of the ith connected component, the indexing starting from zero. 
				foreach (@Connected_v ) {print SIF "$seq_id\tpp\t$PPI_id_to_seq_id_map{$_}\n";}
			}
			close SIF;

			sleep(1);	
			$string_srch_prg=100;Tkx::update();
			$save_result->configure(-state=>"normal");				
			
			Tkx::tk___messageBox(-message => "Run complete.\n \nSave results and parameter options.", -type=>"ok", -title=>"Success"); 	

	});
		

###do_ess_pro_blast		
	$do_ess_pro_blast->configure(-command=>sub{
		if(!$Eproteome_file){
			Tkx::tk___messageBox(-type => "ok",
					-message => "Kindly import essential proteome\n",
					-icon => "error", -title => "Import essential proteome");
			my $Eproteome=Tkx::tk___getOpenFile();
			if(!((-e "$Eproteome\.phr")and(-e "$Eproteome\.psq")and(-e "$Eproteome\.pin"))){
					Tkx::tk___messageBox(-type => "ok",
					-message => "BLAST database (*.phr,*.psq, *.pin) not found for $Eproteome_file. Use formatdb/makeblastdb to create it ",
					-icon => "error", -title => "Input file Error");
				undef $Eproteome	;	
				}
			else{$Eproteome_file=$Eproteome	;}	
			
			}
			while(!$Eproteome_file){ Tkx::update();}
			$do_PPI_search->configure(-state=>"disabled");
			$do_ess_pro_blast->configure(-state=>"disabled");
			$Eproteome_file=~ s{/}{\\}g;   $Eproteome_file='"\"'.$Eproteome_file.'"\"';

			unlink "$L_root_path/essential_protein_blast2.out.txt";
			my($gap_open, $gap_extns)=split /,/,$gap_score_2;
			my $blast2="";
			$blast2="blastall.exe -p blastp -d $Eproteome_file -i $root_path\\accepted_seq_step-3.fasta -e $e_val_2 -m $out_fmt_2 -W $word_size_2 -M $sub_matrix_2 -G $gap_open -E $gap_extns -o $root_path\\essential_protein_blast2.out -f $threshold_2 -a $use_cores"." $extra_params_BLAST2" if $blast_version eq 'old_blastall';
			$blast2="blastp.exe -db $Eproteome_file -query $root_path\\accepted_seq_step-3.fasta -evalue $e_val_2 -outfmt $out_fmt_2 -word_size $word_size_2 -matrix $sub_matrix_2 -gapopen $gap_open -gapextend $gap_extns -out $root_path\\essential_protein_blast2.out -threshold $threshold_2 -num_threads $use_cores"." $extra_params_BLAST2" if $blast_version eq 'blast+';
				`echo echo off > batch.bat`;
				`echo color B0 >> batch.bat`;
				`echo cls >> batch.bat`;				
				`echo echo :: External program BLAST :: >> batch.bat`;
				`echo echo ---------------------------------------------- >> batch.bat`;
				
				`echo echo Parameters >> batch.bat`;
				`echo echo 	Program: blastp >> batch.bat`;
				`echo echo 	Query		: $root_path\\accepted_seq_step-3.fasta >> batch.bat`;
				`echo echo 	Database	: $Eproteome_file >> batch.bat`;
				`echo echo 	E-value	: $e_val_2 >> batch.bat`;
				`echo echo 	Scoring matrix: $sub_matrix_2 >> batch.bat`;
				`echo echo 	Gap-penalty (Open,Extension): $gap_open,$gap_extns >> batch.bat`;
				`echo echo 	Word size: $word_size_2 >> batch.bat`;
				`echo echo 	Threshold for extending hits: $threshold_2 >> batch.bat`;
				`echo echo 	CPUs: $use_cores >> batch.bat`;
				
				`echo echo 	Identity percentage: $perc_identity_2 >> batch.bat`;
				`echo echo Please wait.......... >> batch.bat`;
				`echo $blast2 >> batch.bat`;
				`echo rename $root_path\\essential_protein_blast2.out essential_protein_blast2.out.txt >> batch.bat`;	##mv works
				`echo exit >> batch.bat`;
				if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
				else{system("start batch.bat ");}

				while(!(-e "$L_root_path/essential_protein_blast2.out.txt")){
				#print "$blast_prg2  $sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-orthologous})--\n";
				my $t=$sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-orthologous}+$sequence_summary{-host_orthologs});
				$blast_prg2=blast_progress("$root_path\\essential_protein_blast2.out","$L_root_path/essential_protein_blast2.out",$t); 
				sleep(3);
				Tkx::update(); 
				
				}	##wait till essential_protein_blast2.out.txt is available
				#$front_page_status
				$blast_prg2=100;Tkx::update();
				if(count_blast_hits("$L_root_path/essential_protein_blast2.out.txt")<1){
				Tkx::tk___messageBox(-message => "No BLAST hit found!!!Aborting analysis\nPlease change BLAST parameters and run analysis again.", -type=>"ok", -title=>"Alert",-icon=>'warning'); 
				return ();
				}
				
				my($essential_drug_proteins,$not_drug_like_proteins)=process_host_blast_out("$L_root_path/accepted_seq_step-3.fasta","$L_root_path/essential_protein_blast2.out.txt" , $perc_identity_2);
				$sequence_summary{-drug_target_non_homologs}=scalar @$not_drug_like_proteins;
				$sequence_summary{-drug_target_homologs}=scalar @$essential_drug_proteins;
				my $drug_like_seq_id=fetch_seq_by_id($all_t_seq,$essential_drug_proteins);
				write_fasta_seq($drug_like_seq_id,"$L_root_path/accepted_seq_step-4.fasta");
				my $r=fetch_seq_by_id($all_t_seq,$not_drug_like_proteins);	##
				write_fasta_seq($r,"$L_root_path/excluded_seq_step-4.fasta");
				
				my ($most_drug_like_seq_id,$less_drug_like_seq_id)=filter_by_bit_score("$L_root_path/essential_protein_blast2.out.txt",100);
				$sequence_summary{-less_drug_like_seq_id}=scalar @$less_drug_like_seq_id; ##<100 bit
				$sequence_summary{-putative_drug_targets}=scalar @$most_drug_like_seq_id; ##>100 bit
				my $drug_like_seq_ids=fetch_seq_by_id($all_t_seq,$most_drug_like_seq_id);
				write_fasta_seq($drug_like_seq_ids,"$L_root_path/accepted_seq_step-4_1.fasta");
				my $r=fetch_seq_by_id($all_t_seq,$less_drug_like_seq_id);	##
				write_fasta_seq($r,"$L_root_path/excluded_seq_step-4_1.fasta");
				my $rejected=$sequence_summary{-drug_target_non_homologs}+$sequence_summary{-less_drug_like_seq_id};
				#$entry_non_ess_seq->delete(0, "end"); $entry_non_ess_seq->insert(0, $rejected) ;
				$entry_drug_cand_seq_app->delete(0, "end"); $entry_drug_cand_seq_app->insert(0, $sequence_summary{-putative_drug_targets}) ;
				sleep(1);
				
				$save_result->configure(-state=>"normal");								
				Tkx::tk___messageBox(-message => "Run complete.\n \nSave results and parameter options.", -type=>"ok", -title=>"Success"); 	
	});		
	$do_PPI_search->configure(-state=>"normal");
	$do_ess_pro_blast->configure(-state=>"normal");
	Tkx::update();
	###Enable downstream analysis
		$dwn_str_anal->entryconfigure("Broadspectrum analysis",-state=>"normal",-command =>sub{broad_spect_run($frm,$run_button)});
		$dwn_str_anal->entryconfigure("GO analysis",-state=>"normal");
		$dwn_str_anal->entryconfigure("Sub-cellular localization",-state=>"normal");
		$save_options_but->configure(-state=>"normal");
}




##Menu->Downstream analaysis-> Broad-spectrum analysis 
sub broad_spect_run
{
	my $frm=shift;
	my $run_button = shift;
	
	my $frm_top=$$frm->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 500,-padding => "0 0 0 0");
	$frm_top->g_grid(-column=>0,-row=>0,-sticky=>"nswe");
		
	Tkx::image_create_photo( "BANER", -file => "banner.gif");
	($frm_top->new_ttk__label(-image=>'BANER'))->g_grid(-column=>0, -row=>0,-sticky=>"nwes",-columnspan=>2, -padx=>60);	
	my $heading = $frm_top->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool",-justify=>"center",-foreground=>"blue",-font => "Helvetica 16 bold underline");
	$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>50);
	
	$frm_top->new_ttk__label(-text=>"BROAD-SPECTRUM TARGET ANALYSIS",-justify=>"left",-foreground=>"darkgreen",-font => "Helvetica 12")->g_grid(-column=>0,-row=>1,-sticky=>"s",-padx=>0);
	
	my ($input_seq,$brd_blast_db); 		##
	$input_seq = ($Tproteome_file?"$L_root_path/accepted_seq_step-4_1.fasta":"");
	my $new_frm = $frm_top->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 200,-padding => "10 20 0 0");
	$new_frm->g_grid(-column=>0,-row=>2,-sticky=>"nswe");
	
	$new_frm->new_ttk__label(-text=>"Input target sequences")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$new_frm->new_ttk__label(-text=>"Output folder")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	$new_frm->new_ttk__label(-text=>"Use Settings (Alt+s) Menu for parameter setting",-foreground=>"red",-justify=>"right")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	
	$new_frm ->new_ttk__entry(-textvariable => \$input_seq,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-columnspan=>2);
	$new_frm ->new_ttk__entry(-textvariable => \$root_path,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$input_seq = Tkx::tk___getOpenFile(-parent=>$mw);$mw->g_raise();  		## Assuming that no project defined;
	if($input_seq){
		$sequence_summary{-putative_drug_targets}=count_fasta_seq($input_seq); 	## reset drug target counts
		$input_seq=~ s{/}{\\}g; $input_seq='"'.$input_seq.'"'; 					##Convert to windows format
		}	
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$root_path = Tkx::tk___chooseDirectory(-parent=>$mw);$mw->g_raise();
	
	if(!$root_path){$root_path=win_path($L_root_path);}  							##if cancel pressed ; fail safe
	else{
		$root_path =~s/"//g;														##incase of Default Desktop path
		$L_root_path=$root_path;													##preserve UNIX format
		$root_path=win_path($root_path); 						##Convert to windows format	
	}
	
	})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1,-sticky=>"wn");	
	
	my $canvas = $new_frm->new_tk__canvas(-scrollregion => "0 0 1000 1000",-width=>400, -height=>200);
	$canvas->g_grid(-column=>1, -row=>5, -padx=>5, -pady=>5,-sticky=>"nw",-columnspan=>4);
	my $hscroll = $new_frm->new_tk__scrollbar(-orient => "horizontal", -command => [$canvas, "xview"]);
	my $vscroll = $new_frm->new_tk__scrollbar(-orient => "vertical", -command => [$canvas, "yview"]);
	$hscroll->g_grid(-column => 1, -row => 6, -sticky => "we",-columnspan=>4);
	$vscroll->g_grid(-column => 5, -row => 5, -sticky => "ns");
	$new_frm->new_ttk__sizegrip()->g_grid(-column => 5, -row => 6, -sticky => "se");
	$canvas->configure(-yscrollcommand => [$vscroll, "set"], -xscrollcommand => [$hscroll, "set"]);
	
	$new_frm->new_ttk__label(-text=>"Broad-spectrum analysis:")->g_grid(-column=>0,-row=>7,-padx=>1,-pady=>1,-sticky=>"w");
	my $prg_blast_brd_sp=$new_frm->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$blast_prg3);
	$prg_blast_brd_sp->g_grid(-column=>1,-row=>7,-padx=>3,-pady=>1,-sticky=>"w");
	$new_frm->new_ttk__label(-text=>"Minimum no. of conserved species:")->g_grid(-column=>0,-row=>8,-padx=>1,-pady=>1,-sticky=>"w");
	$new_frm->new_ttk__entry(-width=>5,-textvariable=>\$broad_spe_species_per_query)->g_grid(-column=>1,-row=>8,-padx=>1,-pady=>1,-sticky=>"w");
	
	my $filter_broad_spe_but=$new_frm->new_button(-text=>"Apply filter",-width=>10, -state=>"disabled");
	$filter_broad_spe_but->g_grid(-column=>2, -row=>8,-padx=>2,-pady=>1,-sticky=>"w");
	my $save_result_but=$new_frm->new_button(-text=>"Save results",-width=>10, -state=>"disabled");
	$save_result_but->g_grid(-column=>3, -row=>8,-padx=>2,-pady=>2,-sticky=>"e");
	

	$$run_button->configure(-state=>"normal", -command =>sub
	{
	$$run_button->configure(-state=>"disabled");
	unlink "$L_root_path/broad_spe_blast3.out.txt" if (-e "$L_root_path/broad_spe_blast3.out.txt");
	my($gap_open, $gap_extns)=split /,/,$gap_score_3;
	my @broad_spectrum_pathogen_db_list = split /,/,$broad_spectrum_pathogen_db_list;
	if(scalar@broad_spectrum_pathogen_db_list <1){ Tkx::tk___messageBox(-message => "ERROR: no database selected. This coould be due to none of the pathogenes selected to be compared or no pathogen protome has been added.\n Go to  Setting->Downstream analysis, for selecting pathogenes.\nUse Utility menu to add pathogen proteomes ", -type=>"ok", -title=>"Alert", -icon=>"warning" ); return();}
	
	`echo echo off > batch.bat`;
	`echo color 90 >> batch.bat`;
	`echo 0 > brd_run.txt`;
	open (BAT, ">> batch.bat") or die "$! batch.bat";
	my $org=1;
	foreach my $pathogen_db(@broad_spectrum_pathogen_db_list){
			$org;
			my $tot=$#broad_spectrum_pathogen_db_list+1;
			my $blast3="";
			$blast3="blastall.exe -p blastp -d $pathogen_db -i $input_seq -e $e_val_3 -m $out_fmt_3 -W $word_size_3 -M $sub_matrix_3 -G $gap_open -E $gap_extns -o $root_path\\broad_spe_blast3.out1 -a $use_cores -f $threshold_3"." $extra_params_BLAST3" if ($blast_version eq 'old_blastall');
			$blast3="blastp -db $pathogen_db -query $input_seq -evalue $e_val_3 -outfmt $out_fmt_3 -word_size $word_size_3 -matrix $sub_matrix_3 -gapopen $gap_open -gapextend $gap_extns -out $root_path\\broad_spe_blast3.out1 -num_threads $use_cores -threshold $threshold_3"." $extra_params_BLAST3" if ($blast_version eq 'blast+');
			
				print BAT "cls\n";
				print BAT "echo :: External program BLAST :: \n";
				print BAT "echo ----------------------------------------------\n";
				print BAT "echo $org > brd_run.txt\n";
				print BAT "echo Parameters\n";
				print BAT "echo 	Program: blastp\n";
				print BAT "echo 	Query		: $input_seq \n";
				print BAT "echo 	Database	: $pathogen_db (see locale_dat dir.) \n";
				print BAT "echo 	E-value	: $e_val_3\n";
				print BAT "echo 	Scoring matrix: $sub_matrix_3 \n";
				print BAT "echo 	Gap-penalty (Open,Extension): $gap_open,$gap_extns \n";
				print BAT "echo 	Word size : $word_size_3 \n";
				print BAT "echo 	Threshold for extending hits : $threshold_3\n";
				print BAT "echo 	CPUs : $use_cores \n";
					
				print BAT "echo $org of $tot processing \n";
				print BAT "echo Please wait.......... \n";
				print BAT "$blast3\n";
				#`echo rename $root_path\\broad_spe_blast3.out broad_spe_blast3.out1 >> batch.bat`;	##mv works
				if($org>1){ 
				print BAT "type $root_path\\broad_spe_blast3.out1 >> $root_path\\broad_spe_blast3.out.f \n";
				}
				else {print BAT "type $root_path\\broad_spe_blast3.out1 > $root_path\\broad_spe_blast3.out.f\n";}
				
				#print BAT "echo \necho \necho \n";
				
				$org++;
														##wait till essential_protein_blast2.out.txt is available	
	}
	close BAT;
	
	`echo rename $root_path\\broad_spe_blast3.out.f broad_spe_blast3.out.txt >> batch.bat`;	##mv works
	`echo exit >> batch.bat`;
	
	if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
	else{system("start batch.bat ");}
	
	$blast_prg3=0;
	while(!(-e "$L_root_path/broad_spe_blast3.out.txt")){
		open(K,"< brd_run.txt") or die "$! brd_run.txt";		
		my $blast_prg_o=<K>;
		close K;
		chomp($blast_prg_o);
		$blast_prg3=($blast_prg_o/$#broad_spectrum_pathogen_db_list)*100;
		#blast_progress("$root_path\\broad_spe_blast3.out","$L_root_path/broad_spe_blast3.out",$sequence_summary{-putative_drug_targets}*10 ); ##as this blast is not with one hit per one query, so miscalculation happens, so multipled 10,  10 hits per query
		$blast_prg3=96 if $blast_prg3>96;		##fail safe
		sleep(3);Tkx::update(); 	
	}
	$blast_prg3=100;  Tkx::update();
	if(count_blast_hits("$L_root_path/broad_spe_blast3.out.txt")<1){
		Tkx::tk___messageBox(-message => "No BLAST hit found!!!Aborting analysis\nPlease change BLAST parameters and run analysis again.", -type=>"ok", -title=>"Alert",-icon=>'warning'); 
		return ();
	}
		my $ref_broad_spec_counts=process_broad_spe_BLAST_out("$L_root_path/broad_spe_blast3.out.txt",$perc_identity_3); 				##\%hash
		my %broad_spec_counts;					
		foreach my $a( keys %$ref_broad_spec_counts ){
			$broad_spec_counts{$a}=scalar @{$ref_broad_spec_counts->{$a}};
		}
		my (@query_id,@cons_counts,@bar_dclrs);
		my $total_hits=0;
		foreach my $t(sort { $broad_spec_counts{$b} <=> $broad_spec_counts{$a} } keys %broad_spec_counts ){
			push (@query_id,$t);
			push (@cons_counts,$broad_spec_counts{$t});
			push @bar_dclrs,'black';
		}
		
		$sequence_summary{-broad_spectrum}=scalar @query_id;	##update later on applying filter
		
		if(!(scalar @query_id)){ Tkx::tk___messageBox(-message => "None of the queries are conserved in any of the species", -type=>"ok", -title=>"Alert", -icon=>"warning" ); return();}
		
		my $tot_cons_count=0; my @cons_counts_perc; my $max_cons=-1;
		foreach(@cons_counts){ my $x=$_;$max_cons=($x>$max_cons?$x:$max_cons); $tot_cons_count+=$x;}
		#print "max : $max_cons\n";
		#foreach my $i(0.. scalar@cons_counts-1){$cons_counts_perc[$i]=($cons_counts[$i]/$tot_cons_count)*100;}
		my @data=(\@query_id,\@cons_counts);
		my($img_width,$img_height)=(1000, 600);
		my $graph = GD::Graph::bars->new(800, 400);
		$graph->set( 
			x_label           => 'Putative drug targets',
			y_label           => 'Conserved in # of species',
			title             => 'Broadspectrum analysis',
			y_max_value       => $max_cons+5,
			#y_tick_number     => 8,
			#y_label_skip      => 2,
			bar_width			=>20,
			bar_spacing			=>0,
			x_labels_vertical	=>1,
			transparent			=>0,
			accentclr		=> 'black',
			#boxclr			=> 'gray',
			axislabelclr		=>'black',
			dclrs 				=>\@bar_dclrs,
		) or die $graph->error;
		my $gd = $graph->plot(\@data,correct_width => 0) or die $graph->error;
		my $img_broad_spec;
		
			$img_broad_spec="$L_root_path/broad_spec.gif";
			open(IMG, ">$img_broad_spec") or die $!;
			binmode IMG;
			print IMG $gd->gif;
			close IMG;
		
	  Tkx::image_create_photo( "BROAD_SPE", -file => "$img_broad_spec");
	  $canvas->create_image(400, 0, -image=>"BROAD_SPE", -anchor =>'n' );
	  
	  $filter_broad_spe_but->configure(-state=>"normal");
	  my(@filter_ids,@not_filter_ids);
	  $filter_broad_spe_but->configure(-command=>sub{
			undef @filter_ids; undef @not_filter_ids;
			$save_result_but->configure(-state=>"disabled");
			#$$drug_bank_blast_but->configure(-state=>"disabled");
			foreach my $t(sort keys %broad_spec_counts ){
				if($broad_spec_counts{$t}>=$broad_spe_species_per_query){push (@filter_ids,$t);}
				else{push (@not_filter_ids,$t); 
				}
			}
			$sequence_summary{-broad_spectrum}=scalar @filter_ids;
			Tkx::tk___messageBox(-message => "$sequence_summary{-broad_spectrum} queries filtered!!!\nNow click 'Save results' to save sequences.",-type=>"ok", -title=>"Success"); 
			my $r=fetch_seq_by_id(read_fasta_sequence(unix_path($input_seq)),\@filter_ids);	##
			write_fasta_seq($r,"$L_root_path/accepted_seq_step-5.fasta");
			my $r=fetch_seq_by_id(read_fasta_sequence(unix_path($input_seq)),\@not_filter_ids);	##
			write_fasta_seq($r,"$L_root_path/excluded_seq_step-5.fasta");
			#$$drug_bank_blast_but->configure(-state=>"normal");
			$save_result_but->configure(-state=>"normal");
			
			
				$save_result_but->configure(-command=>sub{
				my $save_result = Tkx::tk___getSaveFile();
				my $r=fetch_seq_by_id(read_fasta_sequence(unix_path($input_seq)),\@filter_ids);	##
				write_fasta_seq($r,"$save_result") if $save_result;
				});
		
		});
		
	Tkx::tk___messageBox(-message => "Run complete.\nChoose minimum number of conserved species and Click on Apply button.",-type=>"ok", -title=>"Success"); 
	}); ##END RUN button
	
}	## END BROAD spect

##Menu->Downstream analaysis-> Compare with known target
sub comp_known_DT
{
	my $frm=shift;
	my $run_button = shift;	
	my $frm_top=$$frm->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 500,-padding => "0 0 0 0");
	$frm_top->g_grid(-column=>0,-row=>0,-sticky=>"nswe");
		
	Tkx::image_create_photo( "BANER", -file => "banner.gif");
	($frm_top->new_ttk__label(-image=>'BANER'))->g_grid(-column=>0, -row=>0,-sticky=>"nwes",-columnspan=>2, -padx=>60);	
	my $heading = $frm_top->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool",-justify=>"center",-foreground=>"blue",-font => "Helvetica 16 bold underline");
	$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>50);
	$frm_top->new_ttk__label(-text=>"COMPARE TARGETS WITH KNOWN DRUG-TARGETS",-justify=>"left",-foreground=>"darkgreen",-font => "Helvetica 12")->g_grid(-column=>0,-row=>1,-sticky=>"s",-padx=>0);
	
	my $input_seq;
	$input_seq = ($Tproteome_file?"$L_root_path/accepted_seq_step-4_1.fasta":"");	##Skipp if project defined;
		
	my $new_frm = $frm_top->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 200,-padding => "0 0 50 0");
	$new_frm->g_grid(-column=>0,-row=>2,-sticky=>"nswe");
	
	$new_frm->new_ttk__label(-text=>"Input target sequences")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$new_frm->new_ttk__label(-text=>"Output folder")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	$new_frm->new_ttk__label(-text=>"Use Settings (Alt+s) Menu for parameter setting",-foreground=>"red",-justify=>"right")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw", -columnspan=>2);
	
	$new_frm ->new_ttk__entry(-textvariable => \$input_seq,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-columnspan=>2);
	$new_frm ->new_ttk__entry(-textvariable => \$root_path,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
		
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$input_seq = Tkx::tk___getOpenFile(-parent=>$mw);$mw->g_raise();  		## Assuming that no project defined;
	if($input_seq){
		$sequence_summary{-putative_drug_targets}=count_fasta_seq($input_seq); 	## reset drug target counts; not called if in project call
		$input_seq=~ s{/}{\\}g; $input_seq='"'.$input_seq.'"'; 					##Convert to windows format
	  }	
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
		$root_path = Tkx::tk___chooseDirectory(-parent=>$mw);$mw->g_raise();		
		if(!$root_path){$root_path=win_path($L_root_path);}  							##if cancel pressed ; fail safe
		else{
		$root_path =~s/"//g;														##incase of Default Desktop path
		$L_root_path=$root_path;													##preserve UNIX format
		$root_path=win_path($root_path); 						##Convert to windows format	
		}		
	})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1,-sticky=>"wn");	
	
	my $canvas = $new_frm->new_tk__canvas(-scrollregion => "0 0 500 500",-width=>400, -height=>200);
	$canvas->g_grid(-column=>1, -row=>5, -padx=>5, -pady=>5,-sticky=>"nw",-columnspan=>4);

	$new_frm->new_ttk__label(-text=>"Run progress:")->g_grid(-column=>0,-row=>7,-padx=>1,-pady=>1,-sticky=>"w");
	my $prg_blast_brd_sp=$new_frm->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$blast_prg4);
	$prg_blast_brd_sp->g_grid(-column=>1,-row=>7,-padx=>3,-pady=>1,-sticky=>"w");
	
	my $Save_novel_targets_predicted_but=$new_frm->new_button(-text=>"Save novel targets",-width=>15, -state=>"disabled");
	$Save_novel_targets_predicted_but->g_grid(-column=>0, -row=>8,-padx=>2,-pady=>1,-sticky=>"w");
	my $Save_Known_targets_predicted_but=$new_frm->new_button(-text=>"Save known targets",-width=>15, -state=>"disabled");
	$Save_Known_targets_predicted_but->g_grid(-column=>1, -row=>8,-padx=>2,-pady=>2,-sticky=>"w");
	my $Save_known_targets_predicted_tsv_but=$new_frm->new_button(-text=>"Known targets annotations",-width=>20, -state=>"disabled");
	$Save_known_targets_predicted_tsv_but->g_grid(-column=>0, -row=>9,-padx=>2,-pady=>2,-sticky=>"w");
	
	
	$$run_button->configure(-state=>"normal", -command =>sub {
	
		if(scalar @$ref_drug_db_array<1){Tkx::tk___messageBox(-message => "ERROR: No Drug target database found. YOu may add a new drug target database using Utility menu",-type=>"ok", -title=>"Alert",-icon=>"warning"); return();}
			$$run_button->configure(-state=>"disabled");
			unlink "$L_root_path/drug_target_blast4.out.txt";
			my($gap_open, $gap_extns)=split /,/,$gap_score_1;
			my $blast4="";
			$blast4="blastall.exe -p blastp -d $drug_blast_db_names -i $input_seq -e $e_val_4 -m $out_fmt_4 -W $word_size_4 -M $sub_matrix_4 -G $gap_open -E $gap_extns -o $root_path\\drug_target_blast4.out -a $use_cores -f $threshold_4"." $extra_params_BLAST4" if ($blast_version eq 'old_blastall');
			$blast4="blastp -db $drug_blast_db_names -query $input_seq -evalue $e_val_4 -outfmt $out_fmt_4 -word_size $word_size_4 -matrix $sub_matrix_4 -gapopen $gap_open -gapextend $gap_extns -out $root_path\\drug_target_blast4.out -num_threads $use_cores -threshold $threshold_4"." $extra_params_BLAST4" if ($blast_version eq 'blast+');
			
				`echo echo off > batch.bat`;
				`echo color 90 >> batch.bat`;
				`echo cls >> batch.bat`;
				`echo echo :: External program BLAST :: >> batch.bat`;
				`echo echo ---------------------------------------------- >> batch.bat`;
				
				`echo echo Parameters >> batch.bat`;
				`echo echo 	Program: blastp >> batch.bat`;
				`echo echo 	Query		: $input_seq >> batch.bat`;
				`echo echo 	Database	: drug_blast_db_names (see locale_dat\ dir.) >> batch.bat`;
				`echo echo 	E-value	: $e_val_4 >> batch.bat`;
				`echo echo 	Scoring matrix: $sub_matrix_4 >> batch.bat`;
				`echo echo 	Gap-penalty (Open,Extension): $gap_open,$gap_extns >> batch.bat`;
				`echo echo 	Word size : $word_size_4 >> batch.bat`;
				`echo echo 	Threshold for extending hits : $threshold_4 >> batch.bat`;
				`echo echo 	CPUs : $use_cores >> batch.bat`;
				
				
				`echo echo Please wait.......... >> batch.bat`;
				`echo $blast4 >> batch.bat`;
				`echo rename $root_path\\drug_target_blast4.out drug_target_blast4.out.txt >> batch.bat`;	##mv works
				`echo exit >> batch.bat`;
				
				
				
				if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
				else{system("start batch.bat ");}
				
				while(!(-e "$L_root_path/drug_target_blast4.out.txt")){
				$blast_prg4=blast_progress("$root_path\\drug_target_blast4.out","$L_root_path/drug_target_blast4.out",$sequence_summary{-putative_drug_targets} ); sleep(3); Tkx::update(); 
				}	##wait till blast4.out.txt is available
				$blast_prg4=100;  Tkx::update();
				if(count_blast_hits("$L_root_path/drug_target_blast4.out.txt")<1){
					Tkx::tk___messageBox(-message => "No BLAST hit found!!!Aborting analysis\nPlease change BLAST parameters and run analysis again.", -type=>"ok", -title=>"Alert",-icon=>'warning'); 
					return ();
				}
				
				my($known_drug_targets,$novel_drug_targets)=process_host_blast_out(unix_path($input_seq),"$L_root_path/drug_target_blast4.out.txt",0);
				
				
				my @all_known_drug_targets=keys %{$drug_target_annot};
				my ($not_identified_drugTarget,$identified_drugTarget,$identified_novel_drugTarget)=(scalar@all_known_drug_targets-scalar@$known_drug_targets,scalar@$known_drug_targets,scalar@$novel_drug_targets);
				
				$canvas->create_rectangle(100, 110, 50, 10, -fill => "red",);	#x1, y1, x2, and y2
				$canvas->create_rectangle(10, 60, 100, 110, -fill => "yellow",-stipple=>"gray75");		#x1, y1, x2, and y2
				$canvas->create_text(25, 70, -text=>"$identified_novel_drugTarget");		#x1, y1,  
				$canvas->create_text(65, 90, -text=>"$identified_drugTarget");		#x1, y1,  
				$canvas->create_text(65, 20, -text=>"$not_identified_drugTarget");		#x1, y1,  

				$canvas->create_rectangle(150, 10, 160, 20, -fill => "red",);	#x1, y1, x2, and y2
				$canvas->create_rectangle(150, 30, 160, 40, -fill => "yellow",-stipple=>"gray75");		#x1, y1, x2, and y2
				$canvas->create_text(210, 15, -text=>"Known drugtarget db", );		#x1, y1,  
				$canvas->create_text(200, 35, -text=>"Identified targets",);		#x1, y1, 
				
				##venn diagram
				$Save_Known_targets_predicted_but->configure(-state=>"normal",-command=>sub{
					my $save_result = Tkx::tk___getSaveFile();
					#$crt_win->g_raise();
					my $r=fetch_seq_by_id(read_fasta_sequence(unix_path($input_seq)),$known_drug_targets);	##
					write_fasta_seq($r,"$save_result") if $save_result ;
				});
				$Save_novel_targets_predicted_but->configure(-state=>"normal",-command=>sub{
					my $save_result = Tkx::tk___getSaveFile();
					#$crt_win->g_raise();
					my $r=fetch_seq_by_id(read_fasta_sequence(unix_path($input_seq)),$novel_drug_targets);	##
					write_fasta_seq($r,"$save_result") if $save_result ;;
				});
				$Save_known_targets_predicted_tsv_but->configure(-state=>"normal",-command=>sub{
					my $save_result = Tkx::tk___getSaveFile();
					#$crt_win->g_raise();
					if ($save_result){
						open (O,">$save_result") or die "$! $save_result\n";
						my $c=0;
						print O "Sl no\tQuery id\tTaeget hit id\tDescription\tDrug bound\tUniprot accn\n";
						open (D,"<$L_root_path/drug_target_blast4.out.txt") or die "$! $L_root_path/drug_target_blast4.out.txt\n";
						while(<D>){
							chomp;
							++$c;
							my @l=split /\s+/,$_;
							##sp|Q8GBW6|12S_PROFR	Methylmalonyl-CoA carboxyltransferase 12S subunit	DB04045; DB04183	Propionibacterium freudenreichii subsp. shermanii	Q8GBW6
							print O "$l[0]\n";
							print O "$c\t$l[0]\t$l[1]\t$drug_target_annot->{$l[1]}->[1]\t$drug_target_annot->{$l[1]}->[2]\t$drug_target_annot->{$l[1]}->[3]\t$drug_target_annot->{$l[1]}->[4]\n";		
						}
						close O; close D;
					}
				});
	
	});	##END of run button sub
}

##Menu->Downstream analaysis-> GO enrichment analysis
#Args##
#Returns
#http://www.uniprot.org/uniprot/?query=gene%3aseca+AND+organism%3a83333&format=tab&columns=id,protein%20names,genes,go,go-id
#http://geneontology.org/page/download-ontology
	##Extract go_daily-termdb-tables.gz and copy the term.txt
	
##E coli reference proteome as standard
			#ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625_83333.fasta.gz

sub GO_analysis{
	
	my $frm=shift;
	my $run_button = shift;	
	my $frm_top=$$frm->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 500,-padding => "0 0 0 0");
	$frm_top->g_grid(-column=>0,-row=>0,-sticky=>"nswe");
		
	Tkx::image_create_photo( "BANER", -file => "banner.gif");
	($frm_top->new_ttk__label(-image=>'BANER'))->g_grid(-column=>0, -row=>0,-sticky=>"nwes",-columnspan=>2, -padx=>60);	
	my $heading = $frm_top->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool",-justify=>"center",-foreground=>"blue",-font => "Helvetica 16 bold underline");
	$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>50);
	$frm_top->new_ttk__label(-text=>"GENE ONTOLOGY ENRICHMENT ANALYSIS",-justify=>"left",-foreground=>"darkgreen",-font => "Helvetica 12")->g_grid(-column=>0,-row=>1,-sticky=>"s",-padx=>0);
	
	my ($input_seq,$background_seq);
	$input_seq = ($Tproteome_file?"$L_root_path/accepted_seq_step-4_1.fasta":"");	##Skipp if project defined;
	$background_seq = ($Tproteome_file?$Tproteome_file:"");							##Skipp if project defined;
		
	my $new_frm = $frm_top->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 200,-padding => "0 0 50 0");
	$new_frm->g_grid(-column=>0,-row=>2,-sticky=>"nswe");
	
	$new_frm->new_ttk__label(-text=>"Input background sequences")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");	
	$new_frm->new_ttk__label(-text=>"Input target sequences")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	$new_frm->new_ttk__label(-text=>"Output folder")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$new_frm->new_ttk__label(-text=>"Use Settings (Alt+s) Menu for parameter setting",-foreground=>"red",-justify=>"right")->g_grid(-column=>0,-row=>3,-padx=>2,-pady=>0,-sticky=>"nw", -columnspan=>2);
	
	$new_frm ->new_ttk__entry(-textvariable => \$background_seq,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-columnspan=>2);
	$new_frm ->new_ttk__entry(-textvariable => \$input_seq,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	$new_frm ->new_ttk__entry(-textvariable => \$root_path,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2);
		
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$background_seq = Tkx::tk___getOpenFile(-parent=>$mw);$mw->g_raise();  		## Assuming that no project defined;
	if($background_seq){
		$sequence_summary{-total_seq}=count_fasta_seq($background_seq); 					## reset total background_seq counts; not called if in project call
		$background_seq=~ s{/}{\\}g; $background_seq='"'.$background_seq.'"'; 					##Convert to windows format
	}	
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$input_seq = Tkx::tk___getOpenFile(-parent=>$mw);$mw->g_raise();  		## Assuming that no project defined;
	if($input_seq){
		$sequence_summary{-putative_drug_targets}=count_fasta_seq($input_seq); 	## reset drug target counts; not called if in project call
		$input_seq=~ s{/}{\\}g; $input_seq='"'.$input_seq.'"'; 					##Convert to windows format
	  }	
	})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
		$root_path = Tkx::tk___chooseDirectory(-parent=>$mw);$mw->g_raise();		
		if(!$root_path){$root_path=win_path($L_root_path);}  							##if cancel pressed ; fail safe
		else{
			$root_path =~s/"//g;														##incase of Default Desktop path
			$L_root_path=$root_path;													##preserve UNIX format
			$root_path=win_path($root_path); 						##Convert to windows format				
			}		
	})->g_grid(-column=>4,-row=>2,-padx=>2,-pady=>1,-sticky=>"wn");	
	
	my $notebook = $new_frm->new_ttk__notebook;
	$notebook->g_grid(-column=>0, -row=>5, -padx=>5, -pady=>0,-sticky=>"nw",-columnspan=>4, -rowspan=>5,);
	my $f1 = $notebook->new_ttk__frame; # first page, which would get widgets gridded into it
	my $f2 = $notebook->new_ttk__frame; # second page
	my $f3 = $notebook->new_ttk__frame; # Third page
	$notebook->add($f1, -text => "Biological Process");
	$notebook->add($f2, -text => "Molecular Function");
	$notebook->add($f3, -text => "Cellular Component");

	my $canvas1 = $f1->new_tk__canvas(-scrollregion => "0 0 1000 1000",-width=>400, -height=>250);
	$canvas1->g_grid(-column=>1, -row=>5, -padx=>5, -pady=>1,-sticky=>"nw",-columnspan=>4);
	my $hscroll1 = $f1->new_tk__scrollbar(-orient => "horizontal", -command => [$canvas1, "xview"]);
	my $vscroll1 = $f1->new_tk__scrollbar(-orient => "vertical", -command => [$canvas1, "yview"]);
	$hscroll1->g_grid(-column => 1, -row => 6, -sticky => "we",-columnspan=>4);
	$vscroll1->g_grid(-column => 5, -row => 5, -sticky => "ns");
	$f1->new_ttk__sizegrip()->g_grid(-column => 5, -row => 6, -sticky => "se");
	$canvas1->configure(-yscrollcommand => [$vscroll1, "set"], -xscrollcommand => [$hscroll1, "set"]);
	
	my $canvas2 = $f2->new_tk__canvas(-scrollregion => "0 0 1000 1000",-width=>400, -height=>250);
	$canvas2->g_grid(-column=>1, -row=>5, -padx=>5, -pady=>1,-sticky=>"nw",-columnspan=>4);
	my $hscroll2 = $f2->new_tk__scrollbar(-orient => "horizontal", -command => [$canvas2, "xview"]);
	my $vscroll2 = $f2->new_tk__scrollbar(-orient => "vertical", -command => [$canvas2, "yview"]);
	$hscroll2->g_grid(-column => 1, -row => 6, -sticky => "we",-columnspan=>4);
	$vscroll2->g_grid(-column => 5, -row => 5, -sticky => "ns");
	$f2->new_ttk__sizegrip()->g_grid(-column => 5, -row => 6, -sticky => "se");
	$canvas2->configure(-yscrollcommand => [$vscroll2, "set"], -xscrollcommand => [$hscroll2, "set"]);
	
	my $canvas3 = $f3->new_tk__canvas(-scrollregion => "0 0 1000 1000",-width=>400, -height=>250);
	$canvas3->g_grid(-column=>1, -row=>5, -padx=>5, -pady=>1,-sticky=>"nw",-columnspan=>4);
	my $hscroll3 = $f3->new_tk__scrollbar(-orient => "horizontal", -command => [$canvas3, "xview"]);
	my $vscroll3 = $f3->new_tk__scrollbar(-orient => "vertical", -command => [$canvas3, "yview"]);
	$hscroll3->g_grid(-column => 1, -row => 6, -sticky => "we",-columnspan=>4);
	$vscroll3->g_grid(-column => 5, -row => 5, -sticky => "ns");
	$f3->new_ttk__sizegrip()->g_grid(-column => 5, -row => 6, -sticky => "se");
	$canvas3->configure(-yscrollcommand => [$vscroll3, "set"], -xscrollcommand => [$hscroll3, "set"]);

	$new_frm->new_ttk__label(-text=>"Run progress (background):")->g_grid(-column=>5,-row=>5,-padx=>1,-pady=>5,-sticky=>"w");
	my $prg_blast_brd_sp_1=$new_frm->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$blast_prg5);
	$prg_blast_brd_sp_1->g_grid(-column=>6,-row=>5,-padx=>3,-pady=>1,-sticky=>"w");
	
	my $GO_ana_progress=0;
	$new_frm->new_ttk__label(-text=>"Calculating GO :")->g_grid(-column=>5,-row=>6,-padx=>5,-pady=>1,-sticky=>"e");
	my $prg_blast_brd_sp_2=$new_frm->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate',-variable=>\$GO_ana_progress);
	$prg_blast_brd_sp_2->g_grid(-column=>6,-row=>6,-padx=>3,-pady=>1,-sticky=>"e");
	
	my $seq_used_for_analysis=$new_frm->new_button(-text=>"Sequeces with GO",-width=>15, -state=>"disabled");
	$seq_used_for_analysis->g_grid(-column=>5, -row=>9,-padx=>2,-pady=>5,-sticky=>"w");
	
	#my $save_go_enrich_table=$new_frm->new_button(-text=>"Export result as text",-width=>15, -state=>"disabled");
	#$save_go_enrich_table->g_grid(-column=>6, -row=>9,-padx=>2,-pady=>5,-sticky=>"w");
	
	$$run_button->configure(-state=>"normal", -command =>sub{
			if (!chk_GO_db()){Tkx::tk___messageBox(-message => "Warning:E.coli GO annotation not found. Add GO enrichmentdata using Utility menu\n", -type=>"ok", -title=>"Alert", -icon=>'warning'); return 0; }
			my ($match,$mat_ids)=compare_two_fasta_file(unix_path($input_seq),unix_path($background_seq));
			#print "Percentage of match in input file iwth bck $match\n";
			
			if($match<100){ 
				Tkx::tk___messageBox(-message => "Warning:\nAll the input sequences are not present in background sequence file. $match% of target sequence found in Background sequence file", -type=>"ok", -title=>"Warning", -icon=>'warning'); 
				my $r=fetch_seq_by_id(read_fasta_sequence(unix_path($input_seq)),$mat_ids);	##
				write_fasta_seq($r,"$L_root_path\/target_match_bckg_seq.fasta");
				$input_seq = win_path("$L_root_path\/target_match_bckg_seq.fasta");
				print STDERR "$L_root_path\/target_match_bckg_seq.fasta Created\n";				
			}
			
			$$run_button->configure(-state=>"disabled");
			unlink "$L_root_path/E_coli_ortho_blast5.out.txt";
			my($gap_open, $gap_extns)=split /,/,$gap_score_5;
			my $blast5="";
			$blast5= "blastall.exe -p blastp -d $model_org_ref_proteome -i $background_seq -e $e_val_5 -m $out_fmt_5 -W $word_size_5 -M $sub_matrix_5 -G $gap_open -E $gap_extns -o $root_path\\E_coli_ortho_blast5.out -a $use_cores -f $threshold_5"." $extra_params_BLAST5" if $blast_version eq 'old_blastall';
			$blast5= "blastp -db $model_org_ref_proteome -query $background_seq -evalue $e_val_5 -outfmt $out_fmt_5 -word_size $word_size_5 -matrix $sub_matrix_5 -gapopen $gap_open -gapextend $gap_extns -out $root_path\\E_coli_ortho_blast5.out -num_threads $use_cores -threshold $threshold_5"." $extra_params_BLAST5" if $blast_version eq 'blast+';
			
			
				`echo echo off > batch.bat`;
				`echo color 90 >> batch.bat`;
				`echo cls >> batch.bat`;
				`echo echo :: External program BLAST :: >> batch.bat`;
				`echo echo ---------------------------------------------- >> batch.bat`;
				
				`echo echo Parameters >> batch.bat`;
				`echo echo 	Program: blastp >> batch.bat`;
				`echo echo 	Query		: $background_seq >> batch.bat`;
				`echo echo 	Database	: $model_org_ref_proteome >> batch.bat`;
				`echo echo 	E-value	: $e_val_5 >> batch.bat`;
				`echo echo 	Scoring matrix: $sub_matrix_5 >> batch.bat`;
				`echo echo 	Gap-penalty (Open,Extension): $gap_open,$gap_extns >> batch.bat`;
				`echo echo 	Word size : $word_size_5 >> batch.bat`;
				`echo echo 	Threshold for extending hits : $threshold_5 >> batch.bat`;
				`echo echo 	CPUs : $use_cores >> batch.bat`;
				
				
				`echo echo Please wait.......... >> batch.bat`;
				`echo $blast5 >> batch.bat`;
				`echo rename $root_path\\E_coli_ortho_blast5.out E_coli_ortho_blast5.out.txt >> batch.bat`;	##mv works
				`echo exit >> batch.bat`;
				
				if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
				else{system("start batch.bat ");}
				
				while(!(-e "$L_root_path/E_coli_ortho_blast5.out.txt")){
					$blast_prg5=blast_progress("$root_path\\E_coli_ortho_blast5.out","$L_root_path/E_coli_ortho_blast5.out",$sequence_summary{-total_seq} ); sleep(3); Tkx::update(); 
				}	##wait till blast4.out.txt is available
				if(count_blast_hits("$L_root_path/E_coli_ortho_blast5.out.txt")<1){
				Tkx::tk___messageBox(-message => "No BLAST hit found!!!Aborting analysis\nPlease change BLAST parameters and run analysis again.", -type=>"ok", -title=>"Alert",-icon=>'warning'); 
				return ();
				}
				
				$blast_prg5=100;  Tkx::update();
				my $model_org_ortho = process_GO_BLAST_out("$L_root_path/E_coli_ortho_blast5.out.txt"); ## [-seq_id => Ecoli_id_uniprt_ac]; all genes; bckgrnd
				$GO_ana_progress=10;
				my $p=0;
				my (%background_go_BP,%background_go_MF,%background_go_CC);
				foreach my $g (keys %{$model_org_ortho} ){
					my $go = fetch_GO_db($model_org_ortho->{$g});
					$GO_ana_progress=(++$p/scalar (keys %$model_org_ortho))*50;
					
					foreach my $onto_cat (keys %$go)
					{
						
						if($onto_cat eq 'biological_process') {$background_go_BP{$g}= $go->{biological_process} }
						if($onto_cat eq 'cellular_component') {$background_go_CC{$g}= $go->{cellular_component} }
						if($onto_cat eq 'molecular_function') {$background_go_MF{$g}= $go->{molecular_function} }
						Tkx::update();
					}
					
				}
				
				my @test_gene_set = (keys %{read_fasta_sequence(unix_path($input_seq))});
				my $BP_dat=calc_GO_chi(\%background_go_BP,\@test_gene_set,$L_root_path."/GO_enrich_BP","GO enrichment : Biological Process");
				$GO_ana_progress=70; Tkx::update();
				my @test_gene_set = (keys %{read_fasta_sequence(unix_path($input_seq))});
				my $MF_dat=calc_GO_chi(\%background_go_MF,\@test_gene_set,$L_root_path."/GO_enrich_MF","GO enrichment : Molecular Function");
				$GO_ana_progress=80; Tkx::update();
				my @test_gene_set = (keys %{read_fasta_sequence(unix_path($input_seq))});
				my $CC_dat=calc_GO_chi(\%background_go_CC,\@test_gene_set,$L_root_path."/GO_enrich_CC","GO enrichment : Cellular Componenet");
				$GO_ana_progress=90; Tkx::update();
				
				Tkx::image_create_photo( "BP", -file => $L_root_path."/GO_enrich_BP_enrich.gif");
				$canvas1->create_image(250, 0, -image=>"BP", -anchor =>'n' );
				Tkx::image_create_photo( "MF", -file => $L_root_path."/GO_enrich_MF_enrich.gif");
				$canvas2->create_image(250, 0, -image=>"MF", -anchor =>'n' );
				Tkx::image_create_photo( "CC", -file => $L_root_path."/GO_enrich_CC_enrich.gif");
				$canvas3->create_image(250, 0, -image=>"CC", -anchor =>'n' );				
				Tkx::update();

				$seq_used_for_analysis->configure(-state=>"normal",-command=>sub{
					my $save_result = Tkx::tk___getSaveFile();
					#$crt_win->g_raise();
					my $r=fetch_seq_by_id(read_fasta_sequence(unix_path($input_seq)),\@test_gene_set);	##
					write_fasta_seq($r,"$save_result") if $save_result ;
									
				});
				
				
				$GO_ana_progress=100; Tkx::update();
				
				Tkx::tk___messageBox(-message => "Run completete. See output folder for results.", -type=>"ok", -title=>"Success"); 
								
		});
	
	

}##END GO anal

##Args: a hash-ref of seq_id=>array_ref_of_"term|GO-id", a test seq_id arrayref, and a file prefix/path
##Returns: ref of a hash {_a=>0,_b=>0,_c=>0,_d=>0, _go=>$go_id, _bg_gene=>[], _test_gene=>[]} ; creates a output file with (_BP/CC/MP.txt), a graphics file (_BP/CC/MP.gif))
sub calc_GO_chi
{
	my $hash_ref=shift;
	my $test_id=shift;
	my $out_file=shift;
	my $GO_category=shift || "GO enrichment analysis";
	my $significance_thr=0.5;				#shift || 0.05;		## used to create Bar chart;
	
	my %categories;
	my ($a,$b,$c,$d,$df)=(0,0,0,0,1);
	
	foreach my $gene (keys %{$hash_ref}){
		foreach my $go_cat (@{$hash_ref->{$gene}})
		{
			my($term,$go_id) =split /\|/,$go_cat;
			if(!exists $categories{$term} ){$categories{$term}={_a=>0,_b=>0,_c=>0,_d=>0, _go=>$go_id, _bg_gene=>[], _test_gene=>[]};}
			$categories{$term}->{_a}++;
			push  @{$categories{$term}->{_bg_gene}},$gene;
		}	
	}
	
	
	foreach my $gene (@$test_id)
	{
		foreach my $go_cat (@{$hash_ref->{$gene}})
		{
			my($term,$go_id) =split /\|/,$go_cat;
			if(!exists $categories{$term} ){warn "$term found in test set not in back ground \n"}
			$categories{$term}->{_b}++;
			push  @{$categories{$term}->{_test_gene}},$gene;
		}	
	}
	
	foreach (keys %categories)
	{
		$categories{$_}->{_c} = scalar (keys %$hash_ref) - $categories{$_}->{_a};
		$categories{$_}->{_d} =  scalar @$test_id - $categories{$_}->{_b} ;
	}
	
	open (O, ">$out_file\_enrich.txt") or die "$! $out_file\_enrich.txt";
	print O "#Term\ttest_genes\tno_of_test_gene\tNo_of_bg_genes\tP-value\tBG_genes\n";
	
	my (@data,@categories, @freq, @p_vals);				## For GD plot
		foreach (sort keys %categories)
		{
			next if $categories{$_}->{_b} <1;			## if <1 dont print;

			#Chi-square calculation with Yates-correction
			my $chi_sq_y = chi_squared_y($categories{$_}->{_a},$categories{$_}->{_b},$categories{$_}->{_c},$categories{$_}->{_d});
			my $chi_p_val_y = chisqrprob($df,$chi_sq_y);
			#printf "Chi-square value (Y):%4.3f\nP-value: %4.3f\n\n",$chi_sq_y,$chi_p_val_y,$p_val;
			if($chi_p_val_y < $significance_thr)
			{
				push @categories,$_;
				push @freq,$categories{$_}->{_b};
				push @p_vals,sprintf("%4.3f",$chi_p_val_y) ;
			}	
			print O "$_\t".join(";",@{$categories{$_}->{_test_gene}});
			print O "\t$categories{$_}->{_b}\t$categories{$_}->{_a}\t$chi_p_val_y\t";
			print O "$_\t".join(";",@{$categories{$_}->{_bg_gene}});
			print O "\n";				
		}
		close O;
		if(scalar @categories <1){
			Tkx::tk___messageBox(-message => "No significant (P-value< $significance_thr) GO enrichment found!! ", -type=>"ok", -title=>"Alert", -icon=>'warning'); 
			return 0;
		}
			
		@data =(\@categories,\@freq);

		my $graph_chart = GD::Graph::hbars->new(1000,900);
		my $data_p = GD::Graph::Data->new([\@categories,\@p_vals] );

		$graph_chart-> set (
			x_label=> 'Categories',
			y_label=> "Frequency \n(number at the top of bar is p-value)",
			title=> $GO_category,
			bar_spacing=>1,
			show_values => $data_p,
			correct_width=>1,
			bgclr=>"white",
			accentclr=>"white",
			labelclr=>"black",
			axislabelclr=>"black",
			borderclrs=>"black",
			fgclr =>"black",
			boxclr=>"white",
			transparent=>0,	
			x_label_position =>1/2,
			
		);
		my $win_root= get_windows_root_path('L');
		$graph_chart->set_title_font("$win_root/Fonts/times.ttf",18);
		$graph_chart->set_x_label_font("$win_root/Fonts/times.ttf",14);
		$graph_chart->set_y_label_font("$win_root/Fonts/times.ttf",14);
		$graph_chart->set_x_axis_font("$win_root/Fonts/times.ttf",12);
		$graph_chart->set_y_axis_font("$win_root/Fonts/times.ttf",12);
		
		my $gd =$graph_chart-> plot(\@data) or die $graph_chart-> error;
		my $export="$out_file\_enrich.gif";
		open(IMG, ">$export") or die $!;
		binmode IMG;
		print IMG $gd->gif;
		close IMG;
		return \%categories;
}

##Meta function
##Args: a,b,c,d
##Return: Chi-square value (Yates-correction applied)
sub chi_squared_y {
     my ($a,$b,$c,$d) = @_;
     return 0 if($b+$d == 0);
     my $n= $a + $b + $c + $d;
	return (($n*(abs($a*$d - $b*$c)-($n/2))**2) / (($a + $b)*($c + $d)*($a + $c)*($b + $d)));
}

##Meta function
##Args: degree of freedom, Chi-square value
##Return: P-value as float
sub chisqrprob {
	my ($n,$x) = @_;
	my $p;

	if ($x <= 0) {
		$p = 1;
	} elsif ($n > 100) {
		$p = _subuprob((($x / $n) ** (1/3)
				- (1 - 2/9/$n)) / sqrt(2/9/$n));
	} elsif ($x > 400) {
		$p = 0;
	} else {   
		my ($a, $i, $i1);
		if (($n % 2) != 0) {
			$p = 2 * _subuprob(sqrt($x));
			$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
			$i1 = 1;
		} else {
			$p = $a = exp(-$x/2);
			$i1 = 2;
		}

		for ($i = $i1; $i <= ($n-2); $i += 2) {
			$a *= $x / $i;
			$p += $a;
		}
	}
	return $p;
}
##MEta function; function copied from CPAN Statistics:: Module 
sub _subuprob {
	my ($x) = @_;
	my $p = 0; # if ($absx > 100)
	my $absx = abs($x);

	if ($absx < 1.9) {
		$p = (1 +
			$absx * (.049867347
			  + $absx * (.0211410061
			  	+ $absx * (.0032776263
				  + $absx * (.0000380036
					+ $absx * (.0000488906
					  + $absx * .000005383)))))) ** -16/2;
	} elsif ($absx <= 100) {
		for (my $i = 18; $i >= 1; $i--) {
			$p = $i / ($absx + $p);
		}
		$p = exp(-.5 * $absx * $absx) 
			/ sqrt(2 * PI) / ($absx + $p);
	}

	$p = 1 - $p if ($x<0);
	return $p;
}


## Menu-> Downstram analysis-> Subcellular localization
##ARGS:
##returns:
sub subCellLoc_analysis
{
	my $frm=shift;
	my $run_button = shift;	
	my $frm_top=$$frm->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 500,-padding => "0 0 0 0");
	$frm_top->g_grid(-column=>0,-row=>0,-sticky=>"nswe");
	my $status=(check_internet()?"Active internet": 'No internet');
	Tkx::image_create_photo( "BANER", -file => "banner.gif");
	($frm_top->new_ttk__label(-image=>'BANER'))->g_grid(-column=>0, -row=>0,-sticky=>"nwes",-columnspan=>2, -padx=>60);	
	my $heading = $frm_top->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool",-justify=>"center",-foreground=>"blue",-font => "Helvetica 16 bold underline");
	$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>50);
	$frm_top->new_ttk__label(-text=>"SUBCELLULAR LOCALIZATION ANALYSIS",-justify=>"left",-foreground=>"darkgreen",-font => "Helvetica 12")->g_grid(-column=>0,-row=>1,-sticky=>"s",-padx=>0);
	
	my ($input_seq, $bacteria_strain);
	$input_seq = ($Tproteome_file?"$L_root_path/accepted_seq_step-4_1.fasta":"");	##Skipp if project defined;
			
	my $new_frm = $frm_top->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 200,-padding => "0 0 50 0");
	$new_frm->g_grid(-column=>0,-row=>2,-sticky=>"nswe");
	$new_frm->new_ttk__label(-text=>"Choose bacteria strain")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	$new_frm->new_ttk__combobox(-textvariable => \$bacteria_strain, -values => "None positive negative")->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	$new_frm->new_ttk__label(-text=>"Input target sequences")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	$new_frm->new_ttk__label(-text=>"Output folder")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$new_frm ->new_ttk__entry(-textvariable => \$input_seq,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	$new_frm ->new_ttk__entry(-textvariable => \$root_path,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2);
		
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$input_seq = Tkx::tk___getOpenFile(-parent=>$mw);$mw->g_raise();  		## Assuming that no project defined;
	if($input_seq){
		$sequence_summary{-putative_drug_targets}=count_fasta_seq($input_seq) ; 	## reset drug target counts; not called if in project call
		$input_seq=~ s{/}{\\}g; $input_seq='"'.$input_seq.'"'; 					##Convert to windows format
	}	
	})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
		$root_path = Tkx::tk___chooseDirectory(-parent=>$mw);$mw->g_raise();		
		if(!$root_path){$root_path=win_path($L_root_path);}  							##if cancel pressed ; fail safe
		else{
		$root_path =~s/"//g;														##incase of Default Desktop path
		$L_root_path=$root_path;													##preserve UNIX format
		$root_path=win_path($root_path); 						##Convert to windows format		
		} 
	})->g_grid(-column=>4,-row=>2,-padx=>2,-pady=>1,-sticky=>"wn");	
	
	$new_frm->new_ttk__label(-text=>"Progress:")->g_grid(-column=>0,-row=>3,-padx=>2,-pady=>5,-sticky=>"nw");
	#my $prg=$new_frm->new_ttk__progressbar(-orient => 'horizontal', -value=>0, -length=>100,-mode => 'indeterminate');
	#$prg->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>5,-sticky=>"nw");
	$new_frm->new_ttk__label(-text=>"",-foreground=>"red",-justify=>"right", -textvariable =>\$status)->g_grid(-column=>2,-row=>3,-padx=>2,-pady=>5,-sticky=>"nw", -columnspan=>2);
	
	my $canvas1 = $new_frm->new_tk__canvas(-scrollregion => "0 0 1000 1000",-width=>400, -height=>250);
	$canvas1->g_grid(-column=>0, -row=>5, -padx=>5, -pady=>1,-sticky=>"nw",-columnspan=>5);
	my $hscroll1 = $new_frm->new_tk__scrollbar(-orient => "horizontal", -command => [$canvas1, "xview"]);
	my $vscroll1 = $new_frm->new_tk__scrollbar(-orient => "vertical", -command => [$canvas1, "yview"]);
	$hscroll1->g_grid(-column => 0, -row => 6, -sticky => "we",-columnspan=>5);
	$vscroll1->g_grid(-column => 5, -row => 5, -sticky => "ns");
	$new_frm->new_ttk__sizegrip()->g_grid(-column => 5, -row => 6, -sticky => "se");
	$canvas1->configure(-yscrollcommand => [$vscroll1, "set"], -xscrollcommand => [$hscroll1, "set"]);				
	
	$$run_button->configure(-state=>"normal", -command =>sub {
			if(!$bacteria_strain || !$input_seq){ Tkx::tk___messageBox(-message => "Inputs not found!! ", -icon=>'warning'); 		return 0;}
			if($status eq 'No internet'){ Tkx::tk___messageBox(-message => "Check internet connection ", -icon=>'warning'); 		return 0;}
			$$run_button->configure(-state=>"disabled");
			
			my $size = -s unix_path($input_seq); #print "Size=$size\n";
			my $all_inp_seq = read_fasta_sequence(unix_path($input_seq)); ## each seq assuming 2000bytes * 30 =60000
			my @split_file; my @all_ids= keys %$all_inp_seq;
			if(scalar (keys %$all_inp_seq)>30){  
				my $f = int (scalar (keys %$all_inp_seq)/30); $f=(scalar (keys %$all_inp_seq)%3==0?$f:$f+1);
				for (my  $h=0; $h<$f;$h++)
				{
					my @ids= splice(@all_ids,30*$h,30);
					my $r=fetch_seq_by_id($all_inp_seq,\@ids);	##
					write_fasta_seq($r,"$L_root_path/sub_file.txt");
					psort_subcell($bacteria_strain, "$L_root_path/sub_file.txt", "$L_root_path/psort_subcellular_loc_result.txt");
					$status="$h/$f processed";
					Tkx::update();
				}
				
			}
			else{
			psort_subcell($bacteria_strain, unix_path($input_seq), "$L_root_path/psort_subcellular_loc_result.txt",\$status);
			}
			Tkx::update();
			
			my %distribution;
			my $formated_op="$L_root_path/psort_subcellular_loc_final_result.txt";
			open(O, "<$L_root_path/psort_subcellular_loc_result.txt") or die"$! $L_root_path/psort_subcellular_loc_result.txt\n";
			open(F, "> $formated_op") or die"$! $formated_op\n";
			print F "SeqID	Localization	Score\n";
			while(<O>){
				chomp;
				if(/^(SeqID)|^$/){ next;}
				else{ 
					my @l =split /\t/,$_; 
					next if scalar @l <2;
					$distribution{$l[-2]}=[] if !$distribution{$l[-2]};
					if($l[0]=~m/(\S+)/g){$l[0]=$1};
					push @{$distribution{$l[-2]}},$l[0] ;
					$"="\t";
					print F "@l\n";
				}			
			}
			close F; close O;
			my @categories=keys %distribution;
			my @freq; foreach (@categories) {push @freq, scalar @{$distribution{$_}};}
			my @data =(\@categories,\@freq);

			my $graph_chart = GD::Graph::bars->new(600,600);
			
			$graph_chart-> set (
				x_label=> 'Categories',
				y_label=> "Frequency ",
				title=> 'Subcellular localization analysis',
				bar_spacing=>1,
				correct_width=>1,
				bgclr=>"white",
				accentclr=>"white",
				labelclr=>"black",
				axislabelclr=>"black",
				borderclrs=>"black",
				fgclr =>"black",
				boxclr=>"white",
				transparent=>0,	
				x_label_position =>1/2,
				x_labels_vertical=>1,
			);
			my $win_root= get_windows_root_path('L');
			$graph_chart->set_title_font("$win_root/Fonts/times.ttf",18);
			$graph_chart->set_x_label_font("$win_root/Fonts/times.ttf",14);
			$graph_chart->set_y_label_font("$win_root/Fonts/times.ttf",14);
			$graph_chart->set_x_axis_font("$win_root/Fonts/times.ttf",12);
			$graph_chart->set_y_axis_font("$win_root/Fonts/times.ttf",12);

			my $gd =$graph_chart-> plot(\@data) or die $graph_chart-> error;
			my $export="$L_root_path/Subcellular_enrich.gif";
			open(IMG, ">$export") or die $!;
			binmode IMG;
			print IMG $gd->gif;
			close IMG;	
			Tkx::image_create_photo( "BP", -file => $L_root_path."/Subcellular_enrich.gif");
			$canvas1->create_image(400, 0, -image=>"BP", -anchor =>'n' );
			Tkx::tk___messageBox(-message => "Successful run. \nUse utlity menu to fetch dequnces by id ", -icon=>'info');
		});		
}

##MEta function
##Args; positive/negative, input_se_path_unix, outfile_path_(unix)_appends, status_ref_var
##return :
sub psort_subcell
{
	my $bacteria_strain=shift;
	my $input_seq=shift;
	my $outFile= shift;
	my $status =shift || 0;
	
			my $psort_url = 'http://www.psort.org/psortb/results.pl';
			my $browser = LWP::UserAgent->new;
			
			my $response = $browser->request(
			POST "$psort_url",
			Content_Type => 'form-data',
			Content => [
			'organism'=>'bacteria',
			'gram'=>$bacteria_strain,							#'positive',				#'negative',
			'advancedgram'=>'none',			
			'format'=>'terse',
			'sendresults'=>'display',
			'email'=>'',
			'seqs'=>'',				#$fas,
			'filename'=>[unix_path($input_seq)]
			]);
			if ($response->is_error || $response->content=~m/One or more errors occurred while processing your request/g) {$$status="ERROR in Analysis" if $status; }
			else {
			open (O, ">>$outFile")or die"$! outFile"; 
			print O $response->content;
			close O; $$status="DONE" if $status;
			}	

}

##Resets all defined variables;
##ARGS: 
##Returns:
sub reset_params()
{
	$min_aa_len=50; 
	$cd_hit_identity=60;
	$chk_cdhit=1;
	print STDERR "Reading drug Target data:";
	$drug_db_names=read_drugTarget_db("./local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt");#' {All} {drugBank} {PTTD} ';	##read files to update it;
	$ref_drug_db_array=[];
	open(G,"./local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt") or die"$! ./local_dat/KNOWN_DRUG_TARGETS/drugTarget_db_names.txt"; while(<G>){chomp; push @$ref_drug_db_array,$_;};close G;
	$drug_blast_db_names=create_drugTarget_blast_db($ref_drug_db_array,"./local_dat/KNOWN_DRUG_TARGETS");
	$drug_target_annot=read_drugTarget_annot("./local_dat/KNOWN_DRUG_TARGETS");
	print STDERR "DONE";

	print STDERR "\nReading broad-spectrum data:";
	$broad_spectrum_pathogen_db_sq = win_path($installation_path.'\local_dat\PATHOGENS\pathogen_taxonomy.db');
	##'"\"PATHOGENS\" "PATHOGENS\ACIBC" "PATHOGENS\ACIBS" "PATHOGENS\ACIBT" "PATHOGENS\BURM1\""';##read files to update it;
	$broad_spe_species_per_query=0;
	$tax_level="Family";
	$brd_sp_db_levels;
	$ref_brd_sel_db_array=[]; ## array ref ; array saves  tax_levels sel by user			
	($brd_sp_db_levels,$ref_brd_sel_db_array) = fetch_tax_names($tax_level);		
	$broad_spectrum_pathogen_db_list=create_broad_spe_db_array($ref_brd_sel_db_array,"./local_dat/PATHOGENS");
	print STDERR "DONE\n";

	$PPI_score_cutoff=700;		##STRING score cutoff
	$top_hub_perc=20;
	
	($e_val_1,$out_fmt_1,$sub_matrix_1,$gap_score_1, $extra_params_BLAST1,$word_size_1,$threshold_1,$perc_identity_1,$blast_prg1)=(0.01,8,"BLOSUM62","11,1","-b 1",3,11,20,0);	
	($e_val_2,$out_fmt_2,$sub_matrix_2,$gap_score_2, $extra_params_BLAST2,$word_size_2,$threshold_2,$perc_identity_2,$blast_prg2)=(0.0000000001,8,"BLOSUM62","11,1","-b 1",3,11,0,0);
	($e_val_3,$out_fmt_3,$sub_matrix_3,$gap_score_3, $extra_params_BLAST3,$word_size_3,$threshold_3,$perc_identity_3,$blast_prg3)=(0.01,8,"BLOSUM62","11,1","",3,11,30,0);
	($e_val_4,$out_fmt_4,$sub_matrix_4,$gap_score_4, $extra_params_BLAST4,$word_size_4,$threshold_4,$perc_identity_4,$blast_prg4)=(0.01,8,"BLOSUM62","11,1","-b 1",3,11,0,0);
	($e_val_5,$out_fmt_5,$sub_matrix_5,$gap_score_5, $extra_params_BLAST5,$word_size_5,$threshold_5,$perc_identity_5,$blast_prg5)=(0.01,8,"BLOSUM62","11,1","-b 1",3,11,0,0);

	if ($blast_version eq 'blast+'){
		($out_fmt_1,$out_fmt_2,$out_fmt_3,$out_fmt_4,$out_fmt_5)=(6,6,6,6,6);
		($extra_params_BLAST1,$extra_params_BLAST2,$extra_params_BLAST4,$extra_params_BLAST5)=("-num_alignments 1 ","-num_alignments 1 ","-num_alignments 1 ","-num_alignments 1 ");
	}	
	Tkx::tk___messageBox(-type => "ok",-message => "Parameter reset", -icon => "info", -title => "Success");	
}

### Menu ->File-> Create project function
sub create_project
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Create New Project");
	$root_path=get_my_document_path();		##OS specific; change for linux;; Just resetting
	
	open_tool_window($crt_win,$mw);
	
	
	
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	$frm1->new_ttk__label(-text=>"Create a new project directory that will contain all output files of the current analysis.",)->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>4,); #-foreground=>'red'
	$frm1->new_ttk__label(-text=>"Project directory name		")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => \$project_name,-width=>25,)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	$frm1->new_ttk__label(-text=>"Location of project directory	")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw",);
	$frm1 ->new_ttk__entry(-textvariable => \$root_path,-width=>50,-state=>"disabled",)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	$frm1->new_ttk__button(-text=>"Choose location",-command=>sub{
		$root_path = Tkx::tk___chooseDirectory(-parent=>$crt_win);									##updating root path
		#$root_path.="/$project_name";
		$crt_win->g_raise();
		Tkx::update(); 
	})->g_grid(-column=>3,-row=>2,-padx=>2,-pady=>5);
	
		my $create_proj_but=$frm1->new_button(-text=>"Create project",-width=>18,-command=>sub{
		
		$root_path.="\\$project_name";												##updating root path
		$root_path =~s/"//g;														##incase of Default Desktop path
		$root_path=win_path($root_path);
		$L_root_path=unix_path($root_path);											##preserve UNIX format; fail safe
		if(-d $L_root_path){
		my $dir_chk=Tkx::tk___messageBox(-type => "yesno",
	    -message => "A directory named $project_name already exists!!!\nDo you want to delete and create a new directory",
	    -icon => "question", -title => "Alert", -parent=>$crt_win);
		#print STDERR "YN: ".$dir_chk;
		if($dir_chk eq 'yes'){ rmdir $root_path; }
		else{return 0}		
		}
		$frm1->g_destroy;
		
		my $frm2=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
		$frm2->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
		
		$frm2->new_ttk__label(-text=>"SPECIFY INPUT FILES",-foreground=>'blue')->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>2,-sticky=>"nw");
			
		$frm2->new_ttk__label(-text=>"Pathogen proteome sequence(*.fas)[required]:		")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$Tproteome_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$Tproteome_file = Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1);
			
		$frm2->new_ttk__label(-text=>"Pre-formatted host reference proteome database[required]:	")->g_grid(-column=>0, -row=>2,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$Hproteome_file,-width=>40,-state=>"disabled")->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$Hproteome_file=Tkx::tk___getOpenFile(-parent=>$crt_win,-filetypes =>[['All Files',   '*.psd']]);
		$Hproteome_file=~s/\.psd//g;
		if($Hproteome_file and $Hproteome_file=~m/\s/g){ 
			Tkx::tk___messageBox(-type => "ok",-message => "The file path contains space that create problem with BLAST on Windows OS. Provide a non-space file path", -icon => "info", -title => "Alert", -parent=>$crt_win);
			undef $Hproteome_file;
		}		
		$crt_win->g_raise();
		})->g_grid(-column=>4,-row=>2,-padx=>2,-pady=>1);
		
		$frm2->new_ttk__label(-text=>"Pre-formatted essential protein database:")->g_grid(-column=>0, -row=>3,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$Eproteome_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$Eproteome_file=Tkx::tk___getOpenFile(-parent=>$crt_win,-filetypes =>[['All Files',   '*.psd']]);
		$Eproteome_file=~s/\.psd//g;
		if($Eproteome_file and $Eproteome_file=~m/\s/g){ 
			Tkx::tk___messageBox(-type => "ok",-message => "The file path contains space that create problem with BLAST on Windows OS.\nPlease provide a non-space file path", -icon => "info", -title => "Alert", -parent=>$crt_win);
			undef $Eproteome_file;
		}
		$crt_win->g_raise();
		})->g_grid(-column=>4,-row=>3,-padx=>2,-pady=>1);
		
		$frm2->new_ttk__label(-text=>"Comma-separated protein-protein interaction file (*.csv):		")->g_grid(-column=>0, -row=>4,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$interactome_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>4,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$interactome_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>4,-padx=>2,-pady=>1);
		
		$frm2->new_ttk__label(-text=>"Interaction ID mapping file(*.tsv): 		")->g_grid(-column=>0, -row=>5,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$PPI_id_map_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>5,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$PPI_id_map_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>5,-padx=>2,-pady=>1);
	
		$frm2->new_ttk__label(-text=>"Parameter file (optional):			")->g_grid(-column=>0, -row=>6,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$filter_param_settings_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>6,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$filter_param_settings_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>6,-padx=>2,-pady=>1);
		
		my $frm3=$frm2->new_ttk__frame(-borderwidth=>2);
		$frm3->g_grid(-row=>7,-column=>0,-sticky=>"nsew",-columnspan=>6,-pady=>15);
		
		$frm3->new_ttk__label(-text=>"ADVANCED OPTIONS",-foreground=>'blue')->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>2,-sticky=>"nw");
		$frm3->new_ttk__checkbutton(-text => "Skip CH-HIT",-variable => \$skip_cdHit, -onvalue =>1, -offvalue =>0)->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>2,-sticky=>"nw");
		$frm3->new_ttk__checkbutton(-text => "Skip Host paralog BLAST",-variable => \$skip_host_BLAST, -onvalue =>1, -offvalue =>0)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>2,-sticky=>"nw");
		
		my $species_name="";
		$frm3->new_ttk__label(-text=>"NCBI Taxon id	:")->g_grid(-column=>0, -row=>4,-padx=>2,-pady=>1,-sticky=>"e");
		$frm3 ->new_ttk__entry(-textvariable => \$taxon_id,-width=>20,)->g_grid(-column=>1,-row=>4,-padx=>2,-pady=>0,-sticky=>"w");
		$frm3->new_ttk__button(-text=>"Validate",-width=>10,-command=>sub{
			$species_name="wait..."; Tkx::update();
			my $url = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id='.$taxon_id;
			my $content = get $url;  die "Couldn't get $url . Check Internet" unless defined $content;
			my @content=split /\n/,$content;
			foreach(@content){if(/Taxonomy\s+browser\s+\((.+)\)/g){$species_name="$1";}}
			Tkx::update();		
		})->g_grid(-column=>4,-row=>4,-padx=>2,-pady=>0,-sticky=>"w");
		$frm3->new_ttk__label(-textvariable=>\$species_name)->g_grid(-column=>5, -row=>4,-padx=>2,-pady=>1,-sticky=>"e");
		
		my $data_input_but=$frm2->new_button(-text=>"Ok",-width=>10,-command=>sub{
			if(!$Tproteome_file ||(!$Hproteome_file && !$skip_host_BLAST)){
				Tkx::tk___messageBox(-type => "ok",
					-message => "Probably you forgot to load an input file(s).\n",
					-icon => "error", -title => "Input file Error", -parent=>$crt_win);
					#&create_project();
					return 0;
			}
			
			$crt_win->g_destroy;$mw->configure(-cursor=>"arrow");
			my $batch;
			if($Hproteome_file && !((-e "$Hproteome_file\.phr")and(-e "$Hproteome_file\.psq")and(-e "$Hproteome_file\.pin"))){
			$batch="\"formatdb.exe\" -i \"$Hproteome_file\" -p T  for BLASTALL version\nuse makeblastdb for BLAST+";
				Tkx::tk___messageBox(-type => "ok",
					-message => "BLAST database (*.phr,*.psq, *.pin) not found for $Hproteome_file. Use formatdb.exe/makeblastdb.exe to create it ",
					-icon => "error", -title => "Input file Error");
					#$front_page_status.="\n=> BLAST database (*.phr,*.psq, *.pin) \nnot found for host proteome ";			
			}
			if($batch){
				#$front_page_status.="\n=> Project not created.\n";	
				Tkx::tk___messageBox(-message => "Project $project_name NOT created.\n\nBLAST database (*.phr,*.psq, *.pin) not found for host proteome. $Hproteome_file", -title=>"ERROR", -parent=>$mw);				
				$run_but->configure(-state=>"disabled"); 
			}
			else{
			$run_but->configure(-state=>"normal"); 
			#$front_page_status.="\n=>Project created:$root_path\n=>Target sequence: Tproteome_file\n=>Host proteome: $Hproteome_file\n=>Essential proteome:$Eproteome_file\n=>STRING ID mapping file:$PPI_id_map_file\n=>Select parameters from 'Settings' menu and Click on Run program.\n";
						
			mkdir "$L_root_path", 0755;									##mkdir; WINDOWs specfic way //
			close_tool_window($crt_win,$mw);
			}			
			$dwn_str_anal->entryconfigure("Broadspectrum analysis",-state=>"disabled");   ##$mw is accessible so, 
			$dwn_str_anal->entryconfigure("GO analysis",-state=>"disabled");
			$dwn_str_anal->entryconfigure("Sub-cellular localization",-state=>"disabled");
			Tkx::update(); 
			Tkx::tk___messageBox(-message => "Project $project_name created.\n\nUse settings menu to change default settings and Press Run button.", -title=>"Success", -parent=>$mw);	
		});
		$data_input_but->g_grid(-column=>1, -row=>15,-padx=>5);
		my $close_proj_but=$frm2->new_button(-text=>"Close",-width=>10,-command=>sub{unlink "$root_path//$project_name"; undef $root_path; undef $project_name; Tkx::update(); close_tool_window($crt_win,$mw)});
		$close_proj_but->g_grid(-column=>2, -row=>15,-padx=>5,-sticky=>'w');
	});	
	$create_proj_but->g_grid(-column=>1, -row=>5,-pady=>10,-padx=>10);
	my $close_proj_but=$frm1->new_button(-text=>"Close",-width=>18,-command=>sub{close_tool_window($crt_win,$mw)});
	$close_proj_but->g_grid(-column=>2, -row=>5,-pady=>10);
	
}

##
sub open_project
{
}

##Menu-> Settings-> Pipeline Settings
sub settings
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Pramaeter settings");
	open_tool_window($crt_win,$mw);
	
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	my $CD_hit_identity=$frm1->new_ttk__entry(-textvariable => \$cd_hit_identity, -width=>5,-state=>"normal",);
	$frm1->new_ttk__label(-text=>"Minimum amino acid (translated nucleotide) sequence length:",)->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__entry(-textvariable => \$min_aa_len, -width=>8)->g_grid(-column=>1,-row=>0,-padx=>1,-pady=>5,-sticky=>"nw");
	#$frm1->new_ttk__entry(-textvariable => \$min_ntd_len, -width=>8)->g_grid(-column=>2,-row=>0,-padx=>1,-pady=>5,-sticky=>"w");
	$frm1->new_ttk__checkbutton(-text => "Enable sequence clustering", -command => sub {$CD_hit_identity->configure(-state=>'normal') if $chk_cdhit; $CD_hit_identity->configure(-state=>'disable') if !$chk_cdhit; Tkx::update();},
	    -variable => \$chk_cdhit, -onvalue => 1, -offvalue =>0)->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"       CDHit identity threshold (60%-100%):",)->g_grid(-column=>0,-row=>3,-padx=>2,-pady=>5,-sticky=>"nw");	
	$CD_hit_identity->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>5,-sticky=>"w");
	#$frm1->new_ttk__label(-text=>"%",-foreground=>'red')->g_grid(-column=>2,-row=>3,-padx=>0,-pady=>5,-sticky=>"w");
	$frm1->new_ttk__label(-text=>"BLAST settings against host reference sequence:",)->g_grid(-column=>0,-row=>4,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_button(-text=>"Change/view",-width=>10,-command=>sub{
		view_update_blast_params ("BLASTp with host proteome parameters", \$crt_win, $Tproteome_file, $Hproteome_file, \$e_val_1, \$out_fmt_1, \$word_size_1, \$sub_matrix_1, \$gap_score_1, \$threshold_1,\$perc_identity_1,\$extra_params_BLAST1);
			
	},)->g_grid(-column=>1,-row=>4,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	
	$frm1->new_ttk__label(-text=>"BLAST settings against essential proteome:",)->g_grid(-column=>0,-row=>5,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_button(-text=>"Change/view",-width=>10,-command=>sub {
	
		view_update_blast_params ("BLASTp with host proteome parameters", \$crt_win, "$root_path\\accepted_seq_step-3.fasta", $Eproteome_file, \$e_val_2, \$out_fmt_2, \$word_size_2, \$sub_matrix_2, \$gap_score_2, \$threshold_2,\$perc_identity_2,\$extra_params_BLAST2);
	},)->g_grid(-column=>1,-row=>5,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
		
	$frm1->new_ttk__label(-text=>"Number of processors use for BLAST searches:")->g_grid(-column=>0,-row=>15,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__combobox(-textvariable => \$use_cores, -values=>$CPU_list,-width=>5)->g_grid(-column=>1,-row=>15,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$frm1->new_ttk__label(-text=>"Protein-protein interaction score cutoff:",)->g_grid(-column=>0,-row=>16,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_ttk__entry(-textvariable => \$PPI_score_cutoff, -width=>3,-state=>"normal",)->g_grid(-column=>1,-row=>16,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_ttk__label(-text=>"Top #% of hub proteins:",)->g_grid(-column=>0,-row=>17,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_ttk__entry(-textvariable => \$top_hub_perc, -width=>3,-state=>"normal",)->g_grid(-column=>1,-row=>17,-padx=>2,-pady=>5,-sticky=>"nw");



	
	my $apply_change=$frm1->new_button(-text=>"OK",-width=>8,-command=>sub {$crt_win->g_destroy;close_tool_window($crt_win,$mw);})->g_grid(-column=>0,-row=>25,-padx=>18,-pady=>12,-columnspan=>2,);
	
}


##args:title, parent window ref, query name,db name, ref.e-val, ref.outformat, ref.word size, ref.sub matrix, ref.gap score,ref.Threshold, ref.blast_identy,ref.extra
##returns: none; updates all values parsed to it
#view_update_blast_params ($title, $parent, $q, $db, $e_val, $outfmt, $word_size, $sub_mat, $gap_score, $threshold,$extra)
sub view_update_blast_params 
{
	my ($title, $parent, $q, $db, $e_val, $outfmt, $word_size, $sub_mat, $gap_score, $threshold,$perc_identity,$extra)=@_;	
	my $crt_win_blast2 =$$parent->new_toplevel();
	$crt_win_blast2->g_wm_title("$title");
	open_tool_window($crt_win_blast2,$$parent);
	
	my $frm1=$crt_win_blast2->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	#host, query, e-value, matrix, word size, gap open, gap extend, output format
	$frm1->new_ttk__label(-text=>"BLAST version:")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"$blast_version")->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-sticky=>"nw");
	
	$frm1->new_ttk__label(-text=>"Database:")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => \$db,-width=>35,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1->new_ttk__label(-text=>"Query:")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => \$q,-width=>35,-state=>"disabled",)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1->new_ttk__label(-text=>"E-val threshold:")->g_grid(-column=>0,-row=>3,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => $e_val,-width=>5,)->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Output format:")->g_grid(-column=>0,-row=>4,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => $outfmt,-width=>3,-state=>"disabled",)->g_grid(-column=>1,-row=>4,-padx=>2,-pady=>1,-sticky=>"nw");
	
	$frm1->new_ttk__label(-text=>"Threshold to extend a hit:")->g_grid(-column=>0,-row=>5,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => $threshold,-width=>3,)->g_grid(-column=>1,-row=>5,-padx=>2,-pady=>1,-sticky=>"nw");
		
	#, gap open, gap extend, output format
	$frm1->new_ttk__label(-text=>"Word size:")->g_grid(-column=>0,-row=>6,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => $word_size,-width=>5,)->g_grid(-column=>1,-row=>6,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Substitution matrix:")->g_grid(-column=>0,-row=>7,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1->new_ttk__combobox(-textvariable => $sub_mat,-values=> "BLOSUM62 PAM150 BLOSUM45", -width=>12,-state=>'readonly')->g_grid(-column=>1,-row=>7,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Gap penalty (opening,extension):")->g_grid(-column=>0,-row=>8,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1->new_ttk__combobox(-textvariable => $gap_score,-values=> "11,2 10,2 9,2 8,2 7,2 6,2 13,1 12,1 11,1 10,1 9,1", -width=>12,-state=>'readonly')->g_grid(-column=>1,-row=>8,-padx=>2,-pady=>1,-sticky=>"nw");
	#$frm1 ->new_ttk__entry(-textvariable => \$,-width=>30,-state=>"disabled",)->g_grid(-column=>1,-row=>8,-padx=>2,-pady=>1);
	$frm1->new_ttk__label(-text=>"Percent identity:")->g_grid(-column=>0,-row=>9,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1->new_ttk__entry(-textvariable => $perc_identity, -width=>3,-state=>"normal",)->g_grid(-column=>1,-row=>9,-padx=>2,-pady=>5,-sticky=>"nw");	;
	
	my $BLAST2_params_ele=$frm1->new_ttk__entry(-textvariable => $extra, -width=>50,-state=>"disabled",);
	my $extra_params_BLAST2_chk=0; 
	$frm1->new_ttk__checkbutton(-text => "Apply additional BLAST parameters", -command => sub {
			$BLAST2_params_ele->configure(-state=>"normal") if $extra_params_BLAST2_chk; 
			$BLAST2_params_ele->configure(-state=>"disabled") if !$extra_params_BLAST2_chk; 
			Tkx::update();
			},
	    -variable => \$extra_params_BLAST2_chk, -onvalue => 1, -offvalue =>0)->g_grid(-column=>0,-row=>10,-padx=>2,-pady=>5,-sticky=>"nw");
		
		$frm1->new_ttk__label(-text=>"Enter extra BLAST parameters:")->g_grid(-column=>0,-row=>11,-padx=>2,-pady=>1,-sticky=>"nw");
		$BLAST2_params_ele->g_grid(-column=>0,-row=>12,-padx=>2,-pady=>1,-columnspan=>2);
		
		my $apply_change=$frm1->new_button(-text=>"Close",-width=>5,-command=>sub 
		{	close_tool_window($crt_win_blast2,$$parent);
		})->g_grid(-column=>1,-row=>15,-padx=>2,-pady=>5,-sticky=>"nw");

}


##Menu-> settings->Downstream analysis
##Args:
##return:
sub down_str_anal
{

	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Down-stream analysis Settings");
	open_tool_window($crt_win,$mw);
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	$frm1->new_ttk__label(-text=>"BLAST settings for broad-spectrum analysis:",)->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_button(-text=>"Change/view",-width=>10,-command=>sub{
		view_update_blast_params ("BLASTp broad-spectrum analysis", \$crt_win, "Input fasta file", "broad-spectrum database", \$e_val_3, \$out_fmt_3, \$word_size_3, \$sub_matrix_3, \$gap_score_3, \$threshold_3,\$perc_identity_3,\$extra_params_BLAST3);	
	},)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	
	
	$frm1->new_ttk__label(-text=>"BLAST identity threshold for broad-spectrum analysis(%):",)->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");	
	my $broad_spe_BLAST_identity=$frm1->new_ttk__entry(-textvariable => \$perc_identity_3, -width=>3,-state=>"normal",);
	$broad_spe_BLAST_identity->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");	
	##BROAD_spe DATABASE
	$frm1->new_ttk__label(-text=>"Select Pathogens range",)->g_grid(-column=>0,-row=>3,-padx=>2,-pady=>0,-sticky=>"sw");	
	$frm1->new_ttk__label(-text=>"(Select a taxnomy level click Apply.\nUse Ctrl to select multiples)",)->g_grid(-column=>0,-row=>4,-padx=>2,-pady=>0,-sticky=>"wn");	
	$frm1->new_ttk__combobox(-textvariable =>\$tax_level ,-values=> "Phylum Class Order_ Family Genus", -width=>12,-state=>'readonly')->g_grid(-column=>0,-row=>5,-padx=>2,-pady=>1,-sticky=>"nw");
	#my  $brd_sp_db_levels;		#	' {All} {drugBank} {PTTD} ';
	my $sel_entries;			## stores a ref to array; ech entru cooresponds to a entryin list; so @idx match to list array
	my $u8;						## just a temp variable
	($u8,$sel_entries) = fetch_tax_names($tax_level);		
	$frm1->new_button(-text=>"Apply",-width=>10,-command=>sub 
		{	
			($brd_sp_db_levels,$sel_entries) = fetch_tax_names($tax_level);	
			
			
		})->g_grid(-column=>0,-row=>6,-padx=>1,-pady=>5,-sticky=>"nw");	
	
	#we may show a textbox with number of species selected  scalar @ref_brd_sel_db_array
	#$frm1->new_ttk__entry(-textvariable => $extra, -width=>50,-state=>"disabled",);
	#$broad_spectrum_pathogen_db_list
	my $lbox=$frm1->new_tk__listbox(-height =>4,-listvariable => \$brd_sp_db_levels,-selectmode=>'extended');	
	$lbox->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>5,-rowspan => 4,-columnspan=>3, -sticky => "nsew");	##-selectmode=>'browse'	
	$lbox->selection_set(0,scalar@$ref_brd_sel_db_array);
	my @sel;
	$lbox->g_bind("<<ListboxSelect>>", sub {
			my @idx = split /\s+/,$lbox->curselection;
			if(scalar@idx ==1){@sel=()}						## fix if a single is clicked
			foreach my $u(@idx){push @sel,@{$sel_entries}[$u];}
			$ref_brd_sel_db_array=\@sel;
			$broad_spectrum_pathogen_db_list=create_broad_spe_db_array($ref_brd_sel_db_array,"./local_dat/PATHOGENS");	
		});	
	
	
	$frm1->new_ttk__label(-text=>"Select known drug target database(s)",)->g_grid(-column=>0,-row=>14,-padx=>2,-pady=>0,-sticky=>"sw");	
	$frm1->new_ttk__label(-text=>"(use Ctrl or Select button to select multiple databases)",)->g_grid(-column=>0,-row=>16,-padx=>2,-pady=>0,-sticky=>"n");	
	my $lbox2=$frm1->new_tk__listbox(-height =>4,-listvariable => \$drug_db_names,-selectmode=>'extended');
	$lbox2->g_grid(-column=>1,-row=>14,-padx=>2,-pady=>5,-rowspan => 4,-columnspan=>3, -sticky => "nsew");	##-selectmode=>'browse'
	my $sel_drug_tar_db;
	$lbox2->selection_set(0,scalar@$ref_drug_db_array);
	$lbox2->g_bind("<<ListboxSelect>>", sub {
			my @idx = $lbox2->curselection;
			$sel_drug_tar_db=join(",",@idx);
			$ref_drug_db_array=\@idx;
			$drug_blast_db_names=create_drugTarget_blast_db($ref_drug_db_array,"./local_dat/KNOWN_DRUG_TARGETS");
	
		});
	$frm1->new_ttk__label(-text=>"BLAST settings against drug-target database:",)->g_grid(-column=>0,-row=>20,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_button(-text=>"Change/view",-width=>10,-command=>sub{
		view_update_blast_params ("BLASTp against drug-target database", \$crt_win, "Input.fasta", "Drug target list", \$e_val_4, \$out_fmt_4, \$word_size_4, \$sub_matrix_4, \$gap_score_4, \$threshold_4,\$perc_identity_4,\$extra_params_BLAST4);	
	
	},)->g_grid(-column=>1,-row=>20,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	
	$frm1->new_ttk__label(-text=>"BLAST settings for enrichment analysis:",)->g_grid(-column=>0,-row=>21,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_button(-text=>"Change/view",-width=>10,-command=>sub{
		view_update_blast_params ("BLASTp against E.coli proteome database", \$crt_win, "Input.fasta", "E.coli proteome", \$e_val_5, \$out_fmt_5, \$word_size_5, \$sub_matrix_5, \$gap_score_5, \$threshold_5,\$perc_identity_5,\$extra_params_BLAST5);	
	
	},)->g_grid(-column=>1,-row=>21,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	
	my $ok=$frm1->new_button(-text=>"Close",-width=>10,-command=>sub 
		{
		$crt_win->g_destroy();close_tool_window($crt_win,$mw);
		})->g_grid(-column=>0,-row=>30,-padx=>1,-pady=>5,-sticky=>"ne");	
}


##Menu-> settings-> System preferences; update later with window color
##Args:
##Returns:
sub sys_preferences
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("System preference Settings");
	open_tool_window($crt_win,$mw);
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	
	$frm1->new_ttk__checkbutton(-text => " Hide command prompt on external tool run",-variable => \$cmd_hide, -onvalue =>1, -offvalue =>0)->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>2,-columnspan=>2,-sticky=>"news");
	$frm1->new_ttk__checkbutton(-text => "Show welcome message at program start",-variable => \$wlc_msg_show, -onvalue =>1, -offvalue =>0)->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>2,-columnspan=>2,-sticky=>"news");
	$frm1->new_ttk__label(-text=>"BLAST verion (Not working): ",)->g_grid(-column=>0,-row=>3,-padx=>1,-pady=>5,-sticky=>"nw");	
	my $b_opt=$frm1->new_ttk__combobox(-textvariable =>\$blast_version ,-values=> "old_blastall blast+", -width=>12,-state=>'readonly');
	$b_opt->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>1,-sticky=>"nw");
	$b_opt->g_bind("<<ComboboxSelected>>", sub {
		if ($blast_version eq 'blast+'){
		($out_fmt_1,$out_fmt_2,$out_fmt_3,$out_fmt_4,$out_fmt_5)=(6,6,6,6,6);
		($extra_params_BLAST1,$extra_params_BLAST2,$extra_params_BLAST4,$extra_params_BLAST5)=("-num_alignments 1 ","-num_alignments 1 ","-num_alignments 1 ","-num_alignments 1 ");
		}
		elsif($blast_version eq 'old_blastall'){
		($out_fmt_1,$out_fmt_2,$out_fmt_3,$out_fmt_4,$out_fmt_5)=(8,8,8,8,8);
		($extra_params_BLAST1,$extra_params_BLAST2,$extra_params_BLAST4,$extra_params_BLAST5)=("-b 1 ","-b 1 ","-b 1 ","-b 1 ");
		}
	});
	
	
	my $ok=$frm1->new_button(-text=>"Close",-width=>10,-command=>sub 
		{
		$crt_win->g_destroy();close_tool_window($crt_win,$mw);
		open (O, ">./local_dat/sys_conf") or die "$! ./local_dat/sys_conf\n";
		print O "SYS_CONF_HIDE_CMD_PROMPT= $cmd_hide\n";
		print O "SYS_CONF_SHOW_WELCOME_MSG= $wlc_msg_show\n";
		print O "SYS_CONF_BLAST_VER= $blast_version\n";
		close O;		
		})->g_grid(-column=>0,-row=>30,-padx=>1,-pady=>5,-sticky=>"ne");	
}



################################################################################################
#  Menu-> Utility
################################################################################################

#ARGS:
#Return:
sub util_PPI
{

	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Utility:Create PPI inputs");
	open_tool_window($crt_win,$mw);
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	my ($pathogen_seq, $PPI_seq, $PPI_interaction,$output_dir,$total_seq_count);
	$frm1->new_ttk__label(-text=>"Pathogen proteome")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"PPI proteome sequnnce")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_ttk__label(-text=>"PPI interction file")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Output directory")->g_grid(-column=>0,-row=>3,-padx=>2,-pady=>5,-sticky=>"nw");
	
	
	
	$frm1 ->new_ttk__entry(-textvariable => \$pathogen_seq,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$PPI_seq,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$PPI_interaction,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$output_dir,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>1,-columnspan=>2);
		
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$pathogen_seq = Tkx::tk___getOpenFile(-parent=>$crt_win);
	$total_seq_count =count_fasta_seq($pathogen_seq) if $pathogen_seq; 							## reset drug target counts; not called if in project call
	$pathogen_seq=win_path($pathogen_seq) if $pathogen_seq; 							##Convert to windows format
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$PPI_seq = Tkx::tk___getOpenFile(-parent=>$crt_win);
	
	$PPI_seq=win_path($PPI_seq) if $PPI_seq; 							##Convert to windows format
	if($PPI_seq=~m/\s/g){Tkx::tk___messageBox(-message => "Input file $PPI_seq contains space. Cannot use it in Formatdb/makeblastdb program on Windows OS. Please provide a non-space file path", -icon=>'warning', -type=>"ok", -title=>"Alert");
	$PPI_seq="";
	}
	
	})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$PPI_interaction = Tkx::tk___getOpenFile(-parent=>$crt_win);
	
	})->g_grid(-column=>4,-row=>2,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
		$output_dir = Tkx::tk___chooseDirectory(-parent=>$crt_win);
		$output_dir = win_path($output_dir) if $output_dir;
		
	})->g_grid(-column=>4,-row=>3,-padx=>2,-pady=>1,-sticky=>"wn");	
	
	my $PPI_blast_prg=0;
	my $run_but=$frm1->new_button(-text=>"Run",-width=>8, -command=>sub{
	
		if ((!$pathogen_seq|| !$PPI_seq) || (!$PPI_interaction||!$output_dir)){ Tkx::tk___messageBox(-message => "ERROR: Input file missing", -icon=>'warning', -type=>"ok", -title=>"Alert");  $crt_win->g_destroy();	\&util_PPI(); return 0;}
		## perform a blast to find one-to-one realation
		##create the id mapping file
		## create non-reduntant comma separated PPI file.
		$PPI_blast_prg=5;
		if(!-e unix_path($PPI_seq).".pin") { 
		system("formatdb.exe -p T -i ".$PPI_seq) if ($blast_version eq 'old_blastall');
		system("makeblastdb.exe -dbtype prot -in ".$PPI_seq) if ($blast_version eq 'blast+');
		}
		unlink unix_path($output_dir."\\PPI_blast.out.txt");
		if ($blast_version eq 'old_blastall'){
		my $blastcmd="blastall.exe -p blastp -d $PPI_seq -i $pathogen_seq -m 8 -W 7 -b 1 -a $use_cores -o $output_dir\\PPI_blast.out";		
				`echo echo off > batch.bat`;
				`echo color 90 >> batch.bat`;
				`echo cls >> batch.bat`;
				`echo echo :: External program BLAST :: >> batch.bat`;
				`echo echo ---------------------------------------------- >> batch.bat`;
				
				`echo echo Parameters >> batch.bat`;
				`echo echo 	Program: blastp >> batch.bat`;
				`echo echo 	Query		: $pathogen_seq >> batch.bat`;
				`echo echo 	Database	: $PPI_seq >> batch.bat`;
				`echo echo 	E-value	: $e_val_5 >> batch.bat`;
				`echo echo 	Scoring matrix: $sub_matrix_5 >> batch.bat`;
				
				`echo echo 	Word size : 7 >> batch.bat`;
				`echo echo 	Threshold for extending hits : $threshold_5 >> batch.bat`;
				`echo echo 	CPUs : $use_cores >> batch.bat`;
				
				
				`echo echo Please wait.......... >> batch.bat`;
				`echo $blastcmd >> batch.bat`;
				`echo rename  $output_dir\\PPI_blast.out PPI_blast.out.txt >> batch.bat`;				
				`echo exit >> batch.bat`;				
				if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
				else{system("start batch.bat ");}
		}
		elsif($blast_version eq 'blast+'){
		my $blastcmd="blastp -db $PPI_seq -query $pathogen_seq -outfmt 6 -word_size 7 -num_alignments 1 -num_threads $use_cores -out $output_dir\\PPI_blast.out";
			
				`echo echo off > batch.bat`;
				`echo color 90 >> batch.bat`;
				`echo cls >> batch.bat`;
				`echo echo :: External program BLAST :: >> batch.bat`;
				`echo echo ---------------------------------------------- >> batch.bat`;
				
				`echo echo Parameters >> batch.bat`;
				`echo echo 	Program: blastp >> batch.bat`;
				`echo echo 	Query		: $pathogen_seq >> batch.bat`;
				`echo echo 	Database	: $PPI_seq >> batch.bat`;
				`echo echo 	E-value	: $e_val_5 >> batch.bat`;
				`echo echo 	Scoring matrix: $sub_matrix_5 >> batch.bat`;
				
				`echo echo 	Word size : 7 >> batch.bat`;
				`echo echo 	Threshold for extending hits : $threshold_5 >> batch.bat`;
				`echo echo 	CPUs : $use_cores >> batch.bat`;
				
				
				`echo echo Please wait.......... >> batch.bat`;
				`echo $blastcmd >> batch.bat`;
				`echo rename  $output_dir\\PPI_blast.out PPI_blast.out.txt >> batch.bat`;				
				`echo exit >> batch.bat`;				
				if($cmd_hide){ system("wscript.exe HideCmd.vbs batch.bat ");}
				else{system("start batch.bat ");}		
		}		
		$PPI_blast_prg=50;	Tkx::update(); 
		
		while(!(-e unix_path($output_dir."\\PPI_blast.out.txt"))){sleep(2);}
		
		my $mapped_hash = process_GO_BLAST_out(unix_path($output_dir."\\PPI_blast.out.txt"),98);
		$PPI_blast_prg=55;	Tkx::update(); 
		
		open (O,">",unix_path($output_dir)."/PPI_ID_mapped.txt") or die "$! ".unix_path($output_dir)."/PPI_ID_mapped.txt";
		foreach (keys %$mapped_hash) {print O "$_\t$mapped_hash->{$_}\n"}
		close O;
		$PPI_blast_prg=80;	Tkx::update(); 
		
		open (P,"<$PPI_interaction") or die "$! $PPI_interaction";
		my %PPI_hash;
		while(<P>){ chomp; my ($A,$B,$score)=split /\s+/,$_; if(!$PPI_hash{$A."*".$B} and !$PPI_hash{$B."*".$A}){$PPI_hash{$A."*".$B} = $score}; }
		close P;
		$PPI_blast_prg=95;	Tkx::update(); 
		open (O,">",unix_path($output_dir)."/PPI_interaction_with_score.csv") or die "$! ".unix_path($output_dir)."/PPI_interaction_with_score.csv";
		foreach (keys %PPI_hash) { my ($A,$B)=split /\*/,$_; print O "$A,$B,$PPI_hash{$_}\n";}
		close O;
		$PPI_blast_prg=100;	Tkx::update(); 
		
		close_tool_window($crt_win,$mw);
			
	});
	$run_but->g_grid(-column=>0, -row=>6,-padx=>5,-sticky=>"w");
	my $PPI_blast=$frm1->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$PPI_blast_prg);
	$PPI_blast -> g_grid(-column=>2,-row=>6,-padx=>0,-pady=>1,-columnspan=>2,-sticky=>"w");
}

##
##Args
##returns
sub tool_fetch_id_seq
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Utility:subset fasta file");
	open_tool_window($crt_win,$mw);
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	my ($ids, $fasta_file,$total_seq_count,$tot_id_count);
	
	$frm1->new_ttk__label(-text=>"Paste IDs (one id per line)  using Ctrl+V")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>3,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"(IDs are the first stretch of \naplhanumeric characters without space after >)")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>0,-sticky=>"nw");
	
	my $frm2=$frm1->new_ttk__frame();
	$frm2->g_grid(-row=>0,-column=>1,-sticky=>"nsew",-rowspan=>4);	
	(my $tb = $frm2->new_tk__text(-width => 20, -height => 10))->g_grid(-column => 1, -row => 0, -sticky => "nwes",-padx=>2,-pady=>5,);	
	(my $s = $frm2->new_ttk__scrollbar(-command => [$tb, "yview"], 
        -orient => "vertical"))->g_grid(-column =>2, -row => 0, -sticky => "ns",-rowspan=>4,-padx=>0,-pady=>5 );
	$tb->configure(-yscrollcommand => [$s, "set"]);
	
		
	
	
	$frm1->new_ttk__label(-text=>"FASTA file ")->g_grid(-column=>0,-row=>6,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => \$fasta_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>6,-padx=>2,-pady=>1,-columnspan=>2);
	
	$frm1->new_ttk__button(-text=>"...",-width=>3,-command=>sub{
	$fasta_file = Tkx::tk___getOpenFile(-parent=>$crt_win);
	$total_seq_count =count_fasta_seq($fasta_file) if $fasta_file; 							
	})->g_grid(-column=>3,-row=>6,-padx=>2,-pady=>1,-sticky=>"ns");
	
	
	(my $mat=$frm1->new_ttk__button(-state=>"disabled",-text=>"ID matched sequences",-width=>22,-command=>sub{}))->g_grid(-column=>0,-row=>7,-padx=>5,-pady=>1,-sticky=>"wn");
	(my $Nmat=$frm1->new_ttk__button(-state=>"disabled",-text=>"IDs missing in Fasta",-width=>22,-command=>sub{}))->g_grid(-column=>0,-row=>8,-padx=>5,-pady=>4,-sticky=>"wn");
	
	(my $Mmat=$frm1->new_ttk__button(-state=>"disabled",-text=>"NOT matched sequences",-width=>22,-command=>sub{}))->g_grid(-column=>0,-row=>9,-padx=>5,-pady=>4,-sticky=>"wn");
	
	
	$frm1->new_ttk__button(-text=>"Fetch",-width=>8,-command=>sub{
		$ids = $tb->get("1.0", "end");
		my @ids=split/\n/,$ids;
		if ((scalar@ids<0|| !$fasta_file) && !$total_seq_count ){ Tkx::tk___messageBox(-message => "ERROR: Input file missing", -icon=>'warning', -type=>"ok", -title=>"Alert",-parent=>$crt_win);  return 0;}		
		
		
		my(@matched,@not_matched,@IDS_not_in_tot,%ids);
		my $tot=read_fasta_sequence($fasta_file);
		foreach my $i (@ids)
		{
			next $i =~ /^$/;
			next if $ids{$i};
			if($tot->{$i}){ push @matched,$i; }
			else{push @IDS_not_in_tot,$i;}
			$ids{$i}=1;
		}
		foreach my $i (keys %$tot)
		{
			if(!exists $ids{$i}){ push @not_matched,$i; }
		}
	
		$mat->configure(-state=>"normal", -command=>sub{
			my $file=Tkx::tk___getSaveFile(-parent=>$crt_win);
			if($file)
			{
			my $r= fetch_seq_by_id($tot,\@matched);
			write_fasta_seq($r,$file);
			}			
		});
		$Nmat->configure(-state=>"normal", -command=>sub{
			my $file=Tkx::tk___getSaveFile(-parent=>$crt_win);
			if($file)
			{
			my $r= fetch_seq_by_id($tot,\@not_matched);
			write_fasta_seq($r,$file);
			}			
		});	
		$Mmat->configure(-state=>"normal", -command=>sub{
			my $file=Tkx::tk___getSaveFile(-parent=>$crt_win);
			if($file)
			{
			open (O, ">$file") or die "$! $file\n";
			$"="\n";
			print O "#IDS_not_in_tot\n\n";
			print O "@IDS_not_in_tot\n";
			close O;
			$"=" ";
			}			
		});
			
		my $id_c=scalar (keys %ids); my $m_c=@matched; my $Nm_c=@not_matched; my $mm_c=@IDS_not_in_tot;
		Tkx::tk___messageBox(-message => "Success\n\nSummary:\nID count: $id_c , Total sequences: $total_seq_count\nMatched ids: $m_c\nNot_matched: $Nm_c\n Ids missing in Fasta: $mm_c\n\n", -icon=>'info', -type=>"ok", -title=>"Success",-parent=>$crt_win);
	
	})->g_grid(-column=>2,-row=>16,-padx=>2,-pady=>5,-sticky=>"wn");
	$frm1->new_ttk__button(-text=>"Close",-width=>8,-command=>sub{
	close_tool_window($crt_win,$mw);
	})->g_grid(-column=>3,-row=>16,-padx=>2,-pady=>5,-sticky=>"wn");
}

###
sub add_to_broad_spectrum_db
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Utility:add pathogen(s) to broadspectrum db");
	open_tool_window($crt_win,$mw);
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	
	my ($taxonomy_file, $fasta_dir);
	$frm1->new_ttk__label(-text=>"Taxonomy file")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Directory containing FASTA files")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$frm1 ->new_ttk__entry(-textvariable => \$taxonomy_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$fasta_dir,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$taxonomy_file = Tkx::tk___getOpenFile(-parent=>$crt_win);
	
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
		$fasta_dir = Tkx::tk___chooseDirectory(-parent=>$crt_win);
		$fasta_dir = win_path($fasta_dir) if $fasta_dir;		
	})->g_grid(-column=>4,-row=>1,-padx=>1,-pady=>1,-sticky=>"wn");	
	
	my $brd_proc_prg=0; my $c=0; my $total_lines;
	my $run_but=$frm1->new_button(-text=>"Run",-width=>8, -command=>sub{
		if (!$taxonomy_file || !$fasta_dir){ Tkx::tk___messageBox(-message => "ERROR: Input file missing", -type=>"ok", -title=>"Alert", -icon=>'warning');  $crt_win->g_destroy();	\&add_to_broad_spectrum_db(); return 0;}
	open (FILE, $taxonomy_file) or die "$! $taxonomy_file";
	 $total_lines =@{[<FILE>]};
	close FILE;
		my $BLAST_DB_DIR='PATHOGENS';
		my $database = './local_dat/PATHOGENS/pathogen_taxonomy.db';
		
		open (T, "< $taxonomy_file") or die "$! $taxonomy_file";
		while(<T>){
			chomp;
			 if (/^#/ || /^$/){next;$total_lines--; }
			~s/'//g;
			my @l=split /\t/,$_;
			if(!(-e unix_path($fasta_dir."\\$l[1]"))) { warn unix_path($fasta_dir."\\$l[1]"),"$fasta_dir\\$l[1] FASTA file for $l[0] Not found. Skip\n";--$total_lines; next;}
			open (L,unix_path($fasta_dir."\\$l[1]")) or die "$! $fasta_dir\\$l[1]";
			while(<L>){if(/^>(\S+)/){my ($a,$b,$c)= split /\|/,$1; my ($a,$b) = split /\_/,$c;  if($l[0]=~m/(\S+)\.\w+$/g){ if($1 ne $b){warn "Fasta header does not match file name: $1 ne $b\nTypically Uniprot proteome files expected\nExample: \nFASTA header:tr|Q6FCC9|Q6FCC9_ACIAD\nFile name: ACIAD.fasta. ";}  }  last;} }
			close L;
			##tr|Q6FCC9|Q6FCC9_ACIAD format is fixed;  put exit if not met
			
			copy unix_path($fasta_dir."\\$l[1]"), "./local_dat/PATHOGENS/$l[1]";			
			system ("formatdb.exe -p T -i .\\local_dat\\$BLAST_DB_DIR\\$l[1]") if $blast_version eq 'old_blastall';
			system("makeblastdb.exe -dbtype prot -in .\\local_dat\\$BLAST_DB_DIR\\$l[1]") if ($blast_version eq 'blast+');
			if (scalar @l <5) {warn "insufficient columns in line @l\nSkip.."; next;}
			$l[5]=~s/'/`/g;
		 my ($SuperKingdom,$phylum, $class, $order, $family, $genus) = split /;\s*/,$l[3];
		 my $sql = "INSERT INTO taxonomy (fasta_file, ORG_CODE,species,TaxID,SuperKingdom,phylum, class, order_, family, genus, Description) VALUES ('$l[1]','$l[0]','$l[2]',$l[4],'$SuperKingdom','$phylum', '$class', '$order', '$family', '$genus','$l[5]' )";	
		# print "$sql\n";
					my $driver   = "SQLite";
					my $dsn = "DBI:$driver:dbname=$database";
					my $userid = "";
					my $password = "";
					my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 }) or die $DBI::errstr;
					my $rv = $dbh->do($sql) or sub {warn "$l[0] : ",$DBI::errstr; next;};	
					
		$brd_proc_prg=(++$c/$total_lines)*100;	Tkx::update(); 	
		}
		close T;
		close_tool_window($crt_win,$mw);
		reset_params();
		Tkx::tk___messageBox(-message => "$total_lines imported successfully", -type=>"ok", -title=>"Success", -icon=>'info'); 
	
	});
	$run_but->g_grid(-column=>0, -row=>6,-padx=>5,-sticky=>"w");
	my $brd_sp_add=$frm1->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$brd_proc_prg);
	$brd_sp_add -> g_grid(-column=>2,-row=>6,-padx=>0,-pady=>1,-columnspan=>2,-sticky=>"w");
}


sub add_a_drug_target_db
{

	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Utility:add known drug-target database");
	open_tool_window($crt_win,$mw);
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	
	my ($drug_targ_seq_file, $annot_file,$prefix);
	$frm1->new_ttk__label(-text=>"Drug-target FASTA sequence file")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Drug-target annotation file")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Drug-target Name")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$frm1 ->new_ttk__entry(-textvariable => \$drug_targ_seq_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$annot_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$prefix,-width=>30,)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2,-sticky=>"wn");
	
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$drug_targ_seq_file = Tkx::tk___getOpenFile(-parent=>$crt_win);
	
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$annot_file = Tkx::tk___getOpenFile(-parent=>$crt_win);
	
	})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1,-sticky=>"wn");
	
	my $drg_tar_prg=0; my $c=0;
	my ($seq_count,$annots_count,$seq_no_annot,$annots_no_seq_data)=(0,0,0,0);
	
	my $run_but=$frm1->new_button(-text=>"Run",-width=>8, -command=>sub{
			if ((!$drug_targ_seq_file || !$annot_file)|| !$prefix){ Tkx::tk___messageBox(-message => "ERROR: Input file missing", -type=>"ok", -title=>"Alert", -icon=>'warning');  $crt_win->g_destroy();	\&add_a_drug_target_db(); return 0;}
	
		my $BLAST_DB_DIR='KNOWN_DRUG_TARGETS';
		my $database = 'drugTarget_db_names.txt';
		
		my %annot_seq_id;
		my $all_seq=read_fasta_sequence($drug_targ_seq_file); $seq_count = scalar keys %$all_seq;
		open(S,"<$drug_targ_seq_file") or die "$! $drug_targ_seq_file\n";
			while(<S>){	if(/^>(\S+)/g){ $annot_seq_id{$1}=1;}}
		close S;	
		$drg_tar_prg=10;Tkx::update(); 
		open(A,"<$annot_file") or die "$! $annot_file\n";
		
		my @annotations_with_seq_data;
		while(<A>){	
			chomp;
			next if(/^#/ or !$_);
			my $l=$_;
			$annots_count++;
			my @l=split /\t/,$l;
			if($annot_seq_id{$l[0]}){ $annot_seq_id{$l[0]}=$l; push @annotations_with_seq_data,$l; }
			else{$annots_no_seq_data++}


		}
		close A;	
		$drg_tar_prg=20;Tkx::update(); 
		my @seq_id_with_all_annotations;
		foreach my $t(keys %annot_seq_id)
		{ 
			if (!$annot_seq_id{$t})
			{
			$seq_no_annot++;
			}
			else{push @seq_id_with_all_annotations,$t;}
		}
		$drg_tar_prg=50;Tkx::update(); 
		#print STDERR "Summary: Sequences with no annotation: $seq_no_annot\nAnnotations with no sequence data: $annots_no_seq_data\n";

		my $r=fetch_seq_by_id($all_seq,\@seq_id_with_all_annotations);
		write_fasta_seq($r,"./local_dat/$BLAST_DB_DIR/$prefix\_drug_target_db.fasta");
		$drg_tar_prg=55;Tkx::update(); 
		open(F,">> ./local_dat/$BLAST_DB_DIR/$prefix\_targets_annot.txt") or die "$! ./local_dat/$BLAST_DB_DIR/$prefix\_targets_annot.txt";
		$"="\n";
		print F "@annotations_with_seq_data";
		close F;
		$drg_tar_prg=70;Tkx::update(); 
		open(F,">> ./local_dat/$BLAST_DB_DIR/$database") or die "$! ./local_dat/$BLAST_DB_DIR/$database";
		$"="\n";
		print F "$prefix\n";
		close F;
		$drg_tar_prg=90;Tkx::update(); 
		
		system ("formatdb.exe -p T -i .\\local_dat\\$BLAST_DB_DIR/$prefix\_drug_target_db.fasta") if ($blast_version eq 'old_blastall');
		system ("makeblastdb.exe -dbtype prot -in .\\local_dat\\$BLAST_DB_DIR\\$prefix\_drug_target_db.fasta") if ($blast_version eq 'blast+');
		
		
		$drg_tar_prg=100;Tkx::update(); 
		close_tool_window($crt_win,$mw);
		
		reset_params();		
		Tkx::update(); 
		my $d=(scalar @seq_id_with_all_annotations -$annots_no_seq_data)." records imported sucessfully\n";
		Tkx::tk___messageBox(-message => "$prefix added successfully\n\n\nSummary:\n-total input sequence: $seq_count \nSequences with no annotation: $seq_no_annot \nTotal input annotations: $annots_count\nAnnotations with no sequences: $annots_no_seq_data\n--\nImported: $d\n", -type=>"ok", -title=>"Success", -icon=>'info'); 
		
	
	});
	$run_but->g_grid(-column=>0, -row=>6,-padx=>5,-sticky=>"w");
	my $drg_tar_add=$frm1->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$drg_tar_prg);
	$drg_tar_add -> g_grid(-column=>2,-row=>6,-padx=>0,-pady=>1,-columnspan=>2,-sticky=>"w");
	
}

##Meta function
##Args:
##returns:
sub fetch_tax_names
{
	my $level=shift;
	my $sql="\"SELECT DISTINCT $level from taxonomy\"";
	my $l=`sqlite3.exe $broad_spectrum_pathogen_db_sq $sql`;	
	my @l = split /\n/,$l;
	$l = join('} {',@l); $l = " {".$l."} ";
	return ("$l",\@l);
	
}



##ARGS
##returns:
sub add_ecoli_go_db
{

	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Utility: add GO enrichment database");
	#open_tool_window($crt_win,$mw);
	$mw->g_wm_attributes (-disabled  =>0);
	$crt_win->g_wm_attributes (-toolwindow =>0,-topmost =>1);
	$crt_win->g_focus;
	
	
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	
	my ($tax_id,$format,$term_file,$status)= (83333,'uniprot', 'term.txt'.'');
	my (@ids, %GO);
	$frm1->new_ttk__label(-text=>"Reference bacterial taxonomic ID")->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"GO annotation cross-reference database")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>5,-sticky=>"nw");
	$frm1->new_ttk__label(-text=>"Gene Consortitium term file")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>5,-sticky=>"nw");
	
	$frm1 ->new_ttk__entry(-textvariable => \$tax_id,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>0,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$format,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1 ->new_ttk__entry(-textvariable => \$term_file,-width=>30,)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2,-sticky=>"wn");
	$frm1->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$term_file = Tkx::tk___getOpenFile(-parent=>$crt_win,-filetypes =>[['All Files',   '*.txt']]);
	
	})->g_grid(-column=>3,-row=>2,-padx=>2,-pady=>1,-sticky=>"e");
	my $prg_level=0;
	my $run_but=$frm1->new_button(-text=>"Run",-width=>8, -command=>sub{
			$run_but->configure(-state=>'disabled');
			if (!$term_file|| !-e $term_file ){ Tkx::tk___messageBox(-message => "ERROR: Input files missing.Plese refer to manual.", -type=>"ok", -title=>"Alert", -icon=>'warning',-parent=>$crt_win); return 0;}
			if(!check_internet()){Tkx::tk___messageBox(-message => "Active internet required!!. Check internet connection", -type=>"ok", -title=>"Alert", -icon=>'warning',-parent=>$crt_win); return 0;}
			die "No $term_file" if (!$term_file || !-e $term_file);
		my $model_org_fasta= 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625_83333.fasta.gz';
		$status='downloading E.coli proteome'; Tkx::update();
		my $data   = get($model_org_fasta);
		if (!$data) {  $status='download: proteome download failed'; die "UP000000625_83333.fasta.gz download failed\nCheck internet connection";}
		open(OUT, '>./local_dat/GO/UP000000625_83333.fasta') or die "$!";
		binmode OUT;
		print OUT $data;
		close(OUT);
		sleep 1;
		my $fasta_file='./local_dat/GO/UP000000625_83333.fasta';
		$status='BLAST setup'; Tkx::update();
		system ("formatdb.exe -p T -i .\\local_dat\\GO\\UP000000625_83333.fasta") if ($blast_version eq 'old_blastall');
		system ("makeblastdb.exe -dbtype prot -in .\\local_dat\\GO\\UP000000625_83333.fasta") if ($blast_version eq 'blast+');
				
		open (F,"<$fasta_file") or die "$! $fasta_file\n";
		$status="Reading E.coli proteome";Tkx::update();
		if($format eq 'uniprot'){
		#>sp|A5A605|YKFM_ECOLI Uncharacterized protein YkfM OS=Escherichia coli (strain K12) GN=ykfM PE=4 SV=1
			while(<F>)
			{
				chomp;
				if(/^>\w+\|(\S+)\|\w+\s+/){push @ids,$1; }
			}
		}
		close F;
		#print STDERR "$fasta_file Fasta done\nTax id: $tax_id\n";
		#print STDERR "Fetching GO ids from Uniprot\n";
		#print STDERR "\n";
		Tkx::update();
		my $driver   = "SQLite";
		my $dsn = "DBI:$driver:dbname=local_dat/GO/GO.db";
		my $userid = "";
		my $password = "";
		my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 }) or die $DBI::errstr;
		my $stmt = qq(DELETE FROM ecoli_go;);	#
		my $rv = $dbh->do($stmt) or die $DBI::errstr;
		$status="Importing records";Tkx::update();
		foreach my $i (0..$#ids)
		{		
			my $g = fetch_go_uniprot($ids[$i],$tax_id);		
			if($g=='ERROR') {warn "No GO id fround for $ids[$i]. Skip.. \n"; $status="$ids[$i] error";}
			else{		
				my @go = @{$g}; shift @go;
				foreach my $r(@go)
				{
					my $stmt = qq(INSERT INTO ecoli_go (go_id,ecoli_ac) VALUES ('$r','$ids[$i]'););	#
					my $rv = $dbh->do($stmt) or die $DBI::errstr;
					$i++;	
					$prg_level=((($i/scalar@ids)*100)>80?80:(($i/scalar@ids)*100));	Tkx::update();
					$status=sprintf "Importing E.coli record %s/%s",$i,scalar@ids;Tkx::update();
				}		
			}		
		}
		my $stmt = qq(DELETE FROM go_term;);	#Fartry reset
		my $rv = $dbh->do($stmt) or die $DBI::errstr;
		$status='Adding data to go_term';Tkx::update();
		open (I, "<$term_file") or die "$! $term_file";
		my $i=0;
		while(<I>){
		chomp;
		my @l=split /\t/,$_;
		if($l[3]=~m/^GO\:\d{7}/g){
			my $stmt = qq(INSERT INTO go_term (GO_id,ontology_category,term) VALUES (\"$l[3]\",\"$l[2]\",\"$l[1]\"););	
			my $rv = $dbh->do($stmt) or die $DBI::errstr;
			$i++;
		}
		else{warn "No GO id found\n Skipping @l\n"}
		$prg_level=80+(($i/50000)*20);				## assuming approx 50000 recors, extatly 44000 present
		$status=sprintf "Importing GO terms record %s",$i;Tkx::update();
		}
		close I;
		$status='DONE';Tkx::update();
		$prg_level=100;
		#close_tool_window($crt_win,$mw);
		$crt_win->g_destroy();
		$mw->g_wm_attributes (-disabled  =>0);
		$mw->g_focus;
		Tkx::tk___messageBox(-message => "GO meta-data imported succesfully", -type=>"ok", -title=>"Success", -icon=>'info');
	});		
	
	$run_but->g_grid(-column=>0,-row=>8,-padx=>2,-pady=>1,-sticky=>"wn");
	$frm1->new_ttk__label(-text=>"Progress")->g_grid(-column=>2,-row=>8,-padx=>2,-pady=>5,-sticky=>"ne");
	my $go_data_add=$frm1->new_ttk__progressbar(-orient => 'horizontal', -length => 100, -mode => 'determinate', -variable=>\$prg_level);
	$go_data_add -> g_grid(-column=>3,-row=>8,-padx=>3,-pady=>5,-sticky=>"w",-columnspan=>2);
	$frm1->new_ttk__label(-textvariable=>\$status)->g_grid(-column=>0,-row=>9,-padx=>2,-pady=>5,-sticky=>"w",-columnspan=>2);
	$status=(check_internet()?"Active internet": 'No internet'); Tkx::update();
}

##Meta function
#fetch_go_uniprot('P10408','83333');
##Args: protein uniprot_id, org_id (Ecoli)
## returns: array ref ([uniprot_id,GO id1, GO_id2, what ever comes in fetch ])
sub fetch_go_uniprot
{
my $uniprot_id=shift;
my $org_id =shift;
	my $url ='http://www.uniprot.org/uniprot/?query=accession:'.$uniprot_id.'+AND+organism:'.$org_id.'&format=tab&columns=id,go-id';
	my $content = get $url;
#	die "Couldn't get $url" unless defined $content;
	if(!$content || $content =~m/^400/g) {
		print STDERR "Bad request. There is a problem with your input.\n";
		return 'ERROR';
	  } 
	  else {
		#print "$content\n";
		my ($t,$go)=split /\n/,$content;
		my @g = split /;\s*/,$go;
		return \@g;
	}	
}

##args:	ref of array contain tax level names, L path
##returns: an scalar containting blast db naemss(windows formated)
sub create_broad_spe_db_array
{
	my $array=shift;
	my $path=shift;
	my $db_list;
	$path=~ s{/}{\\}g; 
	my $range = join ('\",\"',@$array); $range='(\"'.$range.'\")';
	my $e= ($taxon_id? " AND TaxID != $taxon_id":"");
	my $sql="\"SELECT fasta_file from taxonomy WHERE $tax_level IN $range $e\"";
	my $l=`sqlite3.exe $broad_spectrum_pathogen_db_sq $sql`;	
	my @l = split /\n/,$l;
	
	if($#l<0){
	print STDERR "ERROR: No pathogen proteome found while reading data for Broad-spectrum analysis\n. Add Few using Utility menu\n";
	$setup_error.="*** No pathogen proteome found while reading data for Broad-spectrum analysis.\nAdd Few using Utility menu\n\n\n";
	}
	
	my @db;
	foreach(@l)
	{
		#my $r='\"'.$path.'\\'.$_.'\"';
		my $r='"'.$path.'\\'.$_.'"';
		push @db,$r;
	}
	
	#$db_list=join(" ",@db); #$db_list.="-"; $db_list=~s/\"-//g;
	$db_list=join(",",@db);
	#$db_list='"'.$db_list.'"';
	#print "$db_list\n";
	return $db_list;
}

##ARGS:
##returns
sub chk_GO_db
{
	my $driver   = "SQLite";
	my $dsn = "DBI:$driver:dbname=$GO_db";
	my $userid = "";
	my $password = "";
	my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 })
	or die $DBI::errstr;

	my $stmt = qq(SELECT COUNT(*) FROM go_term;);	#
	my $sth = $dbh->prepare($stmt);
	$sth->execute(); my $r = $sth->fetch()->[0];
	if($r<1){
	print STDERR "ERROR: No GO annotation found while reading data for go_term table\n. Add Few using Utility menu\n";
	$setup_error.="*** No GO annotation found while reading data \nfor go_term table. Add Few using Utility menu\n\n\n";
	
	}
	
	my $stmt = qq(SELECT  COUNT(*) FROM ecoli_go;);	#
	my $sth = $dbh->prepare($stmt);
	$sth->execute(); my $r = $sth->fetch()->[0];
	if($r <1){
	print STDERR "ERROR: No E.coli GO annotation found while reading data for ecoli_go table \n. Add Few using Utility menu\n";
	$setup_error.="*** No E.coli GO annotation found \nwhile reading data for ecoli_go table. Add using Utility menu\n\n\n";
	return 0;
	}
	else{return 1}
}

##ARgs: go_id
##Returns: ref to a hash
#ecoli_go  go_term
#SELECT go_term.ontology_category, go_term.term. ecoli_go.go_id FROM go_term LEFT JOIN ecoli_go ON #go_term.go_id=ecoli_go.go_id WHERE ecoli_go.ecoli_ac = "P0A6A8";
#molecular_function|acyl binding|GO:0000035
#molecular_function|ACP phosphopantetheine attachment site binding involved in fatty acid biosynthetic process|GO:0000036
#cellular_component|cytosol|GO:0005829
#cellular_component|cytoplasm|GO:0005737
#biological_process|fatty acid biosynthetic process|GO:0006633
#biological_process|lipid biosynthetic process|GO:0008610
#molecular_function|phosphopantetheine binding|GO:0031177
sub fetch_GO_db
{
	my $ecoli_ac = shift;
	
	my $driver   = "SQLite";
	my $dsn = "DBI:$driver:dbname=$GO_db";
	my $userid = "";
	my $password = "";
	my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 })
	or die $DBI::errstr;

	my $stmt = qq(SELECT go_term.ontology_category, go_term.term, ecoli_go.go_id FROM go_term LEFT JOIN ecoli_go ON go_term.go_id=ecoli_go.go_id WHERE ecoli_go.ecoli_ac = \"$ecoli_ac\";);	#
	my $sth = $dbh->prepare($stmt);
	$sth->execute(); 
	my $go_terms;			## ref of hash =
		#{
		#	-molecular_function=>[], ##acyl binding|GO:0000035, phosphopantetheine binding|GO:0031177
		#	-cellular_component=>[],
		#	-biological_process=>[],
		#}	
	my @row;
	while ( @row = $sth->fetchrow_array()) {
		$go_terms->{$row[0]}=[] if !$go_terms->{$row[0]}; 
		push @{$go_terms->{$row[0]}},$row[1].'|'.$row[2];
	}	
	$sth->finish();
	$dbh->disconnect();	
	return $go_terms;		
}





########################################################################################################
#Help menu functions
#########################################################################################################

##Menu-> Help-> Manual
sub manual
{
 system ("start \"\" /max \"User_Manual.pdf\"");
}

sub citation
{

Tkx::tk___messageBox(-type => "ok",
					-message => $citation_text,
					-title => "Citation");

}

##Menu-> Help-> about
sub about
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("About");
	open_tool_window($crt_win,$mw);
	my $heading = $crt_win->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool\nVersion: 1.2.2b\nLast update: $last_update\n",-justify=>"center",-foreground=>"blue");
	$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>10);
	
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>1,-column=>0,-sticky=>"nsew");
	my $t=$frm1->new_tk__text(-width => 40, -height => 10,-wrap=>'word',);
	$t->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>1, -columnspan=>2);
	$t->insert("end", $about_text);
	
	my $hscroll = $frm1->new_tk__scrollbar(-orient => "horizontal", -command => [$t, "xview"]);
	my $vscroll = $frm1->new_tk__scrollbar(-orient => "vertical", -command => [$t, "yview"]);
	$hscroll->g_grid(-column => 0, -row => 3, -sticky => "we",-columnspan=>4);
	$vscroll->g_grid(-column => 4, -row => 2, -sticky => "ns");
	$frm1->new_ttk__sizegrip()->g_grid(-column => 4, -row => 3, -sticky => "se");
	$t->configure(-yscrollcommand => [$vscroll, "set"], -xscrollcommand => [$hscroll, "set"]);
	my $ok=$frm1->new_button(-text=>"Close",-width=>10,-command=>sub 
		{
		$crt_win->g_destroy(); close_tool_window($crt_win,$mw);
		})->g_grid(-column=>0,-row=>6,-padx=>1,-pady=>5,-sticky=>"ne");	
	$frm1->new_ttk__label(-text=>"Report bugs to:\nkcm.eid[at]gmail.com",-foreground=>"red")->g_grid(-column=>0,-row=>8,-padx=>1,-pady=>9,-sticky=>"w");	
}

##Args:
##Returns:
sub welcome_message
{
	Tkx::tk___messageBox(-title=>"Welcome to EDTI",-message => "
  EXOGENOUS DRUG TARGET IDENTIFICATION TOOL
   
EDTI is a pipeline to identify and analyse 
putative bacterial drug-targets.

Use File menu (or Ctrl+n) to create a new 
project.

Use Downstream-analysis menu to analyse 
previous run results.
	
Use Utilities to add/update metadata used 
by the tool.

Use the Run button below to run the program.
	
Presee F1 for help.", -parent=>$mw,-icon=>"info");
}



################################################################################################
#  Meta functions
################################################################################################


##ARgs: blast tab out file, percenta)identy, format
##Returns: ref of hash; [ ]
#sp|A8FJG5|DNAA_CAMJ8	sp|P03004|DNAA_ECOLI	37.65	332	203	2	101	430	132	461	3e-063	 235
sub process_GO_BLAST_out
{
	my $file= shift;
	my $identity = shift || 30;
	my $format =shift || 'uniprot';
	my %h;
	
	open (O, "< $file") or die "$! ";
	while(<O>)
	{
		chomp;
		my @l = split /\s+/,$_;
		if ($l[2]>=$identity)
		{
			my ($a,$ac,$t);
			($a,$ac,$t)= split /\|/,$l[1] if $format eq 'uniprot';
			$ac= $l[1] if $format ne 'uniprot';
			$h{$l[0]}=$ac; 
		}		
	}
	close O;
	return \%h;	
}


##args: nothing
##returns: an ref of hash; ACIAD	Acinetobacter sp. (strain ADP1) 
sub read_broad_spe_codes
{
	
	my %h;
	my $sql="\"SELECT org_code,species from taxonomy\"";
	my $l=`sqlite3.exe $broad_spectrum_pathogen_db_sq $sql`;	
	my @l = split /\n/,$l;
	foreach(@l)
	{
		chomp;
		my @l=split /\|/,$_;
		$h{$l[0]}=$l[1];
	}
	
 return \%h;
}

##args:dat file
##returns: an scalar containing database names: ' {All} {database1} {database2} '
sub read_drugTarget_db{
	my $file=shift;
	my $db=" ";
	open(F,"<$file") or die "$! $file\n  ";
	while(<F>)
	{
		chomp;
		next if(/^#/ or !$_);
		$db.=" {$_}";
	}
	close F;
	return $db;
}

##args:ref of array of db names and database path
##returns: an scalar containing windows BLASTp database(S) input
sub create_drugTarget_blast_db{
	my $arr=shift;
	my $path=shift;
	my $db_list;
	$path=~ s{/}{\\}g; #$root_path='"'.$root_path.'"'; 
	my @db;
	#open(F,"<$file") or die "$! $file\n  ";
	foreach(@{$arr})
	{
		chomp;
		next if(/^#/ or !$_);
		my $r='\"'.$path.'\\'."$_\_drug_target_db.fasta".'\"';
		push @db,$r;
	}
	#close F;
	
	$db_list=join(" ",@db); #$db_list.="-"; $db_list=~s/\"-//g;
	$db_list='"'.$db_list.'"';
	#print "$db_list\n";
	return $db_list;

}

##args:dat file
##returns: an refrence scalar to ana hash; {id=> [annotation array ref]}
sub read_drugTarget_annot{
	#my $file=shift;
	my $path=shift;
	my @annotation_files=<"$path/*_targets_annot.txt">;
	my %hash;
	foreach my $f(@annotation_files){
		open(F,"<$f") or die "$! $f \n  ";
		while(<F>)
			{
				chomp;
				next if(/^#/ or !$_);
				my @l=split /\t/,$_;
				#sp|Q8GBW6|12S_PROFR	Methylmalonyl-CoA carboxyltransferase 12S subunit	DB04045; DB04183	Propionibacterium freudenreichii subsp. shermanii	Q8GBW6
				$hash{$l[0]}=\@l;			
			}
	}
return \%hash;
}

##args: *.fasta
##returns: ref_address of an hash address -id=>{-full_id=>, -seq=>}
sub read_fasta_sequence
{
	my $file=shift;
	my %f;
	open (F, "<$file") or die "$file $!";
	my $id;
	while(<F>)
	{
	chomp;
		if(/^>(\S+)/){$id=$1;die "ERR:$id repeating in fasta file.\n FASTA file is not uniuq\n "if $f{$id};  $f{$id}={-desc=>$_, -seq=>"", -len=>0};   }
		else{$f{$id}->{-seq}.=$_;  $f{$id}->{-len}=length($f{$id}->{-seq}); }	
	
	}
	close F;
	return \%f;
}

##args: address of an hash address {-id=>, -full_id=>, -seq=>}, out file*.fasta
##returns:1 on success 
sub write_fasta_seq
{
	my ($f,$file)=@_;
	open (F, ">$file") or die "$file $!";
	foreach my $id(sort keys %$f)
	{
		print F "$f->{$id}->{-desc}\   |$f->{$id}->{-len}\n$f->{$id}->{-seq}\n";	##spaces added ; BUG fixed
	}

	 close F;
}


##args: address of an hash address {-id=>, -full_id=>, -seq=>}, ref array of ids;
##returns:address of an hash address {-id=>, -full_id=>, -seq=>} of selected ids;
sub fetch_seq_by_id
{
	my ($f,$id_ref)=@_;
	my %r;
	foreach my $i(@$id_ref)
	{
		$r{$i}->{-desc}=$f->{$i}->{-desc};
		$r{$i}->{-seq}=$f->{$i}->{-seq};
		$r{$i}->{-len}=$f->{$i}->{-len};
	}
	return \%r;
}



##args: fasta seq file/Unix path
##returns: number of lines with >
sub count_fasta_seq
{
	my $seq_file=shift;
	my $count=0;
	open (F, "$seq_file") or die "$seq_file $!";
	while(<F>)
	{
		$count++ if /^>\S+/;
	}
	close F;
	return $count;
}

##args: fasta seq file
##returns: returns >(\S+) ids in the fasta file as array ref;only the text after > and before any space
sub ids_in_fasta_seq
{
	my $seq_file=shift;
	my @id;
	open (F, "$seq_file") or die "$seq_file $!";
	while(<F>)
	{
		if (/^>(\S+)\s+/){push @id,$1;}
	}
	close F;
	return \@id;
}


##args: cdhit clstroutput
##returns: reference of two arrays;\@paralogs,\@uniq
sub process_cdHit_clstr			##Correct method, not inclusidn the fist hit of a cluster
{
my $file=shift;
my (@paralogs,@uniq);

open (P, "<$file") or die "$!$file";
my $clstr=0;
my @a;
while(<P>){
		if(/^>Cluster\s+\d+/){  if(scalar @a==1){push @uniq,@a; } elsif(scalar @a>1){ push @paralogs,@a; }  undef @a;}
		elsif(/\d+\s*\d+aa,\s*>(\S+)\.\.\.\s+/){push @a,$1; }
		else {warn "err\n"}
		
	}
	if(scalar @a==1){push @uniq,@a; } elsif(scalar @a>1){ push @paralogs,@a; }	##process last line
	return (\@paralogs,\@uniq);
}

=inp
>Cluster 88
0	662aa, >tr|A8FP24|A8FP24_CAMJ8... *
>Cluster 89
0	659aa, >tr|A8FJV0|A8FJV0_CAMJ8... *
1	651aa, >tr|A8FNM2|A8FNM2_CAMJ8... at 72%
>Cluster 90
0	659aa, >tr|A8FMX5|A8FMX5_CAMJ8... *
parse the lines with ^0 as uniq and other than 0 as paralog
=cut

##args: fasta_file and blast tabular out file and perc_identity_cutoff
##returns: ref of two arrays: 1. blast_hits and non_blast_hits respectively
sub process_host_blast_out
{
	my $input_seq=shift;
	my $file=shift;	##blast_out
	my $perc_identity=shift || 0;
	
	my (@blast_hits,@non_blast_hits);
	open (P, "<$file") or die "$! $file";
	my %a;
	while(<P>){
		my @l=split /\s+/,$_;
		$a{$l[0]}=1 if $l[2] >=$perc_identity ;	##keep if identity is grter than #% identity
	}
	close P;
	@blast_hits= keys %a;
	
	my @b;
	open (P,$input_seq) or die "$!$input_seq";
	while(<P>){
	if(/^>(\S+)/){my $t=$1; $t=~s/\|/\\|/g;  if(!(grep{/^$t$/}@blast_hits)){ $t=~s/\\\|/|/g; push @non_blast_hits,$t; }      }
	#else{}
	}
	close P;
	return (\@blast_hits,\@non_blast_hits);
}


##args:blast_out, identity_cutoff
##returns:address of an has; query_id=> [org1, org2, org3, ]
sub process_broad_spe_BLAST_out
{
	my $file=shift;
	my $homology_cutoff=shift;
#sp|Q46170|ARCD_CLOPE	tr|Q6FCC9|Q6FCC9_ACIAD	44.50	436	229	2	9	439	4	431	2e-081	 297
#sp|Q46170|ARCD_CLOPE	tr|B2HV75|B2HV75_ACIBC	42.63	441	240	2	10	445	5	437	7e-076	 278
	open (F,"<$file") or die ("$! $file");
	my %hash;
	while(<F>)
	{
		chomp;
		next if !$_;
		my @l=split /\s+/,$_;
		$hash{$l[0]}=[] if !$hash{$l[0]};
		next if ($l[2]<$homology_cutoff);
		#my $prot_org;
		my ($x,$y,$prot_org)=split (/\|/,$l[1]);
		my($prot,$org)=split /_/,$prot_org;
		push @{$hash{$l[0]}}, $org if !(grep{/^$org$/}@{$hash{$l[0]}});
	}
	close F;

return \%hash;
}


##args: blast_out filename; and total seq count
#returns:number of lines in it;
sub blast_progress
{
 my $out1=shift;	## path style must be in windows format
 my $out2=shift;	## path style must be in unix format; to check existance
 my $total_seq=shift;
 my $temp_file=$out1.'.tmp';
 copy (unix_path($out1), unix_path($temp_file));	##just to tacle file lock; may increase program time

 
 if(-e $out2){
	open (O, unix_path($temp_file)) or die "$! $temp_file\n";
	my $processed_seq=	@{[<O>]};										#int `findstr /R /N "^" $temp_file | find /C ":"` ; 
	close O;
	return ($processed_seq/$total_seq)*100;
	unlink 	unix_path($temp_file);																#`del $temp_file`;
	}
else{return 0}
	
}

##args:esential_protein_like_aa_blast_out, bit_score_cutoff;
##returns:
sub filter_by_bit_score
{
 my $file=shift;
 my $bit_score=shift;
 my (@blast_hits,@non_blast_hits);
	my %a; my %b;
	open (P, "<$file") or die "$! $file";
	while(<P>){
		my @l=split /\s+/,$_;
		$a{$l[0]}=1 if $l[11]>=$bit_score;	
		$b{$l[0]}=1 if $l[11]<$bit_score;	
	}
	close P;
	@blast_hits= keys %a;
	@non_blast_hits= keys %b;
	return (\@blast_hits,\@non_blast_hits);
}

###Not used; later while reading exisitng project
##Args: a tab deleim file with three # line at the top
##Return: a hash of hashes;
sub process_centrality_measure_file
{
 my $file=shift;
 open (O,$file) or die"$! $file";
 my %h;
 while(<O>){
		next if /^#/;
		my @l =split /\s+/,$_;
		$h{$l[0]}={-PPI_database_id=>$l[1],
		-degree	=>$l[2],
		-radiality	=>$l[3],
		-closeness	=>$l[4],
		-eccentricity	=>$l[5],
		-normalized_total_score=>$l[6],
		};
	}
  close O;
return \%h;
}


##Args: optionally 'L'/'W': default is windows format
##returns: path of My Documents\Desktop; OS specific; windows; 'echo %USERPROFILE%'
sub get_my_document_path
{
my $os=shift || 'W';
my $p=`echo %USERPROFILE%`;
chomp($p);
#$p =~ s{\\}{\\\\}g;
$p.='\\Desktop';		##in users desktop;

if(uc($os) eq 'L'){$p=~s/\\/\//g; return $p;}
else{$p = '"'.$p.'"'; return "$p";}
}

##Args: optionally 'L'/'W': default is windows format
##returns: path of C:\windows OS specific; windows; 'echo %SystemRoot%'
sub get_windows_root_path
{
my $os=shift || 'W';
my $p=`echo %SystemRoot%`;
chomp($p);
if(uc($os) eq 'L'){$p=~s/\\/\//g; return $p;}
else{$p = '"'.$p.'"'; return "$p";}
}



##Args: Win formated path; STRING ("C:\\Users\\SGPGI.SGPGI-PC\\KANHU\\EDIT\\Prev_data";)
##returns: Unix formated path string (C:/Users/SGPGI.SGPGI-PC/KANHU/EDIT/Prev_data)
sub unix_path
{
	my $p = shift;
	#chomp ($p);		##
	$p=~s/\\/\//g;		## \\ to /
	$p=~s/"//g;
	return "$p";		##
}

##Args: Unix formated path; STRING ; C:/Users/SGPGI.SGPGI-PC/KANHU/EDIT/Prev_data
##returns: Win formated path string ("C:\\Users\\SGPGI.SGPGI-PC\\KANHU\\EDIT\\Prev_data")
sub win_path
{
	my $p = shift;
	chomp ($p);			##
	#$p = ~s/"{2,}/"/g;		## sanity checking
	$p=~s/\//\\/g;		##  /to \\
	return '"'.$p.'"';			##
}


##ARgs:
##return:0/1'
sub check_internet
{
my $status=0;
my $ping=`ping -n 1 www.microsoft.com`;
chomp($ping);
my $status = scalar (split /\n/,$ping)-1;
#print STDERR $status;
return $status;
}

##Args:name of executable file
##Retrns: 0/1
sub check_executable_on_PATH
{
	my $p =shift;
	my $r = `$p`; chomp($r);
	if(!$r){return 0;}
	else{return 1}
}

sub count_blast_hits
{
	my $blast_out=shift;
	my $total_lines=0;
	open (O, "$blast_out") or die "$! $$blast_out\n";
	 $total_lines =@{[<O>]};
	close O;
	return $total_lines;

}

##ARGS: two fasta file;file-A file-B; file-A will be checked against B
##return : percentage of ids from A present in B (float); AND a ref to array of ids common in both file
sub compare_two_fasta_file{
	my ($fas_A,$fas_B)=@_;
	my $a = read_fasta_sequence($fas_A);
	my $b = read_fasta_sequence($fas_B);
	my($a_mat,$a_tot)=(0,0);
	my @m;
	foreach my $i(keys %$a)
	{
		$a_tot++;
		if ($b->{$i}){
			$a_mat++;
			push @m,$i;
		}
	}
	undef $a; undef $b;	
	my $m = sprintf("%4.1f",($a_mat/$a_tot)*100);
	return ($m,\@m);
}

##Args: ref of a toplevel, ref of parent toplevel window
##returns: none
sub open_tool_window
{
 my ($crt_win,$mw)=@_;
 $mw->g_wm_attributes (-disabled  =>1);
 $crt_win->g_bind("<Escape>", sub{close_tool_window($crt_win,$mw)});
 $crt_win->g_wm_attributes (-toolwindow =>1,-topmost =>1);
 $crt_win->g_wm_protocol( WM_DELETE_WINDOW =>sub{close_tool_window($crt_win,$mw);});
 $crt_win->g_focus;
 }

##Args: ref of a toplevel, ref of parent toplevel window
##returns: none
sub close_tool_window
{
 my ($crt_win,$mw)=@_;
$crt_win->g_destroy;$mw->g_wm_attributes (-disabled  =>0); $mw->g_raise;$mw->g_focus;
}


Tkx::MainLoop();
