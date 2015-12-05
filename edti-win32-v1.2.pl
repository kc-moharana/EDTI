use strict;
use Tkx;
use GD::Graph::bars;
Tkx::package_require("tile");
Tkx::package_require("BWidget");
#Tkx::package_require('tooltip');
#Tkx::namespace_import("::tooltip::tooltip");
use Cwd;
use DBI;
use Graph;
use Graph::Undirected;
use File::Copy qw(copy);
#use LWP::Simple;



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
our $taxon_id;			## NCBI taxon id; trash

our $project_name="New_Project";	
our $root_path=	get_my_document_path();		#getcwd();##update later
$root_path=~s/\//\\\\/g;					##Win formated
#our $prog_install_path=getcwd();			
our $L_root_path = get_my_document_path('L');	## keeping an extra variable; storing path in Unix format
our $installation_path=getcwd();			##dont update; essential data files and folders;
#$installation_path=~s/\//\\\\/g;
our $last_update="25 June 2015";

our $front_page_status="File > Create a new project. ";
our $filter_param_settings_file;						#$root_path."/param.txt";
print STDERR "Intallation path:$installation_path\n";
our $read_seq_prg=0;
our $cdhit_prg=0;
our $remove_short_seq_prg=0;
our $string_srch_prg=0;
our %sequence_summary=(
	-total_seq=>0,
	-very_short_seq=>0,
	-paralogslogs=>0,
	-host_orthologs=>0,
	-drug_target_homologs=>0,
	-putative_drug_targets=>0,
	
);


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

													#
my $CPU_list;   
foreach(1..$avl_CPU_count){$CPU_list.=$_." ";   } 
my $use_cores=$avl_CPU_count;		##using max no of CPUs by default

our $min_aa_len=50; #after translations
#our $min_ntd_len=100;
our $cd_hit_identity=60;
our $chk_cdhit=1;
print STDERR "Reading drug Target data:";
our $drug_db_names=read_drugTarget_db("$installation_path/local_dat/drugTarget_db_names.txt");#' {All} {drugBank} {PTTD} ';	##read files to update it;
our $ref_drug_db_array=[];
open(G,"$installation_path/local_dat/drugTarget_db_names.txt") or die"$!$installation_path/local_dat/drugTarget_db_names.txt"; while(<G>){chomp; push @$ref_drug_db_array,$_;};close G;
our $drug_blast_db_names=create_drugTarget_blast_db($ref_drug_db_array,"$installation_path/local_dat/KNOWN_DRUG_TARGETS");
our $drug_target_annot=read_drugTarget_annot("$installation_path/local_dat/KNOWN_DRUG_TARGETS");
print STDERR "DONE";
our $enable_broad_spe=1;			##decrypted
print STDERR "\nReading broad-spectrum data:";
our $broad_spectrum_pathogen_db_sq =$installation_path."\\local_dat\\pathogen_taxonomy.db";
$broad_spectrum_pathogen_db_sq=~s/\//\\/g;
##'"\"COMPLETE_PTM_PATHPGENS\" "COMPLETE_PTM_PATHPGENS\ACIBC" "COMPLETE_PTM_PATHPGENS\ACIBS" "COMPLETE_PTM_PATHPGENS\ACIBT" "COMPLETE_PTM_PATHPGENS\BURM1\""';##read files to update it;
our $broad_spe_species_per_query=0;
our $broad_spe_org_code_full_name=read_broad_spe_codes("$installation_path/local_dat/BROADSPECTRUM_CODE_ANNOTATION.txt");#Stores hash ref
our $tax_level="Family";
our $brd_sp_db_levels;
our $ref_brd_sel_db_array=[]; ## array ref ; array saves  tax_levels sel by user			
($brd_sp_db_levels,$ref_brd_sel_db_array) = fetch_tax_names($tax_level);		
our $broad_spectrum_pathogen_db_list=create_broad_spe_db_array($ref_brd_sel_db_array,"$installation_path/local_dat/COMPLETE_PTM_PATHPGENS");
print STDERR "DONE\n";

our $database = $installation_path."\\local_dat\\PPI\\PPI_sqlite3.db";
	$database =~s/\//\\/g;

our $blast_prg1=0; 
our $blast_prg2=0;
our $blast_prg4=0;		## drugBank progress
our $blast_prg3=0;


our ($blast_path,$cdhit_path,$formatdb_path,$sqlite_path)=("$installation_path/executables/blastall.exe","$installation_path/executables/cd-hit.exe","$installation_path/executables/formatdb.exe","$installation_path/executables/sqlite3.exe");


our ($e_val_1,$out_fmt_1,$sub_matrix_1,$gap_score_1, $extra_params_BLAST1,$word_size_1,$threshold_1,$perc_identity_1)=(0.01,8,"BLOSUM62","11,1","-b 1",3,11,20);	#-b 1

our ($e_val_2,$out_fmt_2,$sub_matrix_2,$gap_score_2, $extra_params_BLAST2,$word_size_2,$threshold_2,$perc_identity_2)=(0.0000000001,8,"BLOSUM62","11,1","-b 1",3,11,0);

our ($e_val_3,$out_fmt_3,$sub_matrix_3,$gap_score_3, $extra_params_BLAST3,$word_size_3,$threshold_3,$perc_identity_3)=(0.01,8,"BLOSUM62","11,1","",3,11,30);
our $broad_spe_BLAST_cutoff=40;		##percentidentity
our ($e_val_4,$out_fmt_4,$sub_matrix_4,$gap_score_4, $extra_params_BLAST4,$word_size_4,$threshold_4,$perc_identity_4)=(0.01,8,"BLOSUM62","11,1","-b 1",3,11,0);

our $PPI_score_cutoff=700;		##STRING score cutoff
our $top_hub_perc=20;

our	$brd_spec_chk=1;
our	$drg_tar_blast_chk=1;
our	$KEGG_chk=0;
our	$sub_cel_loc_chk=0;



my $about_text="The availability of complete genome sequences of pathogenic bacteria and their protein complements in public domain has made it possible to determine potential drug targets in these pathogens using computer-based in-silico techniques. Intersection of two datasets, namely \n\t(i)a pathogen's subtractive proteome dataset with the host proteome, and \n\t(ii) the pathogen's minimal essential protein dataset, should represent a set of proteins whose manipulation may reasonably be expected to interfere with the pathogen's survival without adversely affecting the host. These proteins could thus act as potential targets for drugs acting against the particular pathogen.\n\nThis program comes with ABSOLUTELY NO WARRANTY; for details http://www.gnu.org/licenses/.
This is free software, and you are welcome to redistribute it under certain conditions.\n\n";


my $citation_text="Sarangi AN et al., Exogeneous drug target identifiation tool.\nPMID:000000";


if($filter_param_settings_file){
	print STDERR "Loading Parameter file....";
	
	open (P, "$filter_param_settings_file") or die "$!";
	while(<P>)
	{
		if(/^CPU\s+=\s+(\S+)/){  $use_cores=$1}
		if(/^MIN_AA_LEN\s+=\s+(\S+)/){ $min_aa_len=$1}
		if(/^CHK_CD_HIT\s+=\s+(\S+)/){ $chk_cdhit =$1}
		if(/^CD_HIT_IDN\s+=\s+(\S+)/){ $cd_hit_identity=$1}
		#if(/^BRD_SPEC_THR\s+=\s+(\S+)/){ $broad_spe_BLAST_cutoff=$1}
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
       
		if(/^E_VAL_4\s+=\s+(\S+)/){ $e_val_1=$1}
		if(/^OUT_FMT_4\s+=\s+(\S+)/){ $out_fmt_1=$1}
		if(/^SUB_MAT_4\s+=\s+(\S+)/){ $sub_matrix_1=$1}
		if(/^GAP_SCOR_4\s+=\s+(\S+)/){ $gap_score_1=$1}
		if(/^EXTRA_PARAM_4\s+=\s+(\S+)/){ $extra_params_BLAST1=$1}
		if(/^WORD_SIZE_4\s+=\s+(\S+)/){ $word_size_1=$1}
		if(/^THRHOLD_4\s+=\s+(\S+)/){ $threshold_1  =$1}
		if(/^PERC_IDENTY_4\s+=\s+(\S+)/){ $perc_identity_4=$1}
 
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
#my $utils = $menu->new_menu;
my $dwn_str_anal=$menu->new_menu;
my $help = $menu->new_menu;

$menu->add_cascade(-menu => $file, -label => "File",-underline=>0);
$menu->add_cascade(-menu => $settings, -label => "Settings",-underline=>0);
$menu->add_cascade(-menu => $dwn_str_anal, -label => "Downstream analysis",-underline=>0);
$menu->add_cascade(-menu => $help, -label => "Help",-underline=>0,);
my $system = Tkx::widget->new(Tkx::menu($menu->_mpath . ".system"));
$menu->add_cascade(-menu => $system);

if (Tkx::tk_windowingsystem() eq "aqua") {
		$mw->g_bind("<2>", [sub {my($x,$y) = @_; $menu->g_tk___popup($x,$y)}, Tkx::Ev("%X", "%Y")] );
		$mw->g_bind("<Control-1>", [sub {my($x,$y) = @_; $menu->g_tk___popup($x,$y)}, Tkx::Ev("%X", "%Y")]);
	} 
	else {
		$mw->g_bind("<3>", [sub {my($x,$y) = @_; $menu->g_tk___popup($x,$y)}, Tkx::Ev("%X", "%Y")]);
	}
#add menu items 
$file->add_command(-label => "Create new Project",-underline=>1, -command => \&create_project);
#$file->add_command(-label => "Open Project",-underline=>0, -command => sub {});	
#$file->add_command(-label => "Save Project",-underline=>0, -command => sub {});	
$file->add_command(-label => "Quit",-underline=>0, -command => sub {system("del $root_path\\*.tmp"); exit;});

$settings->add_command(-label => "Open settings",-underline=>1, -command =>\&settings);
$settings->add_command(-label => "Down-stream analysis",-underline=>1, -command =>\&down_str_anal);
$settings->add_command(-label => "Executable files path",-underline=>1, -command => \&executable_paths);

$dwn_str_anal->add_command(-label =>"Broadspectrum analysis", -underline=>1, -command =>sub {});	## functions updated later
$dwn_str_anal->add_command(-label =>"Compare with known targets", -underline=>1, -command =>sub {});
$dwn_str_anal->add_command(-label =>"GO analysis",-underline=>1, -command =>sub {});
$dwn_str_anal->add_command(-label =>"Sub-cellular localization",-underline=>1, -command =>sub {});


$help->add_command(-label => "Manual",-underline=>0, -command =>\&manual);	
$help->add_command(-label => "About",-underline=>0, -command =>	\&about);
$help->add_command(-label => "Citation",-underline=>0, -command =>	\&citation);

my $frnt=$mw->new_ttk__frame(-borderwidth=>12,-relief => "sunken",-width => 600, -height => 800,-padding => "0 5 5 5" );
$frnt->g_grid(-column=>0,-row=>0,-sticky => "nwes");

my $frnt_top=$frnt->new_ttk__frame(-borderwidth=>0, -width => 600, -height => 500);
$frnt_top->g_grid(-column=>0,-row=>0,-sticky=>"nswe");

 Tkx::image_create_photo( "BANER", -file => "banner.gif");
($frnt_top->new_ttk__label(-image=>'BANER'))->g_grid(-column=>0, -row=>0,-sticky=>"nwes",-columnspan=>2, -padx=>60);		

 my $heading = $frnt_top->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool",-justify=>"center",-foreground=>"blue",-font => "Helvetica 16 bold underline");
$heading->g_grid(-column=>0,-row=>0,-sticky=>"s",-padx=>50);
 
my $message = $frnt_top->new_ttk__label(-textvariable =>\$front_page_status,-justify=>"left",-foreground=>"red",-font => "Helvetica 12 italic");
$message->g_grid(-column=>0,-row=>2,-sticky=>"wn",-padx=>50,-pady=>10,-rowspan=>10);


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

PPI_THR = $PPI_score_cutoff
TOP_HUB_PERC = $top_hub_perc
 
 ";
Tkx::tk___messageBox(-message => "Parameters used in the current analysis \nwere saved to $L_root_path/Parameters.txt");
 close P;
 });
$save_options_but->g_grid(-column=>2, -row=>0,-sticky=>"e", -padx=>10);

$dwn_str_anal->entryconfigure("Broadspectrum analysis", -command =>sub{ $frnt_top->g_destroy; $frnt->configure(-padding => "0 0 0 0"); broad_spect_run(\$frnt,\$run_but); }); 
$dwn_str_anal->entryconfigure("Compare with known targets", -command =>sub{ $frnt_top->g_destroy; $frnt->configure(-padding => "0 0 0 0"); comp_known_DT(\$frnt,\$run_but); }); 

Tkx::update();

##Args: ref of front panel Frams, ref of Run button, ref of an global variable to store ref of all seq data;
##Returns: nothing
sub main_script
{
	my $frm=shift;
	my $run_button=shift;
	my  $all_t_seq;				##stores all seq
	
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

	$all_t_seq=read_fasta_sequence($Tproteome_file);
	open(R1,"> $L_root_path/excluded_seq_step-1.fasta") or die "$L_root_path/excluded_seq_step-1.fasta";
	open(A1,"> $L_root_path/accepted_seq_step-1.fasta") or die "$L_root_path/accepted_seq_step-1.fasta";
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

	$cdhit_path=~ s{/}{\\}g; $cdhit_path='"'.$cdhit_path.'"';

	if($chk_cdhit)
		{
			my $cdHit_c=$cd_hit_identity/100;

			unlink "$L_root_path/cdHit_out.clstr.txt";		
			#`$cdhit_path -i $root_path/accepted_seq_step-1.fasta -o $root_path/cdHit_out -d 0 -c $cdHit_c -n 3`;
			`echo echo off > batch.bat`;
			`echo color B0 >> batch.bat`;
			`echo cls >> batch.bat`;
			`echo echo :: External program cd-hit.exe :: >> batch.bat`;
			`echo echo ---------------------------------------------- >> batch.bat`;
			
			`echo echo 	Input		: $root_path/accepted_seq_step-1.fasta >> batch.bat`;
			`echo echo 	Sequence identity threshold	: $cd_hit_identity\% >> batch.bat`;
			`echo echo 	Word_length	: 3 >> batch.bat`;
			`echo echo 	Other parameters	: Default (refer to cd-hit website) >> batch.bat`;
			
			
			`echo echo Please wait.......... >> batch.bat`;
			`echo $cdhit_path -i $root_path\\accepted_seq_step-1.fasta -o $root_path\\cdHit_out -d 0 -c $cdHit_c -n 3 >> batch.bat`;
			`echo rename $root_path\\cdHit_out.clstr cdHit_out.clstr.txt >> batch.bat`;	##mv works
			`echo exit >> batch.bat`;
			system("start batch.bat ");
			
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
		
	$blast_path=~ s{/}{\\}g; $blast_path='"'.$blast_path.'"';
	$Hproteome_file=~ s{/}{\\}g;  $Hproteome_file='"\"'.$Hproteome_file.'\""';

	unlink "$L_root_path/host_orthologs_blast1.out.txt";
	my($gap_open, $gap_extns)=split /,/,$gap_score_1;
	my $blast1="$blast_path -p blastp -d $Hproteome_file -i $root_path\\accepted_seq_step-2.fasta -e $e_val_1 -m $out_fmt_1 -W $word_size_1 -M $sub_matrix_1 -G $gap_open -E $gap_extns -o $root_path\\host_orthologs_blast1.out -a $use_cores -f $threshold_1"." $extra_params_BLAST1";
			`echo echo off > batch.bat`;
			`echo color B0 >> batch.bat`;
			`echo cls >> batch.bat`;
			`echo echo :: External program Blastall.exe :: >> batch.bat`;
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
			system("start batch.bat ");

	while(!(-e "$L_root_path/host_orthologs_blast1.out.txt")){
			##wait till human_orthologs_blast1.out.txt is available
		$blast_prg1=blast_progress("$root_path\\host_orthologs_blast1.out","$L_root_path/host_orthologs_blast1.out",$sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-orthologous})); 
		$blast_prg1= 90 if $blast_prg1>90;
		sleep(3);Tkx::update(); 
	}	
	$blast_prg1=100;Tkx::update();
	
	my ($host_like_proteins,$not_host_like_proteins)=process_host_blast_out("$L_root_path/accepted_seq_step-2.fasta","$L_root_path/host_orthologs_blast1.out.txt", $perc_identity_1 );
	$sequence_summary{-host_orthologs}=scalar @$host_like_proteins;
			
	my $non_host_seq_id=fetch_seq_by_id($all_t_seq,$not_host_like_proteins);
	write_fasta_seq($non_host_seq_id,"$L_root_path/accepted_seq_step-3.fasta");
	my $r=fetch_seq_by_id($all_t_seq,$host_like_proteins);
	write_fasta_seq($r,"$L_root_path/excluded_seq_step-3.fasta");
	$entry_host_homolog_seq->delete(0, "end"); $entry_host_homolog_seq->insert(0, $sequence_summary{-host_orthologs}) ;
	$entry_non_host_proteome->delete(0, "end"); $entry_non_host_proteome->insert(0, $sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-paralogslogs}+$sequence_summary{-host_orthologs})) ;
	
	Tkx::update();

	Tkx::tk___messageBox(-message => "Select an target prediction approach\nPerform either 'Sequence-based approach'  or 'PPI network-based approach' ");
	
###do_PPI_search	
	$do_PPI_search->configure(-command=>sub{
			if(!$PPI_id_map_file){
				Tkx::tk___messageBox(-type => "ok",
					-message => "Probably you forgot to load STRING mapping file.\n",
					-icon => "error", -title => "Input file missing");
				my $PPI_id_map_file1=Tkx::tk___getOpenFile();	
				$PPI_id_map_file=$PPI_id_map_file1;
			}
			if(!$interactome_file){
				Tkx::tk___messageBox(-type => "ok",
					-message => "Probably you forgot to load interaction file.\n",
					-icon => "error", -title => "Input file missing");
				my $interactome_file1=Tkx::tk___getOpenFile();	
				$interactome_file=$interactome_file1;
			}
			$do_PPI_search->configure(-state=>"disabled");
			$do_ess_pro_blast->configure(-state=>"disabled");
			
			$string_srch_prg=0;
			$sqlite_path=~ s{/}{\\}g;  #$sqlite_path='"\"'.$sqlite_path.'\""';;
			#$interactome_file=~ s{/}{\\}g;
			#import PPI file to  sql database 25%
			#calculations; 40%
			#id mapping to $L_root_path/accepted_seq_step-3.fasta; 25%
			die "not_host_like_proteins no datya $not_host_like_proteins\n" if (scalar @$not_host_like_proteins <2);
			
			#sorting and mapping to accepted_seq_step-2.fasta ids; 10% $non_host_seq_id
			Tkx::tk___messageBox(-type => "ok",
					-message => "Calculating centrality measures\nWindows may freez for sometime. Please wait...\n",
					-title => "Import essential proteome");
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
			print S "$sqlite_path $database <import.sql\n";
			#print S "copy import.sql  mm.sql\n";
			print S "del import.sql  \n";
			print S "exit  \n";
			close S;
			system("start batch.bat ");
			while(-e "import.sql"){	$string_srch_prg=5;Tkx::update();sleep(1);	}
			$string_srch_prg=35;Tkx::update();

			my $g = Graph::Undirected->new(); # An undirected graph.
			my $driver   = "SQLite";
			my $dsn = "DBI:$driver:dbname=$database";
			my $userid = "";
			my $password = "";
			my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 })
								  or die $DBI::errstr;
			foreach my $i(@$not_host_like_proteins)
			{
				
				my $p=$seq_id_to_PPI_id_map{$i};
				#print "$i --> $p -\n";
				next if !$p;
				my $stmt = qq(INSERT INTO tmp (id) VALUES ('$p'););	
				my $rv = $dbh->do($stmt) or die $DBI::errstr;
			
			}
			
			my $stmt = qq(select PPI.proteinA, PPI.proteinB, PPI.score from PPI INNER JOIN tmp where (tmp.id=PPI.proteinA) AND PPI.score >=$PPI_score_cutoff);	##LIMIT 1500
					
			my $sth = $dbh->prepare( $stmt );
			my $rv = $sth->execute() or die $DBI::errstr;
			if($rv < 0){
					print $DBI::errstr;
				}
			while(my @row = $sth->fetchrow_array()) {
				$g->add_edge($row[0],$row[1]);			
			}
			
		##Calculating nodes		
			$string_srch_prg=40;Tkx::update();
			my @V = $g->vertices;		## nodes in array
			my $V = scalar @V;			##$g->vertices;		## network sizr (n)
			my $gd = $g->diameter;		#diamemter
			
			open (K,">all_node.txt") or die "$!";
			$"="\n";
			print K "@V\n";
			close K;
			$"=" ";
			sleep(2);
			my $apsp = $g->APSP_Floyd_Warshall();
			Tkx::update();
			#$prg_grph_cc->g_destroy();
			#$mw->configure(-cursor=>"arrow");
			#Return the all-pairs shortest path object computed from the graph using Floyd-Warshall's algorithm. The length of a path between two vertices is the sum of weight attribute of the edges along the shortest path between the two vertices. If no weight attribute name is specified explicitly the attribute weight is assumed.
			
			$string_srch_prg=50;Tkx::update();
			my %centrality_scores_per_node;
			my %total_normalized_score_per_node;		##in a separate one as can be sorted easily
			my($tot_radiality,$tot_Ec,$tot_closeness)=( 0, 0, 0);
			my($min_radiality,$min_Ec,$min_closeness,$max_radiality,$max_Ec,$max_closeness)=( 0, 0, 0,0,0,0);
			my %node_componenet_index_tograph;			## keeps componenet idex as key and created undirected graph object as value;later updated to diameter of the object;just to save time at the cost of memory
			
##Cytohubba
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
						system("$sqlite_path $database < import.sql");
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
				foreach (@Connected_v ) {print SIF "$seq_id\t$PPI_id_to_seq_id_map{$_}\n";}
			}
			close SIF;

			sleep(1);	
			$string_srch_prg=100;Tkx::update();
			$save_result->configure(-state=>"normal");				
			
			Tkx::tk___messageBox(-message => "Run complete.\n \nSave results and parameter options."); 	

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
					-message => "BLAST database (*.phr,*.psq, *.pin) not found for $Eproteome_file. Use formatdb.exe to create it ",
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
			my $blast2="$blast_path -p blastp -d $Eproteome_file -i $root_path\\accepted_seq_step-3.fasta -e $e_val_2 -m $out_fmt_2 -W $word_size_2 -M $sub_matrix_2 -G $gap_open -E $gap_extns -o $root_path\\essential_protein_blast2.out -f $threshold_2 -a $use_cores"." $extra_params_BLAST2";
				`echo echo off > batch.bat`;
				`echo color B0 >> batch.bat`;
				`echo cls >> batch.bat`;				
				`echo echo :: External program Blastall.exe :: >> batch.bat`;
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
				system("start batch.bat ");

				while(!(-e "$L_root_path/essential_protein_blast2.out.txt")){
				#print "$blast_prg2  $sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-orthologous})--\n";
				my $t=$sequence_summary{-total_seq}-($sequence_summary{-very_short_seq}+$sequence_summary{-orthologous}+$sequence_summary{-host_orthologs});
				$blast_prg2=blast_progress("$root_path\\essential_protein_blast2.out","$L_root_path/essential_protein_blast2.out",$t); 
				sleep(3);
				Tkx::update(); 
				
				}	##wait till essential_protein_blast2.out.txt is available
				#$front_page_status
				$blast_prg2=100;Tkx::update();
				
				
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
				Tkx::tk___messageBox(-message => "Run complete.\n \nSave results and parameter options."); 	
				

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



 
sub broad_spect_run
{
	my $frm=shift;
	my $run_button = shift;;	
	
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
	$sequence_summary{-putative_drug_targets}=count_fasta_seq($input_seq); 	## reset drug target counts
	$input_seq=~ s{/}{\\}g; $input_seq='"'.$input_seq.'"'; 					##Convert to windows format
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
	$root_path = Tkx::tk___chooseDirectory(-parent=>$mw);$mw->g_raise();
	
	if(!$root_path){$root_path=$L_root_path;}  							##if cancel pressed ; fail safe
	$L_root_path=$root_path;													##preserve UNIX format
	$root_path=~ s{/}{\\}g; $root_path='"'.$root_path.'"'; 						##Convert to windows format	
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
	my $blast3="$blast_path -p blastp -d $broad_spectrum_pathogen_db_list -i $input_seq -e $e_val_3 -m $out_fmt_3 -W $word_size_3 -M $sub_matrix_3 -G $gap_open -E $gap_extns -o $root_path\\broad_spe_blast3.out -a $use_cores -f $threshold_3"." $extra_params_BLAST3";
		`echo echo off > batch.bat`;
		`echo color 90 >> batch.bat`;
		`echo cls >> batch.bat`;
		`echo echo :: External program Blastall.exe :: >> batch.bat`;
		`echo echo ---------------------------------------------- >> batch.bat`;
		
			`echo echo Parameters >> batch.bat`;
			`echo echo 	Program: blastp >> batch.bat`;
			`echo echo 	Query		: $input_seq >> batch.bat`;
			`echo echo 	Database	: broad_spectrum_pathogen_db_list (see locale_dat dir.) >> batch.bat`;
			`echo echo 	E-value	: $e_val_3 >> batch.bat`;
			`echo echo 	Scoring matrix: $sub_matrix_3 >> batch.bat`;
			`echo echo 	Gap-penalty (Open,Extension): $gap_open,$gap_extns >> batch.bat`;
			`echo echo 	Word size : $word_size_3 >> batch.bat`;
			`echo echo 	Threshold for extending hits : $threshold_3 >> batch.bat`;
			`echo echo 	CPUs : $use_cores >> batch.bat`;
			
		
		`echo echo Please wait.......... >> batch.bat`;
		`echo $blast3 >> batch.bat`;
		`echo rename $root_path\\broad_spe_blast3.out broad_spe_blast3.out.txt >> batch.bat`;	##mv works
		`echo exit >> batch.bat`;
		system("start batch.bat ");
		$blast_prg3=0;
		while(!(-e "$L_root_path/broad_spe_blast3.out.txt")){
		
		$blast_prg3=blast_progress("$root_path\\broad_spe_blast3.out","$L_root_path/broad_spe_blast3.out",$sequence_summary{-putative_drug_targets}*10 ); ##as this blast is not with one hit per one query, so miscalculation happens, so multipled 10,  10 hits per query
		$blast_prg3=96 if $blast_prg3>96;		##fail safe
		sleep(3);Tkx::update(); 
		}										##wait till essential_protein_blast2.out.txt is available
		$blast_prg3=100;  Tkx::update();
		my $ref_broad_spec_counts=process_broad_spe_BLAST_out("$L_root_path/broad_spe_blast3.out.txt",$broad_spe_BLAST_cutoff); 				##\%hash
		my %broad_spec_counts;					##stores counts per query;; should be GLOBAL??
		foreach my $a( keys %$ref_broad_spec_counts ){
			$broad_spec_counts{$a}=scalar @{$ref_broad_spec_counts->{$a}};
		}
		my (@query_id,@cons_counts,@bar_dclrs);
		my $total_hits=0;
		foreach my $t(sort { $broad_spec_counts{$b} <=> $broad_spec_counts{$a} } keys %broad_spec_counts ){
		#print "ids: $t".$broad_spe_org_code_full_name->{$t}."\n";
			push (@query_id,$t);
			push (@cons_counts,$broad_spec_counts{$t});
			push @bar_dclrs,'black';
		}
		
		$sequence_summary{-broad_spectrum}=scalar @query_id;	##update later on applying filter
		
		if(!(scalar @query_id)){ Tkx::tk___messageBox(-message => "None of the queries are conserved in any of the species"); return();}
		
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
			Tkx::tk___messageBox(-message => "$sequence_summary{-broad_spectrum} queries filtered!!!\nNow click 'Save results' to save sequences."); 
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
		
	Tkx::tk___messageBox(-message => "Run complete.\nChoose minimum number of conserved species and Click on Apply button."); 
	}); ##END RUN button
	
}	## END BROAD spect

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
	$sequence_summary{-putative_drug_targets}=count_fasta_seq($input_seq); 	## reset drug target counts; not called if in project call
	$input_seq=~ s{/}{\\}g; $input_seq='"'.$input_seq.'"'; 					##Convert to windows format
	})->g_grid(-column=>4,-row=>0,-padx=>2,-pady=>1,-sticky=>"wn");
	
	$new_frm->new_ttk__button(-text=>"...",-width=>5,-command=>sub{
		$root_path = Tkx::tk___chooseDirectory(-parent=>$mw);$mw->g_raise();		
		if(!$root_path){$root_path=$L_root_path;}  							##if cancel pressed ; fail safe
		$L_root_path=$root_path;													##preserve UNIX format
		$root_path=~s{/}{\\}g; $root_path='"'.$root_path.'"'; 						##Convert to windows format	
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
			$$run_button->configure(-state=>"disabled");
			unlink "$L_root_path/drug_target_blast4.out.txt";
			my($gap_open, $gap_extns)=split /,/,$gap_score_1;
			my $blast4="$blast_path -p blastp -d $drug_blast_db_names -i $input_seq -e $e_val_4 -m $out_fmt_4 -W $word_size_4 -M $sub_matrix_4 -G $gap_open -E $gap_extns -o $root_path\\drug_target_blast4.out -a $use_cores -f $threshold_4"." $extra_params_BLAST4";
				`echo echo off > batch.bat`;
				`echo color 90 >> batch.bat`;
				`echo cls >> batch.bat`;
				`echo echo :: External program Blastall.exe :: >> batch.bat`;
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
				
				system("start batch.bat ");
				
				while(!(-e "$L_root_path/drug_target_blast4.out.txt")){
				$blast_prg4=blast_progress("$root_path\\drug_target_blast4.out","$L_root_path/drug_target_blast4.out",$sequence_summary{-putative_drug_targets} ); sleep(3); Tkx::update(); 
				}	##wait till blast4.out.txt is available
				$blast_prg4=100;  Tkx::update();
				
				
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




sub create_project
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Create New Project");
	$root_path=get_my_document_path();		##OS specific; change for linux
	#$mw->configure(-cursor=>"watch");
	$crt_win->g_wm_attributes (-topmost=>1);
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
		$frm1->g_destroy;
		$root_path.="/$project_name";												##updating root path
		$L_root_path=$root_path;													##preserve UNIX format
		$root_path=~ s{/}{\\}g; $root_path='"'.$root_path.'"'; 						##Convert to windows format
		my $frm2=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
		$frm2->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
		
		$frm2->new_ttk__label(-text=>"SPECIFY INPUT FILES",-foreground=>'blue')->g_grid(-column=>0,-row=>0,-padx=>2,-pady=>5,-sticky=>"nw");
			
		$frm2->new_ttk__label(-text=>"Pathogen proteome sequence(*.fas)[required]:		")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$Tproteome_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$Tproteome_file = Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>1,-padx=>2,-pady=>1);
			
		$frm2->new_ttk__label(-text=>"Pre-formatted host reference proteome database[required]:	")->g_grid(-column=>0, -row=>2,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$Hproteome_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$Hproteome_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>2,-padx=>2,-pady=>1);
		
		$frm2->new_ttk__label(-text=>"Pre-formatted essential protein database:")->g_grid(-column=>0, -row=>3,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$Eproteome_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$Eproteome_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>3,-padx=>2,-pady=>1);
		
		$frm2->new_ttk__label(-text=>"Comma-separated protein-protein interaction file (*.csv):		")->g_grid(-column=>0, -row=>4,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$interactome_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>4,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$interactome_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>4,-padx=>2,-pady=>1);
		
		$frm2->new_ttk__label(-text=>"Interaction ID mapping file(*.tsv): 		")->g_grid(-column=>0, -row=>5,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$PPI_id_map_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>5,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$PPI_id_map_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>5,-padx=>2,-pady=>1);
		
=i		$frm2->new_ttk__label(-text=>"Taxon id		")->g_grid(-column=>0, -row=>5,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$taxon_id,-width=>40,)->g_grid(-column=>1,-row=>5,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"Validate",-width=>10,-command=>sub{
		
			my $url = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id='.$taxon_id;
			my $content = get $url;  die "Couldn't get $url" unless defined $content;
			my @content=split /\n/,$content;
			foreach(@content){if(/Taxonomy\s+browser\s+\((.+)\)/g){$taxon_id=$taxon_id.": $1";}}
			Tkx::update();		
		})->g_grid(-column=>4,-row=>5,-padx=>2,-pady=>1);
=cut		
		$frm2->new_ttk__label(-text=>"Parameter file (optional):			")->g_grid(-column=>0, -row=>6,-padx=>2,-pady=>1,-sticky=>"nw");
		$frm2 ->new_ttk__entry(-textvariable => \$filter_param_settings_file,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>6,-padx=>2,-pady=>1,-columnspan=>2);
		$frm2->new_ttk__button(-text=>"...",-width=>5,-command=>sub{$filter_param_settings_file=Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise();})->g_grid(-column=>4,-row=>6,-padx=>2,-pady=>1);
		
		my $data_input_but=$frm2->new_button(-text=>"Ok",-width=>10,-command=>sub{
			if(!$Tproteome_file ||!$Hproteome_file){
				Tkx::tk___messageBox(-type => "ok",
					-message => "Probably you forgot to load an input files.\n",
					-icon => "error", -title => "Input file Error");
					$crt_win->g_destroy;
					&create_project();
					return();
			}
			
			$crt_win->g_destroy;$mw->configure(-cursor=>"arrow");
			my $batch;
			if(!((-e "$Hproteome_file\.phr")and(-e "$Hproteome_file\.psq")and(-e "$Hproteome_file\.pin"))){
			$batch="\"$formatdb_path\" -i \"$Hproteome_file\" -p T";
				Tkx::tk___messageBox(-type => "ok",
					-message => "BLAST database (*.phr,*.psq, *.pin) not found for $Hproteome_file. Use formatdb.exe to create it ",
					-icon => "error", -title => "Input file Error");
					$front_page_status.="\n=> BLAST database (*.phr,*.psq, *.pin) \nnot found for host proteome ";			
			}
			if($batch){
				$front_page_status.="\n=> Project not created.\n";			
				$run_but->configure(-state=>"disabled"); 
			}
			else{
			$run_but->configure(-state=>"normal"); 
			$front_page_status.="\n=>Project created:$root_path\n=>Target sequence: Tproteome_file\n=>Host proteome: $Hproteome_file\n=>Essential proteome:$Eproteome_file\n=>STRING ID mapping file:$PPI_id_map_file\n=>Select parameters from 'Settings' menu and Click on Run program.\n";
			mkdir "$L_root_path", 0755;									##mkdir; WINDOWs specfic way //
			}
			
			$dwn_str_anal->entryconfigure("Broadspectrum analysis",-state=>"disabled");   ##$mw is accessible so, 
			$dwn_str_anal->entryconfigure("GO analysis",-state=>"disabled");
			$dwn_str_anal->entryconfigure("Sub-cellular localization",-state=>"disabled");
		
			
			Tkx::update(); 
		},);
		$data_input_but->g_grid(-column=>1, -row=>10,-pady=>5,);
		my $close_proj_but=$frm2->new_button(-text=>"Close",-width=>10,-command=>sub{unlink "$root_path//$project_name"; undef $root_path; undef $project_name; Tkx::update(); $crt_win->g_destroy;$mw->configure(-cursor=>"arrow");},);
		$close_proj_but->g_grid(-column=>2, -row=>10,-padx=>2,-columnspan=>2,-sticky=>'w');
	});	
	$create_proj_but->g_grid(-column=>1, -row=>5,-pady=>10,-columnspan=>2);
	my $close_proj_but=$frm1->new_button(-text=>"Close",-width=>18,-command=>sub{$crt_win->g_destroy;$mw->configure(-cursor=>"arrow");},);
	$close_proj_but->g_grid(-column=>2, -row=>5,-padx=>10,-columnspan=>2);
}



sub open_project
{




}


sub settings
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Pramaeter settings");
	#$mw->configure(-cursor=>"watch");
	$crt_win->g_wm_attributes (-topmost=>1);
	$crt_win->g_raise();
	
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
	my $apply_change=$frm1->new_button(-text=>"OK",-width=>8,-command=>sub {$crt_win->g_destroy;$mw->configure(-cursor=>"arrow");})->g_grid(-column=>0,-row=>25,-padx=>18,-pady=>12,-columnspan=>2,);
	
}

sub executable_paths
{
	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Binary file path");
	#$mw->configure(-cursor=>"watch");
	$crt_win->g_wm_attributes (-topmost=>1);
	my $frm1=$crt_win->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	
	$frm1->new_ttk__label(-text=>"BLAST bin path:")->g_grid(-column=>0,-row=>1,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => \$blast_path,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>1,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1->new_ttk__button(-text=>"Choose location",-command=>sub{$blast_path = Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise(); })->g_grid(-column=>3,-row=>1,-padx=>2,-pady=>1);
	
	$frm1->new_ttk__label(-text=>"CDhit bin path:")->g_grid(-column=>0,-row=>2,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => \$cdhit_path,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>2,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1->new_ttk__button(-text=>"Choose location",-command=>sub{$cdhit_path = Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise(); })->g_grid(-column=>3,-row=>2,-padx=>2,-pady=>1);
	$frm1->new_ttk__label(-text=>"SqLite bin path:")->g_grid(-column=>0,-row=>3,-padx=>2,-pady=>1,-sticky=>"nw");
	$frm1 ->new_ttk__entry(-textvariable => \$sqlite_path,-width=>40,-state=>"disabled",)->g_grid(-column=>1,-row=>3,-padx=>2,-pady=>1,-columnspan=>2);
	$frm1->new_ttk__button(-text=>"Choose location",-command=>sub{$sqlite_path = Tkx::tk___getOpenFile(-parent=>$crt_win);$crt_win->g_raise(); })->g_grid(-column=>3,-row=>3,-padx=>2,-pady=>1);
	
	my $apply_change=$frm1->new_button(-text=>"Apply",-width=>5,-command=>sub 
		{
		##Chek if exe files exits
			if(!$blast_path ||!$cdhit_path||!$sqlite_path){
				Tkx::tk___messageBox(-type => "ok",
					-message => "Probably you forgot to specify a path.\n",
					-icon => "error", -title => "Input file missing");
					$crt_win->g_destroy;
					&executable_paths();
					return();
				}
			my $flg=1;
			if(!(-e $blast_path)){Tkx::tk___messageBox(-type => "ok",-message => "$blast_path does not exits",
						-icon => "error", -title => "File doesn't exits");$flg=0; }
			if(!(-e $cdhit_path)){Tkx::tk___messageBox(-type => "ok",-message => "$cdhit_path does not exits",
						-icon => "error", -title => "File doesn't exits");$flg=0;  }			
			if(!(-e $sqlite_path)){Tkx::tk___messageBox(-type => "ok",-message => "$sqlite_path does not exits",
						-icon => "error", -title => "File doesn't exits");$flg=0;  }
			if($flg){$crt_win->g_destroy;$mw->configure(-cursor=>"arrow");}
			else{$crt_win->g_raise;}
			
		})->g_grid(-column=>1,-row=>20,-padx=>2,-pady=>5,-sticky=>"nw");
	my $apply_change=$frm1->new_button(-text=>"Reset",-width=>5,-command=>sub 
		{
		##Chek if exe files exits
			my $flg=1;
			($blast_path,$cdhit_path,$formatdb_path,$sqlite_path)=("$installation_path/executables/blastall.exe","$installation_path/executables/cd-hit.exe","$installation_path/executables/formatdb.exe","$installation_path/executables/sqlite3.exe");
			
		})->g_grid(-column=>2,-row=>20,-padx=>2,-pady=>5,-sticky=>"nw");	
}


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



sub about
{

	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("About");
	$crt_win->g_raise();
	$crt_win->g_wm_attributes (-topmost=>1);
	my $heading = $crt_win->new_ttk__label(-text=>"Exogeneous Drug Target Identification Tool\nVersion: 1.0\nLast update: $last_update\n",-justify=>"center",-foreground=>"blue");
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
		$crt_win->g_destroy();
		})->g_grid(-column=>0,-row=>6,-padx=>1,-pady=>5,-sticky=>"ne");	
	
}

##args:title, parent window ref, query name,db name, ref.e-val, ref.outformat, ref.word size, ref.sub matrix, ref.gap score,ref.Threshold, ref.blast_identy,ref.extra
##returns: none; updates all values parsed to it
#view_update_blast_params ($title, $parent, $q, $db, $e_val, $outfmt, $word_size, $sub_mat, $gap_score, $threshold,$extra)
sub view_update_blast_params 
{
	my ($title, $parent, $q, $db, $e_val, $outfmt, $word_size, $sub_mat, $gap_score, $threshold,$perc_identity,$extra)=@_;	
	my $crt_win_blast2 =$$parent->new_toplevel();
	$crt_win_blast2->g_wm_title("$title");
	#$mw->configure(-cursor=>"watch");
	$crt_win_blast2->g_wm_attributes (-topmost=>1);
	$crt_win_blast2->g_raise();
	
	my $frm1=$crt_win_blast2->new_ttk__frame(-borderwidth=>2,-relief=>'sunken',);
	$frm1->g_grid(-row=>0,-column=>0,-sticky=>"nsew");
	#host, query, e-value, matrix, word size, gap open, gap extend, output format
	
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
		
		$frm1->new_ttk__label(-text=>"Enter extra blastall parameters:")->g_grid(-column=>0,-row=>11,-padx=>2,-pady=>1,-sticky=>"nw");
		$BLAST2_params_ele->g_grid(-column=>0,-row=>12,-padx=>2,-pady=>1,-columnspan=>2);
		
		my $apply_change=$frm1->new_button(-text=>"Close",-width=>5,-command=>sub 
		{
			$crt_win_blast2->g_destroy;$$parent->g_raise();
		})->g_grid(-column=>1,-row=>15,-padx=>2,-pady=>5,-sticky=>"nw");

}


##Args:
##return:
sub down_str_anal
{

	my $crt_win =$mw->new_toplevel();
	$crt_win->g_wm_title("Down-stream analysis Settings");
	#$mw->configure(-cursor=>"watch");
	
	$crt_win->g_raise();
	$crt_win->g_wm_attributes (-topmost=>1);
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
	$frm1->new_ttk__combobox(-textvariable =>\$tax_level ,-values=> "Phylum Class Order Family Genus", -width=>12,-state=>'readonly')->g_grid(-column=>0,-row=>5,-padx=>2,-pady=>1,-sticky=>"nw");
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
			$broad_spectrum_pathogen_db_list=create_broad_spe_db_array($ref_brd_sel_db_array,"$installation_path/local_dat/COMPLETE_PTM_PATHPGENS");
	
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
			$drug_blast_db_names=create_drugTarget_blast_db($ref_drug_db_array,"$installation_path/local_dat/KNOWN_DRUG_TARGETS");
	
		});
	$frm1->new_ttk__label(-text=>"BLAST settings against drug-target database:",)->g_grid(-column=>0,-row=>20,-padx=>2,-pady=>5,-sticky=>"nw");	
	$frm1->new_button(-text=>"Change/view",-width=>10,-command=>sub{
		view_update_blast_params ("BLASTp against drug-target database", \$crt_win, "Input.fasta", "Drug target list", \$e_val_4, \$out_fmt_4, \$word_size_4, \$sub_matrix_4, \$gap_score_4, \$threshold_4,\$perc_identity_4,\$extra_params_BLAST4);	
	
	},)->g_grid(-column=>1,-row=>20,-padx=>2,-pady=>5,-sticky=>"nw",-columnspan=>2);
	
	
	
	
	my $ok=$frm1->new_button(-text=>"Close",-width=>10,-command=>sub 
		{
		$crt_win->g_destroy();
		})->g_grid(-column=>0,-row=>30,-padx=>1,-pady=>5,-sticky=>"ne");	
}


sub fetch_tax_names
{
	my $level=shift;
	my $sql="\"SELECT DISTINCT $level from taxonomy\"";
	my $l=`executables/sqlite3.exe $broad_spectrum_pathogen_db_sq $sql`;	
	my @l = split /\n/,$l;
	$l = join('} {',@l); $l = " {".$l."} ";
	return ("$l",\@l);
	
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
	my $sql="\"SELECT fasta_file from taxonomy WHERE $tax_level IN $range\"";
	my $l=`executables/sqlite3.exe $broad_spectrum_pathogen_db_sq $sql`;	
	my @l = split /\n/,$l;
	
	my @db;
	foreach(@l)
	{
		my $r='\"'.$path.'\\'.$_.'\"';
		push @db,$r;
	}
	
	$db_list=join(" ",@db); #$db_list.="-"; $db_list=~s/\"-//g;
	$db_list='"'.$db_list.'"';
	#print "$db_list\n";
	return $db_list;
}



##args:dat file;  code \t full name; strictly tab separated
##returns: an ref of hash
sub read_broad_spe_codes
{
	my $file=shift;
	my %h;
	open(F,"<$file") or die "$! $file\n  ";
	while(<F>)
	{
		chomp;
		next if(/^#/ or !$_);
		my @l=split /\t/,$_;
		$h{$l[0]}=$l[1];
	}
	close F;
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

##args: cdhit clstroutput
##returns: reference of two arrays;\@paralogs,\@uniq

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
sub process_cdHit_clstr_1			##Incorrect method, including the fist hit of a cluster
{
my $file=shift;
my (@paralogs,@uniq);

open (P, "<$file") or die "$!$file";
my $clstr=0;
my @a;
while(<P>){
		if(/^>Cluster\s+\d+/){  next;}
		elsif(/(\d+)\s*\d+aa,\s*>(\S+)\.\.\.\s+/){if($1==0){push @uniq,$2;}else{push @paralogs,$2} }
		else {warn "err in CDHIT cluster\n"}
		
	}
	#if(scalar @a==1){push @uniq,@a; } elsif(scalar @a>1){ push @paralogs,@a; }	##process last line
	return (\@paralogs,\@uniq);
}

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
 `copy $out1 $temp_file`;	##just to tacle file lock; may increase program time

 
 if(-e $out2){
	my $processed_seq=int `findstr /R /N "^" $temp_file | find /C ":"` ; 
	return ($processed_seq/$total_seq)*100;
	`del $temp_file`;
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


##Args: Win formated path; STRING ("C:\\Users\\SGPGI.SGPGI-PC\\KANHU\\EDIT\\Prev_data";)
##returns: Unix formated path string (C:/Users/SGPGI.SGPGI-PC/KANHU/EDIT/Prev_data)
sub unix_path
{
	my $p = shift;
	#chomp ($p);		##
	$p=~s/\\/\//g;		## \\ to /
	$p=~s/"//g;
	return "$p";			##
}

Tkx::MainLoop();
