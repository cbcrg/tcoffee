#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"
#include "dev1_lib_header.h"  //JM_STRAT

#define ACTION(x) ((n_actions>=(x+1))?action_list[x]:NULL)
#define ACTION2(x,y) ((n_actions>=(x+1))?action_list[x]:y)
#define ATOI_ACTION(x) ((ACTION(x)!=NULL)?(atoi(ACTION(x))):0)


/**
 * \file reformat.c
 * Auxiliary functions to format sequences etc.
 *
 */
/**************************************************************************************************/
/*****************************    SEQ_REFORMAT     ******************************************/
/**************************************************************************************************/
int output_transitions(char *outfile, Alignment *A);
static int output_age_matrix ( char *outfile, int val);
int SeqGCGCheckSum(char *seq, int len);
static Sequence *seq2year ( Sequence *S, int modulo);
static Sequence* output_n_pavie_age_channel (Sequence *S, char *name, int n);
static Sequence* output_pavie_age_channel (Sequence *S, char *name, int modulo);

static int output_seq2struc(char *outfile, Alignment *A);
void output_conservation_statistics ( char *file, Alignment *A);

/*-- Maria added this because it was needed for the compiler to now before the linkage that these functions exist and are declared somewhere--*/
/*--The first 3 are declared at the end of this file--*/
int is_stockholm_aln (char *file);
int fast_format_determination  ( char *in_f);
int output_header_mat (int **mat, char *fname);

int tree2nnode_unresolved (NT_node R, int *l);      			/* Is declared in tree_util.c */
NT_node redundate (Sequence* S,NT_node T, char *seq, char *tree);	/* Is declared in util_make_tree_.c */
int  sp_triplet_coffee_evaluate_output ( Alignment *IN,Constraint_list *CL, char *fname);  /* Is declared in evaluate.c */
/**************************************************************************************************/
/*****************************    SEQ_REFORMAT     ******************************************/
/**************************************************************************************************/
int big()
{
  if (getenv ("BIG_4_SEQ_REFORMAT"))return 1;
  return 0;
}
int seq_reformat ( int argc, char **in_argv)
 	{

        Sequence_data_struc *D1=NULL;
	Sequence_data_struc *D2=NULL;
	Sequence_data_struc *D_ST=NULL;
	Action_data_struc  *RAD;



	int a, b;

 	char *in_format;
	char *in2_format;
	char *out_format;
 	char *in_file;
 	char *in2_file;
 	char *out_file;
 	char *out2_file;
 	char *struc_in_format;
 	char *struc_out_format;
	char *struc_in_file;
 	char *struc_out_file;
	char**action_list;
	char **action;
	char *rename_file;
	char *cache;
	char ***rename_list=NULL;
	int code=CODE;
	char **argv;
	int maxlen;
	
	int n_actions=0;
	int print_format=0;
	/*INITIALIZATIONS*/

	RAD=(Action_data_struc*) vcalloc ( 1, sizeof ( Action_data_struc));
	RAD->keep_case=1;
	declare_name (cache);sprintf ( cache, "use");
	declare_name(in_file);
	declare_name(in2_file);
	declare_name(out_file);
	declare_name(out2_file);
	declare_name(struc_in_format);
	declare_name(struc_out_format);
	declare_name(RAD->coor_file);

	declare_name(struc_in_file);
	declare_name(struc_out_file);
	declare_name(in_format);
	declare_name(in2_format);
	declare_name(out_format);
	declare_name(rename_file);

	maxlen=0;
	for (a=0; a<argc; a++)
	  if (strlen (in_argv[a])>maxlen)maxlen=strlen (in_argv[a]);
	argv=break_list ( in_argv, &argc, "=;, \n");

	action_list=declare_char (argc+1,maxlen+1);

/*END INITIALIZATION*/

 	addrandinit ( (unsigned long) 500);

 	if ( argc==1 || strm6 ( argv[1], "h", "-h", "help", "-help", "-man", "?"))
 		{

		fprintf ( stdout, "\n%s (%s,%s,%s [%s])\n",PROGRAM, VERSION,AUTHOR, DATE, URL);
		fprintf ( stdout, "\n***********     MINIMUM SYNTAX        *****************");
		fprintf ( stdout, "\nseq_reformat -in <in_file> -output <out_format>");
		fprintf ( stdout, "\nSome File formats are automatically recognised");
		fprintf ( stdout, "\nSee Format section");
		fprintf ( stdout, "\n");
		fprintf ( stdout, "\n***********        MAIN FLAGS              ******************");
		fprintf ( stdout, "\n-in     name........Name of the file read");


		fprintf ( stdout, "\n-input  format......Name of the format read, see Input Format Section");
		fprintf ( stdout, "\n...................Automatic detection, except for seqs of numbers");
		fprintf ( stdout, "\n...................-input number_aln | number_fasta");
		fprintf ( stdout, "\n-in2    fname......Second alignment");
		fprintf ( stdout, "\n-input2 format.....See -input");
		fprintf ( stdout, "\n-exon_boundaries   obj file");
 		fprintf ( stdout, "\n-out    fname......Output file (default is STDOUT");
 		fprintf ( stdout, "\n-output format.....Output Format, default is fasta_aln");
		fprintf ( stdout, "\n-struc_in   name...File containing a coded aln");
		fprintf ( stdout, "\n-struc_in_f format.See -input and output format section");
 		fprintf ( stdout, "\n-struc_out  fname..Name of the output structure");
		fprintf ( stdout, "\n-struc_out_f symbol");
		fprintf ( stdout, "\n-keep_case=on|off..keep case, On by default");
		fprintf ( stdout, "\n-action +ac1 +ac2..See the action section");
		fprintf ( stdout, "\n-rename <file>.....Rename the sequences following <file> indications");
		fprintf ( stdout, "\n...................File Format: One couple <name1><space><name2>/line");
		fprintf ( stdout, "\n...................Rename order <name1> into <name2>");
		fprintf ( stdout, "\n...................code file: -output code_name");
		fprintf ( stdout, "\n-code   <file>     Rename file <name1> to <name2>");
		fprintf ( stdout, "\n-decode <file>     Rename file <name2> to <name1>");
		fprintf ( stdout, "\n-no_warning........Suppresses all warnings");
		fprintf ( stdout, "\n-cache.............use,ignore,update,local, DirectoryName");
		fprintf ( stdout, "\n-treemode..........mode used to generate any required tree");
		fprintf ( stdout, "\n-keep..............Number of top entries to be kept when trimming");
		


		fprintf ( stdout, "\n");

		fprintf ( stdout, "\n***********     REFORMAT ACTIONS               *****************");
		fprintf ( stdout, "\n     +Xaction.............Specifies which file undergoes the action");
		fprintf ( stdout, "\n     +Xaction.............X=1: -in");
		fprintf ( stdout, "\n     +Xaction.............X=2: -in2");
		fprintf ( stdout, "\n     +Xaction.............X=3: -struc_in");
		fprintf ( stdout, "\n     +name2unique_name....replace duplicated name with name_#");
		fprintf ( stdout, "\n     +swap_header........,swapp comments: replace comments/name in 1 by in 2");
		fprintf ( stdout, "\n     +swap_lib_header.F...Replace the sequences in the tc_lib (-in) with those in F");
		fprintf ( stdout, "\n     .....................F is a legal FASTA file");


		fprintf ( stdout, "\n     +translate[0-2]......Translate on Frame 0, 1, 2 ");
		fprintf ( stdout, "\n     +translate[3]........longuest ORF on direct strand");
		fprintf ( stdout, "\n     +translate[4]........longuest ORF on direct+complementary strand");


		fprintf ( stdout, "\n     +add_scale..<offset>.addscale below aln");

		fprintf ( stdout, "\n     +rm_gap n ...........Removes col with n%% gap [n=100]");
		fprintf ( stdout, "\n     +rm_lower c..........Removes lower case positions and replaces them with c [-]");
		
		fprintf ( stdout, "\n     +rmgap_col SEQ1:SEQ2.Removes column with a gap in SEQ [#] ");

		fprintf ( stdout, "\n     +backtranslate.......Random Backtranslation");
		fprintf ( stdout, "\n     +complement..........Produces the reverse complement");
		fprintf ( stdout, "\n     +shuffle.<N>.<n>.<m>.Produces N replicates named <n>.<N>.shuffled.fa");
		fprintf ( stdout, "\n     .....<m>: seq or aln.Replicates keep gaps (aln) or not (seq). Default is seq");
		fprintf ( stdout, "\n     +reorder.............Reorders sequences of <in> according to <in2>");
		fprintf ( stdout, "\n     .........random......Random_order");
		fprintf ( stdout, "\n     .........tree........Tree Order (in2)");
		fprintf ( stdout, "\n     +reorder_column.....Reorders sequences of <in> according to <in2>");
		fprintf ( stdout, "\n     .........random......Random_order");
		fprintf ( stdout, "\n     .........tree..mode..Tree Order (comuted with mode: sarmat, idmat, blosum62mt...");
		fprintf ( stdout, "\n     +aln2random_aln SCR..Randomize the aln, S: swap sequences names");
		fprintf ( stdout, "\n     .....................Swap residues within colums");
		fprintf ( stdout, "\n     .....................Swap residues across the aln");
		fprintf ( stdout, "\n     +aln2sample......N......");
		fprintf ( stdout, "\n     +aln2bootstrap...N......");


		fprintf ( stdout, "\n     +chain...............Identifies all the intermediate sequences from <-in>");
		fprintf ( stdout, "\n     .....................needed to join every sequence pair in <-in2>");

		fprintf ( stdout, "\n     +aln2cons  mat_name..Ouputs a consensus sequence");
		fprintf ( stdout, "\n     .....................The consensus is determined using mat");
		fprintf ( stdout, "\n     .....................By Default, mat=blosum62mt, name=Cons");
		fprintf ( stdout, "\n     +aln2resindex........Prints the sequence index of each residue in -in for each -in2 sequence");
		fprintf ( stdout, "\n     +collapse_aln <new name> <seq1> <seq2...> | file name");
                fprintf ( stdout, "\n     .....................Replaces a group of sequences with its consensus");
		fprintf ( stdout, "\n     .....................The replacement sequence is named <new_seq>");
		fprintf ( stdout, "\n     .....................List of sequences can be provided via a file");
		fprintf ( stdout, "\n     .....................File:>new_name seq1 seq2 seq3....");
		fprintf ( stdout, "\n     +original_seqnos.....Keep original seqnos [SWITCH]");
		fprintf ( stdout, "\n     +seqnos..............Print Seqnos [SWITCH]");
		fprintf ( stdout, "\n     +code_dna_aln........Undocumented")  ;
		fprintf ( stdout, "\n     +grep..[NAME|SEQ|COMMENT]..[KEEP|REMOVE]..[string]......");
		fprintf ( stdout, "\n     .....................Keeps or Removes Sequences matching string");
		fprintf ( stdout, "\n     +extract_block <seq> <start> <end> | <seq> <pos> |<filename>");
		fprintf ( stdout, "\n     .....................Extract column pos OR [start to end[");
		fprintf ( stdout, "\n     .....................<filename> Format");
		fprintf ( stdout, "\n     .......................seq start end | seq pos");
		fprintf ( stdout, "\n     .......................# for comments");
		fprintf ( stdout, "\n     .......................! seq offset_value (0 by default)");
		fprintf ( stdout, "\n     .....................Can extract as many positions as needed");
		fprintf ( stdout, "\n     .....................seq=cons: measure positions on the full aln");
		fprintf ( stdout, "\n     +cat_aln.............Concatenates the alignments input via -in and -in2");
		fprintf ( stdout, "\n     +cat_aln.............-if no -in2, -in is expected to be a list of alignments to concatenate");
		fprintf ( stdout, "\n     +orthologous_cat..<mode>: mode=voronoi or nothing");
		fprintf ( stdout, "\n     ......................-in: sequences from different species");
		fprintf ( stdout, "\n     ..................... -in2: list of species in fasta");
		fprintf ( stdout, "\n     ..................... sequence must be named: <species>_<genename>");
		fprintf ( stdout, "\n     ..................... all paralogues will be concatenated");

		fprintf ( stdout, "\n     +aln2replicate N name");
		fprintf ( stdout, "\n     ..................... Generates N replicates in Fasta");
		fprintf ( stdout, "\n     ..................... Voronoi weights can be used");

		fprintf ( stdout, "\n     +msalist2cat_pwaln.min..max");
		fprintf ( stdout, "\n     .....................extract all pw projections and conctaenates those\n");
		fprintf ( stdout, "\n     .....................where id>=min and id<=max\n");
		fprintf ( stdout, "\n     .....................min and max can be omitted (min=0, max=100)\n");

		fprintf ( stdout, "\n     +seq2blast...........gather all possible homologues from +db using blastp\n");
		fprintf ( stdout, "\n     .....................with +num_iterations using +thread\n");
		fprintf ( stdout, "\n     .....................the output are put in +outdir and can be +compressed with gzip\n");
		fprintf ( stdout, "\n     .....................+outdir can be used as a cache by psicoffee\n");
		fprintf ( stdout, "\n     .....................The folowing flgas must be set BEFORE +seq2blast\n");
		fprintf ( stdout, "\n     .....................BLAST output already in +outdir will NOT be recomputed\n");
		
		
		fprintf ( stdout, "\n     +thread..<int>.......Number of procs used (0 for all, Default=1)\n");       
		fprintf ( stdout, "\n     +db......</path/db>..Database used by Blast. Compulsory\n");
		fprintf ( stdout, "\n     +num_iterations.<int>Number of Psi-Blast iteations (Def=1)\n");
		fprintf ( stdout, "\n     +outfmt..<int>.......BLAST output format\n");
		fprintf ( stdout, "\n     +outdir..<path>......BLAST output dir (one per seq, Def=./\n");
		fprintf ( stdout, "\n     +compress............gzip output\n");
		
		fprintf ( stdout, "\n     +seq2prf.............runs seq2blast and extracts profiles into +outdir/<seqname>.prf\n"); 
		fprintf ( stdout, "\n     .....................each profileis a BLAST stack on the query where columns gapped in query are removed\n");
		fprintf ( stdout, "\n     .....................the blast seq kepped are thopse with >+prot_min_cov (40), <+prot_max_sim (50) and >+prot_min_sim(100)\n");
		fprintf ( stdout, "\n     .....................The remaining sequences are trimmed to the +trimseq (100) representatives using +trimseq_mode (regtrim)\n");
		fprintf ( stdout, "\n     .....................while using a tree as defined in regtrim_tree(codnd)\n");
		fprintf ( stdout, "\n     +prot_min_cov.<int>..minimum coverage on query for BLAST hit to be kept Def: 40%\n");
		fprintf ( stdout, "\n     +prot_min_sim.<int>..minimum sim with query for BLAST hit to be kept Def: 50%\n");
		fprintf ( stdout, "\n     +prot_max_sim.<int>..max sim with query for BLAST hit to be kept Def: 90%\n");
		fprintf ( stdout, "\n     +psitrim.<int>.......max number of hits in final psicoffee profile Def: 100\n");
		fprintf ( stdout, "\n     +psitrim_tree.<treemode|tree>..tree used to psitrim Def: codnd\n");
		fprintf ( stdout, "\n     +psitrim_mode.<treemode|tree>..mode used to psitrim Def: regtrim\n");
		
		
		
		
 
		
		
		
		fprintf ( stdout, "\n     .....................with +num_iterations using +thread\n");
		fprintf ( stdout, "\n     .....................the output are put in +outdir and can be +compressed with gzip\n");
		fprintf ( stdout, "\n     .....................+outdir can be used as a cache by psicoffee\n");
		fprintf ( stdout, "\n     .....................The folowing flgas must be set BEFORE +eq2blast\n");
		
		fprintf ( stdout, "\n     +thread..<int>.......Number of procs used (0 for all, Default=1)\n");       
		fprintf ( stdout, "\n     +db......</path/db>..Database used by Blast. Compulsory\n");
		fprintf ( stdout, "\n     +num_iterations.<int>Number of Psi-Blast iteations (Def=1)\n");
		fprintf ( stdout, "\n     +outfmt..<int>.......BLAST output format\n");
		fprintf ( stdout, "\n     +outdir..<path>......BLAST output dir (one per seq, Def=./\n");
		fprintf ( stdout, "\n     +compress............gzip output\n");
		
		
		
		
		
		fprintf ( stdout, "\n     +seq2msa <matrix>....makes a standard progressive alignment using matrix");
		fprintf ( stdout, "\n     +realign_block <c1> <c2> <pg>");
		fprintf ( stdout, "\n     .....................Realign column c1 to c2 (non inc.) with pg)");
		fprintf ( stdout, "\n     .....................pg reads fasta and outputs fasta");
		fprintf ( stdout, "\n     .....................pg -infile=<infile> -outfile=<outfile>");
		fprintf ( stdout, "\n     +extract_seq seq_name (start end seq_name start end...) | filename");
		fprintf ( stdout, "\n     .....................seq_name='*': every seq");
		fprintf ( stdout, "\n     .....................start='*'   : real start");
		fprintf ( stdout, "\n     .....................end='*'     : real end");
		fprintf ( stdout, "\n     .....................filename: fasta format");
		fprintf ( stdout, "\n     +extract_seq_list name1 name2");
		fprintf ( stdout, "\n     .....................Extracts entire sequences");
		fprintf ( stdout, "\n     +remove_seq sn1 sn2..Removes sequences sn1, sn2...");
		fprintf ( stdout, "\n     +remove_seq empty....Removes empty sequences (gap only)");
		fprintf ( stdout, "\n     +remove_seq unique...Remove all multiple occurences except the first");
		fprintf ( stdout, "\n     +thread_profile_on_msa <file>");
		fprintf ( stdout, "\n     .....................Threads a list of profiles on corresponding seq");
		fprintf ( stdout, "\n     .....................File: >seqname _R_ <msa file> [nlines]");

		fprintf ( stdout, "\n     +thread_dna_on_prot_aln");
		fprintf ( stdout, "\n     .....................-in DNA.seq and -in2 AA.aln");
                fprintf ( stdout, "\n     +thread_struc_on_aln");
		fprintf ( stdout, "\n     .....................-in structure and -in2 aln");
		fprintf ( stdout, "\n     +use_cons............Use the consensus for n[SWITCH]");
		fprintf ( stdout, "\n     +upper.n|[n1-n2].....n omitted sets everything to upper case");
		fprintf ( stdout, "\n     .....................To use n: provide a number_aln via:");
		fprintf ( stdout, "\n     .....................-struc_in <number_file> -struc_in_f number_aln");
		fprintf ( stdout, "\n     .....................if use_cons is set n, is read on the cons");
		fprintf ( stdout, "\n     .....................n: will upper every residue with a value of n in struc_in");
		fprintf ( stdout, "\n     .....................[n1-n2]: upper residues between n1 and n2");
		fprintf ( stdout, "\n     +lower  n|[n1-n2]....See +upper");
		fprintf ( stdout, "\n     +switchcase  n|[n1-n2]See +upper");
		fprintf ( stdout, "\n     +color_residue <seq> <pos> <color> | file");
		fprintf ( stdout, "\n     .....................File: seq_name pos color");
		fprintf ( stdout, "\n     .....................color: 0-9");
		fprintf ( stdout, "\n     +edit_residue <seq> <pos> <edit> | file");
		fprintf ( stdout, "\n     .....................File: seq_name pos color");
		fprintf ( stdout, "\n     .....................edit: upper|lower|symbol");



		fprintf ( stdout, "\n     +keep   n|[n1-n2]....Only keep residues that have a score between n1 and n2");

		fprintf ( stdout, "\n     +invert..............Inverts the sequences: CAT => TAC");
		fprintf ( stdout, "\n     +rotate name         Rotate an MSA, names each sequence name_col#");
		fprintf ( stdout, "\n     +convert n|[n1-n2] s1 s2 ....");
		fprintf ( stdout, "\n     +merge_annotation.... ");

		fprintf ( stdout, "\n     .....................Converts residues with your alignment");
		fprintf ( stdout, "\n     .....................similar to upper");
		fprintf ( stdout, "\n     .....................s1: ABCDe turns every ABCD into e");
		fprintf ( stdout, "\n     .....................s1: #e turns any residue into e");
		fprintf ( stdout, "\n     aln2short_aln L C S..Turns sequences into shorter sequences");
		fprintf ( stdout, "\n     .....................L: list of residues to keep");
		fprintf ( stdout, "\n     .....................S: Size of Streches replaced by symbol C");


		fprintf ( stdout, "\n     +random n l..........Generates N random sequences of len l");
		fprintf ( stdout, "\n     .....................You must provide a file with -in");
		fprintf ( stdout, "\n     +count n|[n1-n2] s1 s2....");
		fprintf ( stdout, "\n     .....................Counts residues with your alignment");
		fprintf ( stdout, "\n     .....................similar to convert");
		fprintf ( stdout, "\n     +print_format........prints the format name");
		fprintf ( stdout, "\n     +keep_name...........Keep the original sequence name on extraction");

		fprintf ( stdout, "\n     +remove_aa pos Ml Ncycle Random_len");
		fprintf ( stdout, "\n     .....................Randomly modifies an alignment");
		fprintf ( stdout, "\n     .....................pos=0: chosen randomly");
		fprintf ( stdout, "\n     .....................MaxLen of the deletions, Ncycle: number of cycles");
		fprintf ( stdout, "\n     .....................Random_len: 0 sets the len to maxlen, 1 to a random value");
		fprintf ( stdout, "\n     +tree..gap <F|def=0.5> replicates<D|column|def=1> mode <nj|upgma|def=nj>");
		
		fprintf ( stdout, "\n     +remove_nuc.x........Remove Position 1, 2 or 3 of every codon");
		fprintf ( stdout, "\n     +evaluate3D..........strike|distances|contacts");
		fprintf ( stdout, "\n     .....................Uses the -in2 contact_lib or the +pdb2contacts +seq2contacts ");
		fprintf ( stdout, "\n     .....................If none, uses +seq2contacts ");
		fprintf ( stdout, "\n     .....................-output score_ascii, score, score_html,fasta_tree ");
		fprintf ( stdout, "\n     .....................+tree [1] trigers tree computation for -output fasta_tree");
		fprintf ( stdout, "\n     .....................distances..<Max|Def=15 >");
		fprintf ( stdout, "\n     .....................contacts...<Max|Def=1.2> <nb|Def=3>");
		fprintf ( stdout, "\n     .....................contacts|distances provided via -in2=contact_lib");
		fprintf ( stdout, "\n     .....................contacts|distances estimated via +pdb2contacts");
		fprintf ( stdout, "\n     +evaluateTree.<group>.Group is a fasta sequence file");
		
		fprintf ( stdout, "\n     +evaluate matrix..gop..gep");
		fprintf ( stdout, "\n     .....................Make a similarity evaluation with matrix");
		fprintf ( stdout, "\n     .....................use -output=score_ascii, or score_html.");
		fprintf ( stdout, "\n     .....................You can filter on the values");
		fprintf ( stdout, "\n     +evaluate matrix..gop..gep");
		fprintf ( stdout, "\n     .....................Make an SP evaluation with matrix");
		fprintf ( stdout, "\n     .....................Uses Natural Gap penalties");
		fprintf ( stdout, "\n     .....................gop and gep must be negative");
		fprintf ( stdout, "\n     .....................use -output=color_ascii, color_html to get a color display");

		fprintf ( stdout, "\n.....+evaluate_lat........Make a lateral evaluation with matrix");
		fprintf ( stdout, "\n     +msa_weight proc.....Computes weights using the procedure");
		fprintf ( stdout, "\nRNA analysis Post Processing___________________________________________________");
		fprintf ( stdout, "\n     +aln2alifold.........Turns the MSA into a consensus structure");
		fprintf ( stdout, "\n     +add_alifold.........adds an alifold consensus structure");

		fprintf ( stdout, "\n     +alifold2analyze.mode..mode=stat_cache_list_aln_color_html_ps_usegap");
		fprintf ( stdout, "\n     .......................stat: compile Number of compensated mutations");
		fprintf ( stdout, "\n     .......................cache: ascii-code compensated mutations on aln");
		fprintf ( stdout, "\n     .......................html: color-code compensated mutations on aln");
		fprintf ( stdout, "\n     .......................aln: mark compensated mutations on stockholm aln");
		fprintf ( stdout, "\n     .......................usegap: do not ignore positions with gaps");

		fprintf ( stdout, "\n     +RNAfold_cmp.........compares the sec struc of in1 and in2 (computes them with alifold if missing)");

		fprintf ( stdout, "\nMSA Post Processing___________________________________________________");
		fprintf ( stdout, "\n     +force_aln filename|seq1 res1 seq2 res2");
		fprintf ( stdout, "\n     .....................Forces residue 1 of seq1 to be aligned with res2 of seq 2");
		fprintf ( stdout, "\n     .....................In a file, there must be one pair of interaction/line");
		fprintf ( stdout, "\n     +sim_filter[_aln_Ix_iy_Cz_cw <seq>");
		fprintf ( stdout, "\n     ....................._<unaln or aln>, aln is assumed");
		fprintf ( stdout, "\n     ....................._I max identity to seq");
		fprintf ( stdout, "\n     ....................._i min identity to seq");
		fprintf ( stdout, "\n     ....................._C max cov on seq");
		fprintf ( stdout, "\n     ....................._c min cov on seq");
		fprintf ( stdout, "\n     +phylotrim [N|N%%|split] [nj|phyml|trmsd] [trmsd template file]");
		fprintf ( stdout, "\n     +regtrim N[%]........Keep N (or N%) sequences from <-in1> seq using <-in2> tree.");
		fprintf ( stdout, "\n     .....................The first -keep sequences are kept\n");
		fprintf ( stdout, "\n     .....................if no <-in2> tree is provided <-treemode> is used to generate it (def=codnd)");
		
		fprintf ( stdout, "\n     ....................._<seq or aln>, aln is assumed");
		fprintf ( stdout, "\n     ....................._%%%%<max/min_percent_similarity>");
		fprintf ( stdout, "\n     ....................._max Or _min <keep sequences for which sim is the max or the min [Def: _max>");
		fprintf ( stdout, "\n     ....................._cov Or _sim Filter according to the coverage [Def: _sim]");
		fprintf ( stdout, "\n     ....................._n<max_number_of_sequence>       ");
		fprintf ( stdout, "\n     ....................._N<percent_of_sequences_to_keep>");
		fprintf ( stdout, "\n     ....................._T Reorder the sequences according to a tree BEFORE triming");
		fprintf ( stdout, "\n     ....................._Fn Keep only sequences that have AT LEAST ONE residue aligned");
		fprintf ( stdout, "\n     ......................in the n first and n last columns. ");
		fprintf ( stdout, "\n     ....................._O<min sim> Remove outlayers that have less than min average sim with other sequences");
		fprintf ( stdout, "\n     ....................._Kn Forces the n top sequences to be kept");
		fprintf ( stdout, "\n     ....................._P_ Print a summary in stderr");


		fprintf ( stdout, "\n     .....................Keeping Sequences: Sequences provided via -in2 will be kept");

		fprintf ( stdout, "\n     .....................Keeping Sequences: Sequences whose name contains <string> in field fS will be kept");
		fprintf ( stdout, "\n     ....................._f<NAME|SEQ|COMMENT> designates a field");
		fprintf ( stdout, "\n     .....................<string> is a Perl regular expression");
		fprintf ( stdout, "\n     +aln2unalign Mode Penalty Threshold");
		fprintf ( stdout, "\n     .....................Identifies all the streches less conserved than than the average");
		fprintf ( stdout, "\n     .....................Mode: lower|number|unalign Act on all the resiues withs score<Thres");
		fprintf ( stdout, "\n     .....................Penalty: FSA penalty align2unalign, Def=90");
		fprintf ( stdout, "\n     .....................Threshold: Fraction of unaligned residues(0-9) Def=2");

		fprintf ( stdout, "\n     +clean_cdna..........Undocumented");
		fprintf ( stdout, "\n     +clean_maln..........Undocumented");
		fprintf ( stdout, "\nTree Analysis___________________________________________________");


		fprintf ( stdout, "\n     +mafftnewick2newick..replaces names by the index value in -in2=<seq> (1..N)");
		fprintf ( stdout, "\n     +newick2mafftnewick..replaces index values by the names in -in2=<seq> (1..N)");
					  
		fprintf ( stdout, "\n     +newick_suffle.....<N>..Randomly swp left and right node when righting N times the -in Tree");
		fprintf ( stdout, "\n     +newick_randomize..<N>..Randomize leafs of -in Tree...produces N replicates");
		fprintf ( stdout, "\n     +tree................Passes information to +evaluate");
		fprintf ( stdout, "\n     .....................Default: +evaluate returns a tree");
		fprintf ( stdout, "\n     .....<N>.............+evaluate makes N replicates");
		fprintf ( stdout, "\n     .....replicates <N>..+evaluate makes N replicates");
		fprintf ( stdout, "\n     .....mode <nj|upgma>.+evaluate makes N replicates");
		fprintf ( stdout, "\n     .....gap <N float>...+evaluate ignores colum with N float gap");
		fprintf ( stdout, "\n     +columns4tree <file>  Provoide a file where each line is a pair of columns [c1 c2]. The list can be bootstrapped\n");         
		fprintf ( stdout, "\n     +tree2bs.............Add BS support to the original tree drawn from replicates");
		fprintf ( stdout, "\n     +tree2bs.............The replicates must be provided via -in, one newick tree/line");
		fprintf ( stdout, "\n     +print_replicates....Print Replicate trees AFTER the main tree");
		
		
		fprintf ( stdout, "\n     +tree_prune..........Prune the tree -in using the sequences provided via -in2");
		fprintf ( stdout, "\n     +tree_cmp............Compares the tree -in and the tree -in2");
		fprintf ( stdout, "\n     +tree_cmp_list......Compares the tree -in and the tree_list -in2");
		fprintf ( stdout, "\n     .....................Sets the support as boostrap value in the -in tree");

		fprintf ( stdout, "\n     .....................-in and -in2 can contain different taxons");
		fprintf ( stdout, "\n     +tree_scan.P1..P2.....scans alignment <-in> with tree <-in2>)");
		fprintf ( stdout, "\n     ......................+tree_scan help to get P1 information");
		fprintf ( stdout, "\n     ......................+aln2tree help to get P2 information");

		fprintf ( stdout, "\n     .....................-in and -in2 can contain different taxons");
		fprintf ( stdout, "\n     +tree2node.......... Reports the node list along with the split");
		fprintf ( stdout, "\n     ..................... splits can be described with the seq order ");
		fprintf ( stdout, "\n     ..................... provided via -in3=<sequence> ");

		fprintf ( stdout, "\n     +treelist2groups.N....count all topologies within a list of trees");
		fprintf ( stdout, "\n     .....................-in is in fasta format with each name being a newick file");
		fprintf ( stdout, "\n     .....................-in2 can be a list of sequences used to trim the trees");
		fprintf ( stdout, "\n     ......................N can be used to unresolve the trees with Depth N");
		fprintf ( stdout, "\n     +treelist2lti.N.C.....Reports the average stability of each sequence neighborhood");
		fprintf ( stdout, "\n     ......................Species can be selected via -in2 [Fasta file with Taxon names]");
		fprintf ( stdout, "\n     ......................OR the sequences observed in C%% of the files are kept [Def: C=100]");


		fprintf ( stdout, "\n     +treelist2seq.C.......Reports the species observed in C%% of the trees");
		fprintf ( stdout, "\n     +treelist2splits......List and counts all the splits in a list of trees");
		fprintf ( stdout, "\n     ......................splits can be restricted to a list of sequences provided via -in2");
		fprintf ( stdout, "\n     +treelist2dmat.......outputs a distance matrix for a list of trees");

		fprintf ( stdout, "\n     +tree_compute n s....Computes a tree using the MSA provided with -in");
		fprintf ( stdout, "\n     ....................n:0-9, controls the way the MSA is filtered");
		fprintf ( stdout, "\n     ....................s:pam250mt|blosum62mt|categories|entropy");
		fprintf ( stdout, "\n     ....................s:controls the column evaluation in MSA");
		fprintf ( stdout, "\n     +change_distances.f.f:float, sets all the distances to f in the tree");
		fprintf ( stdout, "\n     +change_bootstrap n..:n=0 removes all the bootstrap values");
		fprintf ( stdout, "\n     .....................:n!=0 adds a the value n to every node");
		fprintf ( stdout, "\n     +tree2dpatree........Replaces tree distances with the minimum %%ID in");
		fprintf ( stdout, "\n     .....................the depending subgroup. The ID is measured on an");
		fprintf ( stdout, "\n     .....................-in=TREE -in2=ALN");
		fprintf ( stdout, "\n     +unroot..............Removes the root in the input tree");
		fprintf ( stdout, "\n     +tree2group.N.I.P....Reports all the tree subgroup with at most Nseq");
		fprintf ( stdout, "\n     .....................and at min I%% identity. Output format can be read by");
		fprintf ( stdout, "\n     .....................collapse_tree. New groups are named P_1, P_2...");
		fprintf ( stdout, "\n     +collapse_tree.F.....Collapses trees. F is either a file or a list");
		fprintf ( stdout, "\n     .....................<new name> <seq1> <seq2>...");
		fprintf ( stdout, "\n     +aln2tree............Computes a tree");
		fprintf ( stdout, "\n     ..ktupN|aln|sarmat   ktupN: match size N to estimate distances");
		fprintf ( stdout, "\n     .....................aln: Measures distances on aln");
		fprintf ( stdout, "\n     .....................sarmat: expects in to be a SAR matrix of O and I");
		fprintf ( stdout, "\n     ..nj | cw............Runs Neighbor Joining OR Cw to compute Tree");
		fprintf ( stdout, "\n     ..dpa................Turns the tree into a daptree (+tree2dpatree)");
		fprintf ( stdout, "\n     +node_sort..<name>...Sort leafs of tree n1, by node distance");
		fprintf ( stdout, "\n     +treelist2bs.........Turns a list of tree into a tree with BS support");
		fprintf ( stdout, "\n     .............cons....Runs consense majoprity rule. Distances are the BS values");
		fprintf ( stdout, "\n     .............best....identifies the tree with the best support");
		fprintf ( stdout, "\n     .............first...Reports support for the first tree in the list");  			
		fprintf ( stdout, "\n     +treelist2bs_compare.Measures node support of all T1 trees in T1 and compares with bs support");
		fprintf ( stdout, "\nMatrix Analysis___________________________________________________");
		fprintf ( stdout, "\n     +aln2proba_mat.......Computes the proba of all mutations on MSA provided in fasta format (i.e. each name is a valid MSA file");
		
		fprintf ( stdout, "\n     +aln2mat_diaa........computes a dinucleotide matrix on a list of aln");
		fprintf ( stdout, "\n     +aln2mat.............computes a log odd matrix. Input is a list of MSA provided in fasta format (i.e. each name is a valid MSA file");
		fprintf ( stdout, "\n     +seq2lat_mat.........computes a transition matrix on seq provided via -in");

		fprintf ( stdout, "\nStructure Analysis___________________________________________________");
		
		fprintf ( stdout, "\n     +seq2contacts.mode....Reports contacts or distances in library format");
		fprintf ( stdout, "\n     ......................mode: RNAplfold");
		fprintf ( stdout, "\n     ......................Sequences with a structure are ignored");
		fprintf ( stdout, "\n     ......................OUTPUT: -output=contact_lib");
		fprintf ( stdout, "\n     +struc2nb...D.........Display a list of all the residues D appart");
		fprintf ( stdout, "\n     +pdb2contacts.mode D..Reports contacts or distances in library format");
		fprintf ( stdout, "\n     ......................mode: find_pairs|find_pairs-d|RNAplfold|distances|all|best|count|closest");
		fprintf ( stdout, "\n     ......................If mode is set, sequences without a PDB ignored");
		fprintf ( stdout, "\n     ......................mode [find_pair-p for RNA/PDB-RNAplfold for RNA]");
		fprintf ( stdout, "\n     ......................D: [1.2 | 20 for mode=distances] for proteins");
		
		fprintf ( stdout, "\n     ......................D: probe size for contacts 1.2A by def");
		fprintf ( stdout, "\n     ......................D: MaxDistance for distance 30 A by default");
		fprintf ( stdout, "\n     ......................INPUT: provide template via -in2");
		fprintf ( stdout, "\n     ......................OUTPUT: -output=contact_lib");
		fprintf ( stdout, "\n     +struc2nb...D.........Display a list of all the residues D appart");
		fprintf ( stdout, "\n     +rm_template...V......Removes _[S|G|R]_[template] to sequence names");
		fprintf ( stdout, "\n     ......................V: omitted | sequences <=> Output sequences");
		fprintf ( stdout, "\n     ......................V: template <=> Output templates");

		fprintf ( stdout, "\n     +add_template.F.......Add _[S|G|R]_[template] to sequence names");
		fprintf ( stdout, "\n     ......................F can either be a fasta file or an executable");
		fprintf ( stdout, "\n     ......................F: File: >name _S_ template");
		fprintf ( stdout, "\n     ......................F: executable: pg -infile=<seq> -outfile=<tagged>");
		fprintf ( stdout, "\nMatrix Comparison___________________________________________________");
		fprintf ( stdout, "\n    +mat2cmp...............Returns the correlation coefficient between two matrices");
		fprintf ( stdout, "\n    .......................-in mat1 -input matrix, -in2 mat2 -input2 matrix");
		fprintf ( stdout, "\n***********  INPUT FORMATS: Alignments *****************");
		fprintf ( stdout, "\n     AUTOMATIC RECOGNITION");
		fprintf ( stdout, "\n     perl_xxx:............. runs xxx onto the input file");
		fprintf ( stdout, "\n     xxxx <file> > outfile..xxx reads any formats, outputs fasta");
		fprintf ( stdout, "\n     amps_aln       saga_aln      ");
		fprintf ( stdout, "\n     clustal_aln    fasta_aln     msf_aln  ");
		fprintf ( stdout, "\n     dali_aln       gotoh_aln     pima_aln");
		fprintf ( stdout, "\n     dialign_aln    matrix        conc_aln");
		fprintf ( stdout, "\n     NON AUTOMATIC RECOGNITION (use the -input file to specify the format");
		fprintf ( stdout, "\n     number_aln     newick        mafftnewick");
		
		
		fprintf ( stdout, "\n");
		fprintf ( stdout, "\n***********  INPUT FORMATS: Sequences *****************");
 		fprintf ( stdout, "\n     fasta_seq      dali_seq       pir_seq");
 		fprintf ( stdout, "\n     barton_list_tc amps_sd_scores EST_fasta");
 		fprintf ( stdout, "\n     gor_seq        gor_struc      number_fasta[*]");
		fprintf ( stdout, "\n     tc_lib         pdb_struc");
		fprintf ( stdout, "\n");
		fprintf ( stdout, "\n***********  INPUT FORMATS: Structures   *****************");
 		fprintf ( stdout, "\n    rna_number");
 		fprintf ( stdout, "\n    alifold");
		fprintf ( stdout, "\n***********  OUTPUT FORMATS: Alignments ******************");
		fprintf ( stdout, "\n     compressed_aln saga_aln          clustal_aln");
		fprintf ( stdout, "\n     phylip_aln     relaxed_phylip_aln msf_aln    ");
		fprintf ( stdout, "\n     pir_aln        stockhom_aln      stockholm");
		fprintf ( stdout, "\n     fasta_aln      mfasta_aln");
		fprintf ( stdout, "\n     score....................Tabulated MSA and sequence Score\n");
		fprintf ( stdout, "\n     color_html,color_ps......colored using the struc_in file  ");
		fprintf ( stdout, "\n     color_protogene..........colors codons");
		fprintf ( stdout, "\n     color_exoset.............mixes conservation (gray) and introns (RGB)");
		fprintf ( stdout, "\n     color_pdf      pw_lib_saga_aln tdna_aln");
		fprintf ( stdout, "\n     thread_dna_on_prot_aln");
		fprintf ( stdout, "\n");
		fprintf ( stdout, "\n*********** OUTPUT FORMATS: RNA sequence  ******************");
 		fprintf ( stdout, "\n     vienna2tc_lib  vienna2template  stockholm_aln");
		fprintf ( stdout, "\n*********** OUTPUT FORMATS: sequence  ******************");
 		fprintf ( stdout, "\n     fasta_seq      fasta_seq1     gotoh_seq");
		fprintf ( stdout, "\n     gor_seq        cache_id       fasta_seq2");
		fprintf ( stdout, "\n     tblastx_db1    tblastx_db2    tblastx_db3");

		fprintf ( stdout, "\n*********** OUTPUT FORMATS: weights ******************");
		fprintf ( stdout, "\n     constraints    saga_pw_sd_weights  nseq\n");
 		fprintf ( stdout, "\n");
		fprintf ( stdout, "\n*********** OUTPUT FORMATS: trees   ******************");
		fprintf ( stdout, "\n     newick          dm             newick_dm");
		fprintf ( stdout, "\n     mafftnewick     mafftdndmat");
		
		fprintf ( stdout, "\n     use +print_replicates flag to print the replicates (first line = original)");
		fprintf ( stdout, "\n     with newick_dm, grep \";\" to collect the trees");
 		fprintf ( stdout, "\n");
		fprintf ( stdout, "\n*********** OUTPUT Formats: special  ****************");
		fprintf ( stdout, "\n     len             name               statistics<_hnrglNL>");
		fprintf ( stdout, "\n     sim.............outputs a similarity matrix based on an id comparison of -in");
		fprintf ( stdout, "\n     sim_sarmat......in is sar matrix");

		fprintf ( stdout, "\n     sim_mat_<matname>..Evaluates similarity using the matrix: Blosum62mt, etc");
		fprintf ( stdout, "\n     sim_mat_list.......Produces a list of available matrices\n");
		fprintf ( stdout, "\n     sim_idscore........Makes dp alignment of the sequences using Blosum62mt(prot) or idmat(DNA)");
		fprintf ( stdout, "\n     sim_idscoreDNA.....Makes dp alignment of the sequences using idmat");
		fprintf ( stdout, "\n     sim_sim1(default)..Fraction of identical columns over all columns\n");
		fprintf ( stdout, "\n     sim_sim2...........Fraction of identical columns over shortest seq\n");
		fprintf ( stdout, "\n     sim_sim3...........Fraction of identical columns over longuest seq\n");
		fprintf ( stdout, "\n     sim_gap1...........Fraction of ungapped columns\n");
 
		
			  
		

		fprintf ( stdout, "\n     code_name......Outputs a compact list of names for code/decode");

 		fprintf ( stdout, "\n");


		fprintf ( stdout, "\n");
 	        return EXIT_SUCCESS;
 		}

	argv=standard_initialisation (argv, &argc);


	for ( a=1; a< argc; a++)
 		{
		  if (a==1 && argv[1][0]!='-')
		    {
		      sprintf( in_file, "%s", argv[a]);
		    }
		  else if ( strcmp ( argv[a], "-in_f")==0 ||strm(argv[a],"-input") )
 			{
			if ( strcmp ( argv[a], "-in_f")==0) fprintf ( stdout,"\nWARNING: %s deprecated, use -input instead", argv[a]);

 			sprintf ( in_format, "%s", argv[a+1]);
			a++;
 			}

		else if ( strcmp ( argv[a], "-cache")==0 )
 			{
 			sprintf (cache, "%s", argv[a+1]);

			a++;
 			}


		else if ( strcmp ( argv[a], "-exon_boundaries")==0 )
 			{

			  set_string_variable ("exon_boundaries", argv[a+1]);
			  a++;
			}
		else if ( strcmp ( argv[a], "-overaln_threshold")==0 )
 			{

			  set_int_variable ("overaln_threshold", atoi(argv[a+1]));
			  a++;
			}
		else if ( strcmp ( argv[a], "-overaln_target")==0 )
 			{

			  set_int_variable ("overaln_target", atoi(argv[a+1]));
			  a++;
			}
		else if ( strcmp ( argv[a], "-overaln_P1")==0 )
 			{

			  set_int_variable ("overaln_P1", atoi(argv[a+1]));
			  a++;
			}
		else if ( strcmp ( argv[a], "-overaln_P2")==0 )
 			{

			  set_int_variable ("overaln_P2", atoi(argv[a+1]));
			  a++;
			}
		else if ( strcmp ( argv[a], "-overaln_P3")==0 )
 			{

			  set_int_variable ("overaln_P3", atoi(argv[a+1]));
			  a++;
			}
		else if ( strcmp ( argv[a], "-overaln_P4")==0 )
 			{

			  set_int_variable ("overaln_P4", atoi(argv[a+1]));
			  a++;
			}

		else if ( strcmp ( argv[a], "-in2_f")==0||strm(argv[a],"-input2") )
 			{
			  if ( strcmp ( argv[a], "-in_f")==0) fprintf ( stdout,"\nWARNING: %s deprecated, use -input2 instead", argv[a]);

 			sprintf ( in2_format, "%s", argv[a+1]);
			a++;
 			}
		else if ( strcmp ( argv[a], "-seqnos")==0)
 			{
 			sprintf (action_list[n_actions++], "seqnos");
 			}
		else if ( strcmp( argv[a], "-treemode")==0)
		  {
		    set_string_variable ("treemode", argv[++a]);
		  }
		else if ( strcmp( argv[a], "-keep")==0)
		  {
		    set_int_variable ("keep", atoi(argv[++a]));
		  } 
		else if ( strcmp( argv[a], "-action")==0)
		  {
		    while ((a+1)<argc && argv[a+1][0]!='-')
		      {
			sprintf (action_list[n_actions++], "%s", argv[a+1]);
			a++;
		      }
		  }
		else if ( strcmp ( argv[a], "-keep_case")==0)
		  {
		    if(!NEXT_ARG_IS_FLAG)RAD->keep_case=1;
		    else RAD->keep_case=(strm3(argv[a], "on","ON","On"))?1:0;
		    
		  }
		  
		else if ( strcmp ( argv[a], "-conv")==0)
	                {
			if ( strncmp ( argv[a+1],"set",3)==0)RAD->symbol_list=make_symbols (argv[++a],&(RAD->n_symbol));
			else
			    {
			    RAD->symbol_list=declare_char (STRING, STRING);
			    while(!NEXT_ARG_IS_FLAG)
			          {
				  sprintf ( RAD->symbol_list[RAD->n_symbol], "%s", argv[++a]);
				  RAD->n_symbol++;
				  }
			    }
			}
 		else if ( strcmp ( argv[a], "-struc_in_f")==0 ||strcmp ( argv[a], "-input3")==0 )
 			{
 			sprintf ( struc_in_format, "%s", argv[a+1]);
			a++;
 			}
 		else if ( strcmp ( argv[a], "-out_f")==0 ||strm(argv[a],"-output") )
 			{
			if ( strcmp ( argv[a], "-out_f")==0) fprintf (stdout, "\nWARNING: %s deprecated, use -output instead", argv[a]);
 			sprintf ( out_format, "%s", argv[a+1]);
			a++;
 			}
 		else if ( strm ( argv[a], "-struc_out_f") || strm ( argv[a], "-output_struc") )
 			{
			sprintf ( struc_out_format, "%s", argv[a+1]);
			a++;
 			}
 		else if ( strcmp (argv[a],"-in")==0)
 			{
 			sprintf( in_file, "%s", argv[a+1]);
 			a++;
 			}
		else if ( strcmp (argv[a],"-rename")==0)
 			{
 			sprintf( rename_file, "%s", argv[a+1]);
 			a++;
 			}
		else if ( strcmp (argv[a],"-code")==0)
 			{
			code=CODE;
			sprintf( rename_file, "%s", argv[a+1]);
			a++;
 			}
		else if ( strcmp (argv[a],"-decode")==0)
 			{
			  code=DECODE;
			  sprintf( rename_file, "%s", argv[a+1]);
			  a++;
			}
 		else if ( strcmp (argv[a],"-in2")==0)
 			{
 			sprintf( in2_file, "%s", argv[a+1]);
 			a++;
 			}
		else if ( strcmp (argv[a],"-coor")==0)
 			{
 			sprintf( RAD->coor_file, "%s", argv[a+1]);
			a++;
 			}
 		else if (strcmp (argv[a],"-out")==0)
 			{
 			sprintf (out_file, "%s", argv[a+1]);
 			a++;
 			}
 		else if (strcmp (argv[a],"-out2")==0)
 			{
 			sprintf (out2_file, "%s", argv[a+1]);
 			a++;
 			}
 		else if ( strcmp (argv[a],"-struc_in")==0 || strcmp (argv[a],"-in3")==0 )
 			{
 			sprintf( struc_in_file, "%s", argv[a+1]);
 			a++;
 			}
 		else if (strcmp (argv[a],"-struc_out")==0)
 			{
 			sprintf (struc_out_file, "%s", argv[a+1]);
 			a++;
 			}
 		else if ( strcmp ( argv[a], "-rm_gap")==0)
 			{
 			RAD->rm_gap=1;
 			}
		else if ( strcmp ( argv[a], "-print_format")==0)
 			{
 			print_format=1;
 			}
		else if ( strcmp ( argv[a], "-no_warning")==0)
		        {
			set_warning_mode (NO);
 			}
		else if ( strcmp (argv[a], "-setenv")==0)
		        {
			  while ((a+2)<argc && argv[a+1][0]!='-' && argv[a+2][0]!='-')
			    {
			      cputenv ("%s=%s", argv[a+1], argv[a+2]);
			      a+=2;
			    }
			}
		else if ( strcmp (argv[a], "-big")==0)
		  {
		    cputenv ( "BIG_4_SEQ_REFORMAT=1");
		  }
 		else
 			{
 			fprintf ( stdout, "\nUNKNOWN OPTION: %s", argv[a]);
 			myexit(EXIT_FAILURE);
 			}
 		}
/****************************************************************/
/*                                                              */
/*                          Data Preparation                    */
/*                                                              */
/*                                                              */
/****************************************************************/

	prepare_cache (cache);
/****************************************************************/
/*                                                              */
/*                          INPUT SEQ/ALN                       */
/*                                                              */
/*                                                              */
/****************************************************************/


	if ( strm (out_format, "hasch"))
	  {
	    fprintf ( stdout, "%d\n", (int)hash_file(in_file));
	    return EXIT_SUCCESS;
	  }

	if ( rename_file[0])
	  {
	    rename_list=read_rename_file ( rename_file,code);
	  }
	

	if ((D1=read_data_structure (in_format, in_file,RAD))!=NULL)
	  {
	    in_format=(in_format && in_format[0])?in_format:identify_seq_format(in_file);

	    if (print_format)fprintf ( stdout, "\nFILE:%s FORMAT:%s\n", in_file, in_format);
	  }
	else if ( in_file[0])
	        {
		  fprintf ( stdout, "\nFORMAT of file %s Not Supported (1)[FATAL:%s]\n", in_file, PROGRAM);
		  myexit(EXIT_FAILURE);
		}
	

	if ((D2=read_data_structure (in2_format, in2_file,RAD))!=NULL)
	  {
	    if (print_format)
	      fprintf ( stderr, "\nFILE:%s FORMAT:%s\n", in2_file, (in2_format&&in2_format[0])?in2_format:identify_seq_format(in2_file));
	  }
	
	else if (!D2 && in2_file[0])
	        {
		  fprintf ( stderr, "\nFORMAT of file %s Not Supported (2)[FATAL:%s]\n", in2_file, PROGRAM);
		  myexit(EXIT_FAILURE);
		}

/*STRUCTURE INPUT*/


	if ((D_ST=read_data_structure (struc_in_format, struc_in_file,RAD)))
	    {
	     
	   
	      
	      if ( D_ST->CL)
		{
		  Constraint_list *CL;
		  int *entry;

		  CL=D_ST->CL;

		  while ((entry=extract_entry (CL)))
		    {
		      if ( D_ST->S)(D_ST->S)->seq[entry[SEQ1]][entry[R1]-1]=entry[WE];
		    }
		  thread_seq_struc2aln (D_ST->A, D_ST->S);
		}
	      else if ( name_is_in_list ("cons", ((D_ST)->A)->name, ((D_ST)->A)->nseq, 100)!=-1);
	      else
		{
		  
		  D_ST->A=copy_aln ( D1->A,NULL);
		  thread_seq_struc2aln (D_ST->A, D_ST->S);
		}
	    }
	else if ((strcmp (struc_in_format, "rna_number")==0) && in_file[0])
		{
		D_ST->RNA_ST=read_rna_struc_number((D1->A),struc_in_file);
		}
	else if ( struc_in_format[0] && struc_in_file[0])
	        {

		fprintf ( stderr, "\nSTRUC %s UNKNOWN[FATAL]", struc_in_format);
		myexit(EXIT_FAILURE);
		}
	else
	  {
	    D_ST=(Sequence_data_struc*)vcalloc ( 1, sizeof (Sequence_data_struc));
	  }
	
	action=declare_char(argc+1, maxlen+1);
	for ( a=0; a< n_actions;)
	  {
	   if (action_list[a][0]!='+')
	      {
		fprintf ( stderr, "\nWARNING: Action %s Unknown. Actions start with a +", action_list[a]);
		myexit (EXIT_FAILURE);
	      }
	   else
	     {
	     b=0;
	     sprintf ( action[b++], "%s", action_list[a++]+1);
	     while ( a<n_actions && action_list[a][0]!='+')sprintf ( action[b++], "%s", action_list[a++]);
	     modify_data( D1, D2, D_ST, action,b, RAD);
	     }
	  }

	if (rename_list)
	  {
	    if (D1)D1->A= rename_seq_in_aln(D1->A, rename_list);
	    if (D2)D2->A=rename_seq_in_aln (D2->A, rename_list);
	    if (D_ST)D_ST->A=rename_seq_in_aln (D_ST->A,rename_list);

	    if (D1)D1->T  =rename_seq_in_tree (D1->T, rename_list);
	    if (D2)D2->T  =rename_seq_in_tree (D2->T, rename_list);
	    if (D_ST)D_ST->T=rename_seq_in_tree (D_ST->T,rename_list);
	  }


	if ( !out_format[0] && ! struc_out_format[0])sprintf ( out_format, "%s", (in_format && in_format[0])?in_format:"fasta_aln");
	main_output  ( D1, D2, D_ST, out_format, out_file);
	main_output  ( D1, D2, D_ST, struc_out_format, struc_out_file);
	return EXIT_SUCCESS;
	}




/**************************************************************************************************/
/*****************************    FORMAT GUESSING     ******************************************/
/**************************************************************************************************/
Sequence_data_struc *read_data_structure ( char *in_format, char *in_file,	Action_data_struc  *RAD)

        {
	Sequence_data_struc *D;
	char **seq_name=NULL,
	**sequences=NULL;
	Genomic_info *genome_co = NULL;
	int nseq=0, a;


	D=(Sequence_data_struc*)vcalloc ( 1, sizeof (Sequence_data_struc));

	if (!in_file[0])return NULL;
	if (!in_format[0])
	  {
	    in_format=identify_seq_format(in_file);
	  }
	if (!in_format[0])return NULL;
	
	D->A=declare_Alignment(NULL);
	if ( RAD->keep_case)(D->A)->residue_case=KEEP_CASE;

	D->rm_gap=RAD->rm_gap;
	sprintf ( D->format, "%s", in_format);
	

	if (strlen (in_file)>10000)
	  myexit (fprintf_error (stderr,"Pathname exceeds maximum allowd [10000]:\n%s\n",in_file));
	
	
	sprintf ( D->file, "%s", in_file);
	
	
	if (big())return D;

	if ( strm2(in_format,"saga_aln","clustal_aln"))
		{
		main_read_aln (in_file, D->A);
		D->S=aln2seq(D->A);

		}

	else if (strm  (in_format, "nexus"))
	  {
	    D->S=get_nexus(in_file);
	    D->A=seq2aln(D->S, D->A,NO_PAD);
	  }
	else if ( strm (in_format, "treefile_list") || strm(in_format, "treelist"))
	  {

	    D->S=get_treelist(in_file);
	    D->A=seq2aln(D->S, D->A,NO_PAD);
	  }
	else if ( strm (in_format, "file_list") || strm (in_format, "list"))
	  {
	    D->S=get_file_list(in_file);
	    D->A=seq2aln(D->S, D->A,KEEP_GAP);
	  }
	
	else if ( strm (in_format, "fasta_tree"))
	  {

	    D->S=get_fasta_tree (in_file, NULL);
	    D->A=seq2aln(D->S, D->A,NO_PAD);
	    
	  }
	else if ( strm (in_format, "tree_list") || strm (in_format, "treelist"))
	  {
	    D->S=get_treelist (in_file);
	    D->A=seq2aln(D->S, D->A,NO_PAD);
	  }

	else if (strm (in_format, "matrix"))
	  {
	    D->M=read_matrice (in_file);
	  }
	else if (strm (in_format, "mafft_newick_tree"))
	  {
	   
	    D->T=main_read_tree (in_file);
	    D->T=indextree2nametree (D->S, D->T);
	    //D->S=tree2seq(D->T, NULL);
	    D->A=seq2aln (D->S,D->A, 0);
	  }

	else if (strm4 (in_format, "newick_tree", "newick", "nh", "new_hampshire"))
	  {
	   
	    D->T=main_read_tree (in_file);

	    D->S=tree2seq(D->T, NULL);
	    D->A=seq2aln (D->S,D->A, 0);
	  }
	else if (strm (in_format, "blast_aln"))
	  {
	    if (read_blast_aln (in_file, D->A))
	      {
		D->S=aln2seq(D->A);
	      }
	    else
	      {
		return NULL;
	      }
	  }
	else if ( strm( in_format,"number_aln"))
	  {
	    read_number_aln (in_file, D->A);
	    D->S=aln2seq(D->A);
	  }
	else if ( strm( in_format,"stockholm_aln"))
	  {
	    read_stockholm_aln (in_file, D->A);
	    D->S=aln2seq(D->A);
	  }
	else if ( strm( in_format,"gotoh_aln"))
	  {
	    read_gotoh_aln (in_file, D->A);
	    D->S=aln2seq(D->A);
	  }
	
	else if ( strm ( in_format, "msf_aln"))
		{
		read_msf_aln (in_file, D->A);
		D->S=aln2seq(D->A);
		}
	else if ( strm ( in_format, "amps_aln"))
		{
		read_amps_aln (in_file, D->A);
		D->S=aln2seq(D->A);
		}
     	else if ( strm (in_format, "excel_seq"))
		{
		  D->S=perl_reformat2fasta ("excel2fasta.pl",in_file);
		  (D->S)->contains_gap=0;
		  D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}
	else if ( strm (in_format, "pavie_seq"))
		{
		  D->S=perl_reformat2fasta ("pavie2fasta.pl",in_file);
		  (D->S)->contains_gap=0;
		  D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}
	else if ( strncmp (in_format, "perl_",5 )==0)
		{
		  D->S=perl_reformat2fasta (in_format+5,in_file);
		  (D->S)->contains_gap=0;
		  D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}
	else if ( strm (in_format, "number_fasta"))
		{
		D->S=get_fasta_sequence_num (in_file, NULL);
		(D->S)->contains_gap=0;
		D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}
	else if ( strm (in_format, "raw_fasta"))
		{
		D->S=get_fasta_sequence_raw (in_file, NULL);
		(D->S)->contains_gap=0;
		D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}

	else if ( strm2 (in_format, "fasta_aln", "fasta_seq"))
		{

		D->S=get_fasta_sequence (in_file, NULL);
		if ( strcmp (in_format, "fasta_aln")==0)(D->S)->contains_gap=0;
		D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}
	else if ( strm (in_format, "fasta_tree"))
		{

		D->S=get_fasta_tree (in_file, NULL);
		D->A=seq2aln(D->S, D->A, NO_PAD);
		}

	else if ( strm (in_format, "pdb") || strm (in_format, "pdb_struc"))
		{
		    D->S=get_pdb_sequence (in_file);
		    if ( D->S==NULL)
		      {
			add_warning (stderr, "FAILED TO find PDB File %s", in_file);
			myexit (EXIT_FAILURE);
		      }
		    D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}
	else if ( strm2(in_format, "pir_seq", "pir_aln"))
		{
		D->S=get_pir_sequence ( in_file,NULL );
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
        else if ( strm(in_format, "gor_seq") )
		{
		D->S=get_gor_sequence ( in_file,NULL );
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
	else if ( strm2 ( in_format, "dali_aln", "dali_seq"))
		{
		D->S=get_sequence_dali ( in_file);
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
	else if ( strm (in_format, "barton_list_tc"))
		{
		get_barton_list_tc_seq ( in_file);
		}
	else if ( strm (in_format, "amps_sd_scores"))
		{
		D->W=get_amps_sd_scores ( in_file);
		}

	else if ( strm ( in_format, "pima_aln"))
		{
		D->S=get_pima_sequence ( in_file);
		seq2aln (D->S, D->A, RAD->rm_gap);
		}
	else if ( strm( in_format, "gor_struc"))
	        {
		D->S=get_struc_gor ( in_file);
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
	else if ( strm( in_format, "dialign_aln"))
		{
		D->S=get_dialign_sequence ( in_file);
		seq2aln (D->S, D->A, RAD->rm_gap);
		}
	else if ( strm( in_format, "tc_lib") ||  strm( in_format, "mocca_lib") ||  strm( in_format, "lib"))
	        {
		 
		  D->S=read_seq_in_list (in_file);
		  D->CL=declare_constraint_list ( D->S,NULL, NULL, 0,NULL, NULL);
		  D->CL=read_constraint_list_file(D->CL,in_file);
		  seq2aln (D->S, D->A, RAD->rm_gap);
		  free_char (sequences,-1);
		  free_char (seq_name, -1);
		
		}

	else if  (strm (in_format, "alifold"))
	  {
	    D->S=read_alifold ( in_file);
	    seq2aln (D->S, D->A,0);
	  }
	else
	        {
		return NULL;
		}

	if ( D->A)
	  {
	    for ( a=0; a<(D->A)->nseq; a++)sprintf ( (D->A)->file[a], "%s", in_file);
	  }
	if ( D->S)
	  {
	    for ( a=0; a<(D->A)->nseq; a++)sprintf ( (D->S)->file[a], "%s", in_file);
	  }

	return D;
	}
Sequence *read_sequences (char *name)
{
  return main_read_seq (name);
}
Alignment * alifold2aln  (char *file)
{
  Sequence *S;
  S=read_alifold(file);
  sprintf ( S->seq[0],"%s", S->seq[1]);
  return seq2aln (S, NULL, 0);
}
Sequence  * read_alifold (char *file)
{
  Sequence *S;
  char **list;
  int l;
  S=declare_sequence (1,count_n_char_in_file (file),2);
  list=file2lines (file);

  S->seq[0]=list[1];
  S->seq[1]=list[2];
  substitute (S->seq[0], "\n", "\0");
  substitute (S->seq[0], " ", "\0");
  substitute (S->seq[0], "_", STOCKHOLM_STRING);
  l=strlen (S->seq[0]);
  substitute (S->seq[1], "\n", "\0");
  substitute (S->seq[1], " ", "\0");
  substitute (S->seq[1], ".", STOCKHOLM_STRING);
  S->seq[1][l]='\0';
  sprintf (S->name[0], "cons");
  sprintf (S->name[1], "#=GC SS_cons");
  return S;
}






Sequence  * main_read_seq ( char *name)
       {
       char *format=NULL;
       Sequence *S=NULL;
       Alignment *A=NULL;
       int a;
       register_file4dump (name, "r");/*make sure file is included in dump*/
      

       format=identify_seq_format (name);
       
       if ( getenv4debug ("DEBUG_REFORMAT"))fprintf ( stderr, "\n\nFormat %s\n", format);


       if (format &&strm(format, "fasta_seq"))
	 {
	   S= get_fasta_sequence ( name, NULL);
	 }
       else if (format &&strm(format, "pir_seq"))     S= get_pir_sequence ( name, NULL);
       
       else if (format && strstr (format, "aln"))
	 {
	   A=main_read_aln ( name, NULL);
	   S=aln2seq(A);
	   free_aln(A);
	 }
       else if ( format && strstr (format, "tc_lib"))
	 {
	   S=read_seq_in_list (name);
	 }
       else
	  {
	  /*Use The ClustalW routine*/
	    S=cw_read_sequences (name);
	  }

       for ( a=0; a<S->nseq; a++)sprintf ( S->file[a], "%s", name);
       vfree(format);
       ungap_seq(S);
       S=clean_sequence ( S);
       return S;
       }

Sequence  * quick_read_seq ( char *file)
{
  
  register_file4dump (file, "r");
  Sequence *S=NULL;
  char *dup;
  char *tmp_file;

  
  if (!file || !check_file_exists (file))return NULL;
  else if (format_is_fasta (file))
    {
      S=get_fasta_sequence (file, NULL);
    }
  
  else if (printf_system ( "seq2name_seq.pl %s > %s",file, tmp_file=vtmpnam(NULL))==EXIT_SUCCESS)
    {
      
      S=get_fasta_sequence (tmp_file, NULL);
    }
  ungap_seq(S);
    
  if ( (dup=check_hlist_for_dup( S->name, S->nseq)))
    {
      fprintf ( stderr, "ERROR -- Duplicated Sequences %s", dup);
      myexit(fprintf_error (stderr,"ERROR - quick_read_seq: duplicated sequence in %s  ", file));
    }
  
  return S;
}


int *list2sample(int nseq, int sample)
{
  static int *random;
  static int *sampled;
  int a, tot, tried;
  
  if (nseq==0 || sample==0){vfree(random); random=NULL; return random;}
  
  
  if (!random)random=(int*)vcalloc (nseq, sizeof (int));
  else if (read_array_size_new(random)<nseq)random=(int*)vrealloc (random,sizeof (int)*nseq);

  if (!sampled)sampled=(int*)vcalloc (sample, sizeof (int));
  else if (read_array_size_new(sampled)<sample)sampled=(int*)vrealloc (sampled,sizeof (int)*nseq);
  
  
  if (sample>=nseq)
    {
      for (a=0; a<nseq; a++) sampled[a]=1;
      return sampled;
    }
  
  //try first with random sampling
  tot=0;
  tried=0;
  while (tot<sample && tried <100*sample)
    {
      int i=rand()%nseq;
      tried++;
      if (random[i]==0)
	{
	  sampled[tot]=i;
	  random[i]=1;
	  tot++;
	}
    }
  for ( a=0; a<nseq; a++)random[a]=0;
  
  //If failed, random sort and select
  if (tot<sample)
    {
      int **random2=declare_int (nseq, 2);
      for ( a=0; a<nseq; a++)
	{
	  random2[a][0]=a;
	  random2[a][1]=rand()%nseq;
	}
      sort_int (random2, 2, 1, 0, nseq-1);
      for (a=0; a<sample; a++)
	sampled[a]=random2[a][0];
      free_int (random2, -1);
    }
  return sampled;
}
char *FastaRecord2name (char *record)
{
  static char *name;
  char *value;
  char *p;
  if (!record) return NULL;
  value=p=name=csprintf(name, "%s", record);
  value++;
  while (!isspace(p[0]))p++;
  p[0]='\0';
  
   
  return value;
}

int is_fasta_record(char *r, char *file, int action)
{
  if (!r || r[0]!='>')
    if (action==EXIT_FAILURE)
      myexit(fprintf_error (stderr,"Invalid Fasta Record [%s][%s]",(r)?r:"empty", (file)?file:"empty"));
    else
      return 0;
  
  return 1;
}
char *FastaRecord2comment (char *record)
{
  //comment starts after the first space
  static char *name;
  char *value;
  char *p;
  if (!record) return NULL;
  p=name=csprintf(name, "%s", record);
  
  while (!isspace(p[0]))p++;
  value=p;
  if (p[0]=='\n')value[0]='\0';
  else 
    {
      p++;
      value=p;
      while (p[0]!='\n')p++;
      p[0]='\0';
    }
return value;
}
char *FastaRecord2header (char *record)
{
  //header is the first line
  static char *name;
  char *p;
  char *value;
  
  if (!record) return NULL;
  value=p=name=csprintf(name, "%s", record);
  
  while (p[0]!='\n')p++;
  p[0]='\0';
  return value;
}
char *FastaRecord2seq (char *record)
{
  //Sequence starts just after the first line break
  static char *name;
  char *p, *value;
  if (!record) return NULL;
  p=name=csprintf(name, "%s", record);
  
  while (p[0]!='\n')p++;
  return p+1;
}					       
char *seq2clean(char*seq)
{
  int a=0; 
  int b=0;
  
  while (seq[a]!='\0')
    {
      if (isspace(seq[a]));
      else seq[b++]=seq[a];
      a++;
    }
  seq[b]='\0';
  return seq;
}

char*file2record_it (char *file, int i, long *map)
{
  static FILE *fp;
  static char *seq;
  static char *name;
  long start, end;
  int len;
  
  if (!file) 
    {
      if (fp)vfclose (fp);
      fp=NULL;
      return NULL;
    }
    
  if (name!=file)
    {
      if ( fp)vfclose (fp);
      fp=NULL;
      name=file;
    }
  if (!fp)
    {
      fp=vfopen (name, "r");
    }
 

  if (i<0) return NULL;

  start=map[i];
  end=map[i+1];

 
  len=(int)(end-start);
  if (len<=0) return NULL;
  
  if (!seq)seq=(char*)vcalloc (len+1, sizeof (char));
  else if  (read_array_size_new(seq)<=(len+1))seq=(char*)vrealloc (seq,(len+1)*sizeof(char));
  seq[0]='\0';
  
  fseek (fp, start, SEEK_SET);
  fread (seq,sizeof (char), len, fp);
  
  seq[len]='\0';

  return seq;
}

char* file2record (char *file, int i, long *map)
{
  //Safe but ineficient function: causes the file to be closed and opend
  //Only use if files with same name keep being re-cycled
  char *s=file2record_it(file, i, map);
  file2record_it(NULL, i, map);//close file
  return s;
}

char *file2record_old (char *file, int i, long *map)
{
  FILE *fp;
  static char *seq;
  long start, end;
  int len;
  

  if (i<0) return NULL;

  start=map[i];
  end=map[i+1];

 
  len=(int)(end-start);
  if (len<=0) return NULL;
  
  if (!seq)seq=(char*)vcalloc (len+1, sizeof (char));
  else if  (read_array_size_new(seq)<=(len+1))seq=(char*)vrealloc (seq,(len+1)*sizeof(char));
  seq[0]='\0';
  
  fp=vfopen (file, "r");
  fseek (fp, start, SEEK_SET);
  fread (seq,sizeof (char), len, fp);
  vfclose (fp);
  seq[len]='\0';
  return seq;
}

    

  
long *fasta2map(char *file)
{
  FILE *fp;
  int i=0;
  long pos=0;
  char c,lc;
  long *map;
  int ml=1000;
  int a;
  char buf[VERY_LONG_STRING];
  

  
  if ( !file || !check_file_exists (file) || !(fp=vfopen(file, "r"))) return NULL;
  
  map=(long*)vcalloc (ml, sizeof (long));
  
  
  while ((c=fgetc(fp))!=EOF){;}
  vfclose(fp);
  
  fp=vfopen(file, "r");
  while (fgets(buf,VERY_LONG_STRING,fp));
  
  vfclose(fp);  
  

  fp=vfopen(file, "r");
  pos=0;
  while (fgets(buf,VERY_LONG_STRING,fp))
    {
      int d=0;
      while ((c=buf[d++])!='\0')
	{
	  if (c=='>')
	    {
	      if (i>=ml){ml+=VERY_LONG_STRING; map=(long*)vrealloc (map, ml*sizeof (long));}
	      map[i++]=pos;
	    }
	  pos++;
	}
    }
  map=(long*)vrealloc (map, (i+1)*sizeof (long));//map will be used to estimate nseq
  
  map[i]=pos;
  
  vfclose (fp);
  
  return map;
}
    
int fasta2nseq (char *file)
{
  FILE *fp;
  char c;
  int n=0;
  char buf [1000];
  
  if (!(fp=vfopen (file, "r")))return -1;
  
  while (fgets(buf, 1000, fp))
    {
      int d=0;
      while ((c=buf[d++])!='\0')if ( c=='>')n++;
    }	
  vfclose (fp);
  return n;
}
Alignment *reload_aln(Alignment *A)
{
  int a;
  FILE *fp;
  char *tmp=vtmpnam (NULL);
  dump_msa (A,tmp);
  return quick_read_fasta_aln(A,tmp);
}
Alignment * quick_read_fasta_aln (Alignment *A, char *file)
{
  long *map;
  char *s;
  int i, nseq,a;
  file2record_it(NULL,0, NULL);
  map=fasta2map(file);
  nseq=read_array_size_new (map)-1;

  if (!A)A=declare_aln2(nseq,1);
  else if (A->max_n_seq<nseq)A=realloc_aln2(A,nseq,0);
  
  A->nseq=nseq;
  for(i=0;i<nseq; i++)
    {
      s=file2record_it(file,i, map);
      A->seq_al[i]=csprintf (A->seq_al[i], "%s", seq2clean(FastaRecord2seq(s)));
      A->name[i]=csprintf (A->name[i], "%s", FastaRecord2name(s));
      A->aln_comment[i]=csprintf (A->aln_comment[i], "%s", FastaRecord2comment(s));
      A->file[i]=csprintf (A->file[i], "%s", file);
    }

  file2record_it(NULL,0, NULL);
  vfree (map);
  A->declared_len=A->len_aln=strlen (A->seq_al[0]);
  
  return A;
}
			   

Alignment * quick_read_aln ( char *file)
{
  return main_read_aln (file, NULL);
}

Alignment * main_read_aln ( char *name, Alignment *A)
       {
       int a;
       char *dup;
       static char *format;
       Sequence *S=NULL;
       char*tmp_name;
       
       register_file4dump (name, "r");/*make sure file is in dump*/
       


       if ( !name)return NULL;
       else if (!check_file_exists(name))
	 {
	   if ( !check_file_exists (name+1))return NULL;
	   else if ( name[0]=='A') name++;
	   else if ( name[0]=='S') name++;/*Line Added for the -convert flag of T-Coffee*/
	 }
       
       
       if (format_is_fasta(name))//takes care of FASTA
	 {
	   A=quick_read_fasta_aln (A, name);
	   
	 }
       else if   ( printf_system ( "seq2name_seq.pl %s > %s",name, tmp_name=vtmpnam(NULL))==EXIT_SUCCESS)//takes care of common formats
	 {
	   A=quick_read_fasta_aln (A, tmp_name);
	 }
       else//Now go for the unusual
	 {
	   format=identify_seq_format (name);
	   
	   if (format && strm (format, "conc_aln"))A=input_conc_aln (name,NULL);
	   else if (format &&strm(format, "msf_aln"  ))read_msf_aln ( name, A);
	   else if (format &&strm(format, "blast_aln"))read_blast_aln (name, A);
	   else
	     myexit(fprintf_error (stderr,"ERROR - (main_read_aln): unknown format for %s ",name));
	 }
       if (!is_aligned (A->nseq,A->seq_al))
	 {
	   free_aln (A);
	   return NULL;
	 }
       
       if ( (dup=check_hlist_for_dup( A->name, A->nseq)))
          {
	    fprintf ( stderr, "ERROR -- Duplicated Sequences %s", dup);
	    myexit(fprintf_error (stderr,"ERROR - (main_read_aln): duplicated sequence in File %s ", A->file[0]));
	  }

       //Very important---> make sure A stays in sync when A->S provided
       if (!A->S)A->S=ungap_seq(aln2seq(A));
       A=fix_aln_seq(A, A->S);
       compress_aln (A);
       for ( a=0; a< A->nseq; a++) sprintf ( A->file[a], "%s", name);
       A=clean_aln (A);
       return A;
       }
int is_aligned (int nseq, char **seq)
{
  int a,l;
  
  if (!seq|| !nseq)return 0;
  l=strlen (seq[0]);
  
  for (a=1; a<nseq; a++)
    {
      if (l!=strlen (seq[a]))return 0;
    }
  return 1;
}

char * identify_seq_format ( char *file)
       {
       char *format=NULL;
       /*This function identify known sequence and alignmnent formats*/

       
       if (big())
	 {
	   char c;
	   FILE *fp=vfopen (file, "r");
	   while (isspace ((c=fgetc(fp))) && c!=EOF);
	   if ( c!='>')
	     myexit(fprintf_error (stderr,"ERROR - [-big] can only process FASTA files and %s is not a FASTA [FATAL:%s] ",file, PROGRAM));
	   else
	     format=csprintf (format, "fasta_seq");
	   return format;
	 }
       
       
       if ( format==NULL)format=(char*)vcalloc ( 100, sizeof (char));
       else format[0]='\0';
       
       if (1==2)
	 {
	   HERE ("TEST The various formats");
	   
	   HERE ("test 1: %d", is_nameseq (file));;
	   HERE ("test 2: %d", is_stockholm_aln(file));
	   HERE ("test 3: %d", format_is_msf (file));
	   HERE ("test 4: %d", is_blast_file (file));
	   HERE ("test 5: %d", fast_format_determination(file));
	   HERE ("test 6: %d", format_is_oligo    (file));
	   HERE ("test 7: %d", format_is_saga    (file));
	   HERE ("test 8: %d", format_is_conc_aln    (file));
	   HERE ("test 9: %d", is_lib (file));
	   HERE ("test 10: %d", is_lib_02 (file));
	   HERE ("test 11: %d", is_treelist (file));
	   HERE ("test 12: %d", is_newick(file));
	   HERE ("test 13: %d",  is_mafft_newick(file));
	   HERE ("test 14: %d", is_nexus (file));
	   exit (0);
	 }


       int format_val  = 0;
       if ( !check_file_exists(file))
	 {
	   fprintf (stderr, "ERROR: %s Does Not Exist [FATAL:%s]\n",file, PROGRAM);
	   myexit (EXIT_FAILURE);
	 }
       else
	 
	 
	 // 	if ( format_is_fasta_seq(file))sprintf ( format, "fasta_seq");
	 // 	else if ( format_is_fasta_aln(file,1))
	 // 	sprintf ( format, "fasta_aln");
	 if ( is_newick(file))sprintf ( format, "newick_tree");
	 else if ( is_mafft_newick(file))sprintf ( format, "mafft_newick_tree");
	 else if ( is_stockholm_aln (file))sprintf (format, "stockholm_aln");
         else if (format_is_msf (file))sprintf (format, "msf_aln");
	 else if ( is_blast_file (file))sprintf ( format, "blast_aln");
	 else if ( is_pdb_file(file))sprintf ( format, "pdb_struc");
	 else if (format_val = fast_format_determination(file))
	   {
	     
	     switch( format_val )
	       {
	       case 1 : sprintf ( format, "fasta_seq");
		 break;
	       case 2 : sprintf ( format, "fasta_aln");
		 break;
	       case 3 : sprintf ( format, "pir_seq");
		 break;
	       case 4 : sprintf ( format, "pir_aln");
		 break;
	       }
	   }
	 else if ( format_is_oligo    (file))sprintf ( format, "oligo_aln");
	
	 else if ( format_is_saga     (file))sprintf ( format, "clustal_aln");
	 else if ( format_is_conc_aln (file))sprintf ( format, "conc_aln");
	 else if ( is_lib (file))sprintf ( format, "tc_lib");
	 else if ( is_lib_02 (file))sprintf ( format, "tc_lib_02");
	 else if ( is_treelist (file))sprintf ( format, "treelist");
	 
         else if ( is_nexus (file))sprintf ( format, "nexus");
       
       
	 else
	   {
	     //add_warning ( stderr, "\nThe Format of File: %s was not recognized [SERIOUS:%s]",file, PROGRAM);
	     ;
	   }
       // 	 printf("DETERMINED FORMAT: %s\n", format);
      
       return format;
       }
char **identify_list_format ( char **list, int n)
       {
	   int a;
	   char *name;
	   char *string;
	   char mode;



	   declare_name (name);
	   for ( a=0; a< n; a++)
	       {

		 sprintf (name, "%s", list[a]);
		 string=list[a];
		 if ((mode=identify_format ( &string))!='?')
		   {
		       sprintf ( name, "%s", string);
		       sprintf ( list[a], "%c%s", mode,name);
		   }
	       else
	           {
		       fprintf ( stderr, "\nERROR: %s not recognised [FATAL:%s]", name, PROGRAM);
		   }

	       }

	   vfree(name);
	   return list;
       }

char * name2type_name ( char *name)
{
	/*turns <file> into <Sfile>, <Afile>...*/
	char *new_name;
	char mode;
	


	new_name=(char*)vcalloc ( strlen (name)+2, sizeof (char));
	sprintf ( new_name, "%s", name);

	if (is_in_set (name[0], "ALSMXPRW") && !check_file_exists(name))
	  {
		sprintf ( new_name, "%s", name);
	  }
	else
	  {
	    mode=identify_format (&new_name);
	    sprintf ( new_name, "%c%s", mode,name);
	    
	  }

	return new_name;
}

char identify_format (char **fname)
       {
	   char mode='?';
	   mode=fname[0][0];
	   
	   if ((is_in_set (mode, "ALMSPR") && check_file_exists(fname[0]+1)) ||(mode=='X' && is_matrix ( fname[0]+1)) ||(mode=='M' && is_method(fname[0]+1)) )
	     {

	       fname[0]++;
	     }
	   else if (mode=='W' && !check_file_exists(fname[0])){fname[0]++;}
	   else
	       {

		 /*WARNING: Order matters => internal methods can be confused with files, must be checked last*/
                      if (is_lib(fname[0]))mode='L';
		      else if (is_pdb_file(fname[0]))mode='P';
		      else if (is_seq(fname[0]))mode='S';
		      else if (is_aln(fname[0]))mode='A';
		      else if (is_matrix(fname[0]))mode='X';
		      else if (is_method(fname[0]))mode='M';
		      else mode='?';
		  }
	   return mode;
       }


int is_pdb_name_new ( char *name)
    {
      //for some reason this crashes
      int result;
      
      static char **buf_names;
      static int   *buf_result;
      static int   nbuf;
      static int maxnbuf;
      
      /*Use the look up*/
      if (!buf_names)
	{
	   maxnbuf+=10;
	   buf_names =(char**)vcalloc (maxnbuf,sizeof (char*)); 
	   buf_result=(int*  )vcalloc (maxnbuf,sizeof (int)); 
	}
      else if (nbuf==maxnbuf)
	{
	  maxnbuf+=1000;
	  buf_names =(char**)vrealloc (buf_names ,maxnbuf*sizeof (char*)); 
	  buf_result=(int*  )vrealloc (buf_result,maxnbuf*sizeof (int)); 
	}
      if ((result=name_is_in_list ( name, buf_names,nbuf,100))!=-1)result=buf_result[result];
      else 
	{
	  char *s=printf_system2string ("extract_from_pdb -is_pdb_name \'%s\'", name);

	  if (!s)result=0;
	  else
	    {
	      
	      substitute (s, "\n", "");
	      substitute (s, "\n", "");
	     
	      result=atoi (s);
	      
	      buf_names[nbuf]=csprintf ( buf_names[nbuf], "%s", name);
	      result=buf_result[nbuf++]=(result==1)?1:0;
	      hupdate(buf_names);
	    }
	}
      return result;
    }

int is_pdb_name ( char *name)
    {
      char command[1000];
      int result;
      char *result_file;
      static char **buf_names;
      static int   *buf_result;
      static int   nbuf;
      FILE *fp;


      /*Use the look up*/
      if ( !buf_names)
	{
	  buf_names=declare_char (10000, 100);
	  buf_result=(int*)vcalloc (10000, sizeof (int));
	}
      if ((result=name_is_in_list ( name, buf_names,nbuf,100))!=-1)return buf_result[result];

      

      result_file=vtmpnam (NULL);

      sprintf ( command, "extract_from_pdb -is_pdb_name \'%s\' > %s", name, result_file);
      if ( getenv4debug ("DEBUG_EXTRACT_FROM_PDB"))fprintf ( stderr, "\n[DEBUG_EXTRACT_FROM_PDB:is_pdb_name] %s\n", command);
      my_system ( command);

      fp=vfopen ( result_file, "r");
      fscanf ( fp, "%d", &result);
      vfclose (fp);
      vremove ( result_file);

      sprintf ( buf_names[nbuf], "%s", name);
      result=buf_result[nbuf++]=(result==1)?1:0;
      hupdate(buf_names);
      return result;

    }

char*  get_pdb_id ( char *file)
{
  /*receives the name of a pdb file*/
  /*reads the structure id in the header*/
  /*returns the pdb_id*/
  char *tmp_name;
  char command[10000];
  char cached [1000];
  char fname[1000];
  FILE *fp;
  char *id;
  char buf[1000];


  tmp_name=vtmpnam(NULL);

  sprintf ( cached, "%s/%s", get_cache_dir(),file);
  if ( check_file_exists(cached))sprintf ( fname, "%s", cached);
  else sprintf ( fname, "%s", file);

  sprintf ( command, "extract_from_pdb -get_pdb_id %s > %s",fname, tmp_name);

  if ( getenv4debug ("DEBUG_EXTRACT_FROM_PDB"))fprintf ( stderr, "\n[DEBUG_EXTRACT_FROM_PDB:get_pdb_id] %s\n", command);
  my_system ( command);

  buf[0]='\0';
  fp=vfopen (tmp_name, "r");
  fscanf ( fp, "\n%s\n", buf);
  vfclose (fp);

  if ( getenv4debug ("DEBUG_EXTRACT_FROM_PDB"))fprintf ( stderr, "\n[DEBUG_EXTRACT_FROM_PDB:get_pdb_id]DONE\n");

  id=(char*)vcalloc ( strlen (buf)+1, sizeof (char));
  sprintf ( id, "%s", buf);



  return id;
}


char*  get_pdb_struc(char *in_name, int start, int end)
    {
      char *name1,*name2;
      char command[LONG_STRING];
      char *name;




      name=(char*)vcalloc ( STRING, sizeof (char));
      sprintf ( name, "%s", in_name);

      if ( (name1=is_pdb_struc(name))==NULL && (name[0]=='P' && ((name1=is_pdb_struc (name+1))==NULL)))
	{
	  fprintf ( stderr, "\nERROR Could not download structure %s [FATAL:%s]\n", name, PROGRAM);crash("");
	}
      else if ( (start==0) && (end==0))return name1;
      else
	{
	  declare_name(name2);
	  sprintf ( name2, "%s_%d_%d.pdb", name, start, end);
	  sprintf ( command, "extract_from_pdb -infile \'%s\' -chain FIRST -coor %d %d > %s%s",check_file_exists(name1),start, end, get_cache_dir(),name2);
	  if ( getenv4debug ("DEBUG_EXTRACT_FROM_PDB"))fprintf ( stderr, "\n[DEBUG_EXTRACT_FROM_PDB:get_pdb_struc] %s\n", command);
	  my_system (command);

	  if ( is_pdb_file(name2))return name2;
	  else
	    {
	      fprintf ( stderr, "\nERROR Could not extract segment [%d %d] from structure %s [FATAL:%s]\n",start, end, name, PROGRAM);crash("");
	    }
	  myexit (EXIT_FAILURE);
	}

      return NULL;
    }

char*  seq_is_pdb_struc ( Sequence *S, int i)
{

  if (!S){return NULL;}
  else if ( !S->T[i]){return NULL;}
  else if ( !((S->T[i])->P)){return NULL;}
  else return ((S->T[i])->P)->template_file;
}
char*  is_pdb_struc_strict ( char *name);
char*  is_pdb_struc ( char *iname)
{
  static char **name=(char**)vcalloc(10, sizeof (char*));
  char *r=NULL;
  int a, i;
  
  a=0;
  name[a]=csprintf ( name[a], "%s", iname);a++;
  name[a]=csprintf ( name[a], "%s.pdb", iname);a++;
  if (getenv ("PDB_DIR"))name[a]=csprintf ( name[a], "%s/%s", getenv("PBD_DIE"),iname);a++;
  if (getenv ("PDB_DIR"))name[a]=csprintf ( name[a], "%s/%s.pdb",getenv("PDB_DIR"),iname);a++;
  if (get_cache_dir())name[a]=csprintf ( name[a], "%s/%s", get_cache_dir(),iname);a++;
  if (get_cache_dir())name[a]=csprintf ( name[a], "%s/%s.pdb", get_cache_dir(),iname);a++;
  
  for (i=0; i<a;i++)
    {
      name[i]=substitute(name[i], " ", "");
      if ((r=is_pdb_struc_strict (name[i]))!=NULL)return r;
    }
  add_warning(stderr, "%s Could not be used to find a PDB template",iname);
  return r;
}
char*  is_pdb_struc_strict ( char *name)
   {
     char *r=NULL;
     
     if ( !name || name[0]=='\0')return NULL;
     
     if ((r=get_string_variable(name))!=NULL)return csprintf(NULL, "%s",r);
          
     if (is_pdb_file(name)){r=name;}
     else if (is_pdb_name (name))
       {
	 printf_system ("extract_from_pdb -netfile \'%s\' > %s/%s 2>/dev/null",name, get_cache_dir(),name);
	 if ( is_pdb_file(name))r=name;
	 else r=NULL;
	 
       }

     if (r)
       {
	 set_string_variable (name, r);
	 return csprintf (NULL, "%s", r);
       }
     else return NULL;
   }
char*  is_pdb_struc_strict_old2 ( char *name)
   {
     /*Receives a name
       checks if this is the name of a local file that contains PDB data
       checks if this is the name of a file from a local db
                                            put the file in the cache
       checks if this is a file from a remote db (extract_from_pdb
       return NULL if everything fails
     */

     static char **buf_names;
     static char **buf_result;
     static int   nbuf, s;
     static int max_nbuf;

     char *r=NULL;
     char command[1000];
     char *rvalue;

     if ( !name || name[0]=='\0')return NULL;


     /*Use the look up*/
     if ( !buf_names)
	{

	  buf_names=(char**)vcalloc ( 1000, sizeof (char*));
	  buf_result=(char**)vcalloc ( 1000, sizeof (char*));
	  max_nbuf=1000;
	}

     if (nbuf>=max_nbuf)
       {
	 int i;
	 buf_names= (char**)vrealloc ( buf_names, (max_nbuf+1000)*sizeof (char*));
	 buf_result=(char**)vrealloc ( buf_result,(max_nbuf+1000)*sizeof (char*));
	 
	 for (i=max_nbuf;i<max_nbuf+1000; i++)
	   {
	     buf_names[i]=NULL;
	     buf_result[i]=NULL;
	   }
	 max_nbuf+=1000;
       }
     
     if ( (s=name_is_in_list ( name, buf_names,nbuf,-1))!=-1)return buf_result[s];
     

     r=NULL;
     



     if (is_pdb_file(name)){r=name;}
     else if (is_pdb_name (name))
       {
	 printf_system ("extract_from_pdb -netfile \'%s\' > %s/%s 2>/dev/null",name, get_cache_dir(),name);
	 if ( is_pdb_file(name))r=name;
	 else r=NULL;
	 
       }
    

      /*Fill the buffer*/
     buf_names[nbuf]=(char*)vcalloc ( strlen (name)+1, sizeof (char));
     sprintf ( buf_names[nbuf], "%s", name);
     if ( r)
       {
	 buf_result[nbuf]=(char*)vcalloc ( strlen (r)+1, sizeof (char));
	 sprintf (buf_result[nbuf], "%s", r);
       }
     else buf_result[nbuf]=NULL;
     nbuf++;
     return r?csprintf (NULL, "%s", r):NULL;
     
   }

char *atom2pdbF (char*in)
{
  char *tmp=vtmpnam(NULL);
  printf_system ("extract_from_pdb %s > %s", in,tmp);
  return tmp;
}
char *fix_pdb_file ( char *in1)
{
  char*in=NULL;
  char*tmp=NULL;
  
  
  
  //1 - Make sure we have the proper file
  //    The file may be in PDB_DIR, 
  if (check_file_exists (in1))//This will get the cached name if any
    {
      in=csprintf (in,"%s",check_file_exists (in1));
    }
  else if (getenv ("PDB_DIR") && !check_file_exists (in1))//look at provided address first and THEN on PDB_DIR
    {
      in=csprintf (in,"%s/%s",getenv("PDB_DIR"), in1);
      if ( check_file_exists (in));
      else in=csprintf (in,"%s",in1);
    }
  else in=csprintf (in,"%s", in1);
  

  


  if ( !in || !check_file_exists (in))
    {
      tmp=NULL;
    }
  else if (!is_pdb_file(in))
    {
      if (!pdb_has_atom(in))
	{
	  add_warning ( stderr, "File %s does not contain any ATOM field [%s:WARNING]",in, PROGRAM);
	  tmp=NULL;
	}
      else if (is_pdb_file (tmp=atom2pdbF(in)))
	add_warning ( stderr, "File %s is a partial PDB File - regenerated from ATOM with extract_from_pdb [%s:WARNING]",in, PROGRAM);
      
    }
  else if ( !seqres_equal_atom(in))
    {
      tmp=atom2pdbF(in);
      if (!is_pdb_file (tmp))return NULL;
      else
	  add_warning ( stderr, "Could not use SEQRES field in %s -- used ATOM instead [%s:WARNING]",in, PROGRAM);
    }
  else
    {
      tmp=csprintf (tmp, "%s", in);
    }
  vfree (in);
  return tmp;
}

int is_sap_file ( char *name)
	{
	FILE *fp;
	if (!name)return 0;
	if (!check_file_exists(name))return 0;
	
	if ((fp=find_token_in_file (name, NULL, "Percent"))!=NULL)
	  {
	    if ((fp=find_token_in_file (name,fp, "Percent"))!=NULL)
	      {
		vfclose (fp);
             	return 1;
	      }
	    else
	      {
		return 0;
	      }
	  }
	else
	  {
	    return 0;
	  }
	}


int is_blast_file ( char *name)
       {
	 if ( !check_file_exists(name) ) return 0;
	 else if (token_is_in_file (name, "<SequenceSimilaritySearchResult>"))
	   {
	     return BLAST_XML;
	   }
	 else
	   {
	     if (token_is_in_file (name, "Lambda") && token_is_in_file (name, "Altschul,"))
	       {
		 return BLAST_TXT;
	       }
	     else
	       {
		 return 0;
	       }
	   }
	 return 0;
       }
int is_simple_pdb_file ( char *name)
{
  FILE *fp;
  if ((fp=find_token_in_file (name, NULL, "SIMPLE_PDB_FORMAT"))!=NULL){vfclose (fp);return 1;}
  return 0;
}

/**
 * \brief Checks if a given filename belongs to a PDB.
 * \param fname The filename
 * \return 1 if the file exists and the keywords "HEADER", "SEQRES" and "ATOM" appear (in this order) at the beginning of a line. 0 otherwise
 */
int pdb_has_atom ( char *name)
       {
	 FILE *fp;
	 int ispdb=0;

	 if ( name==NULL) return 0;
	 if (!check_file_exists (name))return 0;

	

	 if ((fp=find_token_in_file (name, NULL, "\nATOM"))!=NULL)
	   {
	     vfclose (fp);
	     ispdb++;
	   }
	 else
	   {
	     ispdb=0;
	   }
	 return ispdb;
       }
int is_pdb_file ( char *name)
{
  //everything except ATOM can be regenerated
  return pdb_has_atom (name);
}
int is_pdb_file_old ( char *name)
       {
	 FILE *fp;
	 int ispdb=0;
	 int *a;
	 if ( name==NULL) return 0;
	 if (!check_file_exists (name))return 0;

	 if ((fp=find_token_in_file (name, NULL, "\nHEADER"))!=NULL)
           {

	     vfclose (fp);
	     ispdb++;
	   }

	 
	 if ((fp=find_token_in_file (name, NULL, "\nSEQRES"))!=NULL)
           {

	     vfclose (fp);
	     ispdb++;
	   }

	 if ((fp=find_token_in_file (name, NULL, "\nATOM"))!=NULL)
	   {
	     vfclose (fp);
	     ispdb++;
	   }
	 else
	   {
	     ispdb=0;
	   }


	 if ( ispdb>=2)return 1;
	 else return 0;
       }


// int is_pdb_file ( char *fname)
//        {
// 	 FILE *fp=NULL;
// 	 int ispdb=0;
//
// 	 if ( fname==NULL) return 0;
// 	 if (!check_file_exists (fname))return 0;
//
// // 	static char *name;
// 	int token_len;
//
// 	int only_start;
// 	const int LINE_LENGTH=1000;
// 	char line[LINE_LENGTH];
//
//
// 	/*Note: Token: any string
// 	If Token[0]=='\n' Then Token only from the beginning of the line
// 	*/
//
// 	if (!fp && !file_exists("CACHE",fname))
// 		return NULL;
//
// 	if (!fp)
// 	{
// 		fp=vfopen ( fname, "r");
// 	}
//
// 	line[LINE_LENGTH-2]='\0';
// 	while (fgets(line, LINE_LENGTH, fp)!=NULL)
// 	{
// 		if ((line[LINE_LENGTH-2]!='\0') && (line[LINE_LENGTH-2]!='\n'))
// 			line[LINE_LENGTH-2]='\0';
// 		else
// 		{
// 			if (!strncmp(line,"HEADER",6))
// 			{
// 				++ispdb;
// 				break;
// 			}
// 		}
// 	}
// 	while (fgets(line, LINE_LENGTH, fp)!=NULL)
// 	{
// 		if ((line[LINE_LENGTH-2]!='\0') && (line[LINE_LENGTH-2]!='\n'))
// 			line[LINE_LENGTH-2]='\0';
// 		else
// 		{
// 			if (!strncmp(line,"SEQRES",6))
// 			{
// 				++ispdb;
// 				break;
// 			}
// 		}
// 	}
// 	while (fgets(line, LINE_LENGTH, fp)!=NULL)
// 	{
// 		if ((line[LINE_LENGTH-2]!='\0') && (line[LINE_LENGTH-2]!='\n'))
// 			line[LINE_LENGTH-2]='\0';
// 		else
// 		{
// 			if (!strncmp(line,"ATOM",4))
// 			{
// 				++ispdb;
// 				break;
// 			}
// 		}
// 	}
//
// 	return (ispdb>=2?1:0);
// }


int is_seq ( char *name)
       {
	 char *format;

	 if ( !check_file_exists(name))return 0;

	 format= identify_seq_format(name);
	 if(!format || format[0]=='\0'){vfree (format);return 0;}
	 else if (strstr(format, "seq")){vfree (format);return 1;}
	 else return 0;
       }
int is_aln ( char *name)
       {
       char *format;
       if ( !check_file_exists       (name))return 0;
       
       format= identify_seq_format(name);
       if ( !format || format[0]=='\0'){vfree (format);return 0;}
       else if (strstr(format, "aln")){vfree (format); return 1;}
       else return 0;
       }

int is_matrix (char *name)
       {
       int **m;
       if (!name) return 0;
       if ((m=read_matrice (name))!=NULL){free_int (m, -1); return 1;}
       return 0;
       }
int is_treelist (char *name)
{
  int c;
  FILE *fp;
  static char *buf=NULL;
  int cont=1;
  int n=0;
  fp=vfopen (name, "r");
  while ((buf=vfgets(buf, fp)) && cont && n<2)
    {
      cont=0;
      if (buf[0]=='(')
	{
	  int l=strlen (buf)-1;
      	  while (isspace(buf[l]))l--;
	  if (buf[l]==';'){n++;cont=1;}
	}
    }
  vfclose (fp);
  if (n>1) return 1;
  return 0;
}
int is_mafft_newick (char *name)
   {
     if (file2firstchar (name)=='(' && file2lastchar (name)!=';')return 1;
     return 0;
   }
int is_newick (char *name)
   {
     if (file2firstchar (name)=='(' && file2lastchar (name)==';')return 1;
     return 0;
   }
int is_clustalw_matrix ( char *name)
{

  FILE *fp;


       if ( (fp=find_token_in_file (name, NULL, "CLUSTALW_MATRIX"))!=NULL){vfclose(fp);return 1;}
       else return 0;
}
int is_pavie_matrix ( char *name)
{

  FILE *fp;


       if ( (fp=find_token_in_file (name, NULL, "PAVIE_MATRIX"))!=NULL){vfclose(fp);return 1;}
       else return 0;
}
int is_distance_matrix_file (char *name)
{
  FILE *fp;
  if ( (fp=find_token_in_file (name, NULL, "TC_DISTANCE_MATRIX_FORMAT_01"))!=NULL){vfclose(fp);return 1;}
  else return 0;
}
int is_similarity_matrix_file (char *name)
{
  FILE *fp;
  if ( (fp=find_token_in_file (name, NULL, "TC_SIMILARITY_MATRIX_FORMAT_01"))!=NULL){vfclose(fp);return 1;}
  else return 0;
}

int is_blast_matrix ( char *name)
{

  FILE *fp;
  int r=0;
  if (check_file_for_token   (name,"BLAST_MATRIX"))r=1;
  else if (check_file_for_token(name,"BZX*"))r=1;
  else if (check_file_for_token(name,"B  Z  X  *"))r=1;
  else if (check_file_for_token(name,"B   Z   X   *"))r=1;
  else r=0;
  return r;

}
int is_single_seq_weight_file ( char *name)
{

  
  return token_is_in_file ( name, "SINGLE_SEQ_WEIGHT_FORMAT_01");

}
int is_nameseq (char *file)
{
  return token_is_in_n_lines (file,"#NAMESEQ_01",2);
}
int is_nexus (char *file)
{
  return token_is_in_n_lines (file,"#NEXUS",10);
  
}
int is_stockholm_aln (char *file)
{
  
  return token_is_in_n_lines (file,"STOCKHOLM",2);
}

int is_lib ( char *name)
{
  return is_lib_01(name);
}

int is_lib_02 ( char *name)
{

  return token_is_in_file ( name, "TC_LIB_FORMAT_02");

}

int is_lib_01 (char *name)
       {


	 if ( token_is_in_file ( name, "TC_LIB_FORMAT_01")) return 1;
	 else if (token_is_in_file ( name, "T-COFFEE_LIB_FORMAT_01"))return 1;
	 else if (token_is_in_file (name, "SEQ_1_TO_N"))return 1;
	 else return 0;
       }

int is_lib_list ( char *name)
{
  if ( !check_file_exists (name))return 0;
  if ( token_is_in_file ( name, "TC_LIB_LIST_FORMAT_01")) return 1;
  return 0;
}
int is_method ( char *file)
    {
	char new_file[200];


	sprintf ( new_file, "%s", file);
	if ( (token_is_in_file(new_file, "TC_METHOD_FORMAT_01"))){return 1;}
	if ( is_in_pre_set_method_list(new_file))
	    {

		vremove ( new_file);
		return 1;
	    }
	else
	  {

	    return 0;
	  }
    }

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              SEQUENCE FORMAT IDENTIFIERS                                */
/*                                                                                         */
/***************************************************************************************** */
int type_is_exon_boundaries(char **seq, int n)
{
  int a, l, b;
  for (a=0; a<n; a++)
    {
      l=strlen (seq[a]);
      for (b=0; b<l; b++)
	if ( strchr ("bojBOJ", seq[a][b]))return 1;
    }
  return 0;
}

int format_is_oligo(char *file)
    {
    char *buf=NULL;
    FILE *fp;
    int r=0;

    fp=vfopen ( file, "r");
    buf=vfgets(buf,fp);
    vfclose ( fp);


    if ( strm (buf, "ALPHABET"))r=1;

    vfree (buf);

    return r;
    }
int format_is_msf ( char *file)
    {
      
      
      if (!token_is_in_n_lines (file,"MSF:",5)) return 0;
      else if (!token_is_in_n_lines (file,"Type:",5))return 0;
      else return 1;

    }

//Fasta and PIR
int format_is_not_fasta (char *name)
{
  FILE *fp;
  char c;
  
  fp=vfopen (name, "r");
  c=fgetc(fp);
  vfclose (fp);
  if ( c!='>')return 1;
  else return 0;
}

int format_is_fasta_aln ( char *file, int i_know_that_it_not_seq)

{
	if( i_know_that_it_not_seq && format_is_fasta(file))
		return 1;
	else
		if ( format_is_fasta(file) && !format_is_fasta_seq(file))
			return 1;
	return 0;
}


/*
 * checks if the file is in fasta or pir format
 * returns:
 *   0 nothing of these
 *   1 fasta sequences
 *   2 fasta alignemnt
 *   3 pir sequences
 *   4 pir alignment
 */
int fast_format_determination  ( char *in_f)
{
	const unsigned int READ_LENGTH = 401;
	char line[READ_LENGTH];
	FILE *in_F = fopen(in_f, "r");
	int is_pir = 1;
	int is_aln = 1;
	char *tmp;
	if (fgetc(in_F) != '>')
	{
		fclose(in_F);
		return 0;
	}
	ungetc('>',in_F);
	int last_length = -1;
	int current_length = -1;
	unsigned int n_seqs = 0;
	char last='/';
	int has_gap = 0;
	int line_len;
	while (fgets(line, READ_LENGTH, in_F) != NULL)
	{
		if (line[0] == '>')
		{
			++n_seqs;
			if (!is_pir_name(&line[1]))
				is_pir = 0;
			if ((last_length != -1))
			{
				if (last != '*')
					is_pir = 0;
				if (current_length != last_length)
					is_aln = 0;
			}
			last_length = current_length;
			current_length = 0;
		}
		else
		{
			tmp = &(line[0]);
			while (*tmp != '\0')
			{
				if ((*tmp != '\n') && (*tmp != ' '))
				{
					if ((*tmp == '-') || (*tmp == '.')|| (*tmp == '*')|| (*tmp == '#')|| (*tmp == '~'))
						has_gap = 1;
					++current_length;
					last = *tmp;
				}
				++tmp;
			}
		}
	}

	if (has_gap)
		is_aln = 1;
	else
		is_aln = 0;
	if (last != '*')
		is_pir = 0;


	fclose(in_F);
	if (!is_pir)
	{
		if (is_aln && (n_seqs >1))
			return 2;
		return 1;
	}
	if (is_pir)
	{
		if (is_aln && (n_seqs >1))
			return 4;
		return 3;
	}
}


int format_is_fasta_seq  ( char *file)
    {
      int a, l1, l2,l;
      Sequence *S;

      if ( format_is_fasta (file))
	{
	S=get_fasta_sequence (file, NULL);
	if (!S) return 0;
	else if ( !S->seq[0]){free_sequence (S, S->nseq); return 1;}
	l=strlen ( S->seq[0]);
	for ( a=0; a< S->nseq; a++)if(strlen(S->seq[a])!=l){free_sequence (S, S->nseq);return 1;}
	for ( a=0; a< S->nseq; a++)
	  {
	    l1=strlen ( S->seq[a]);
	    ungap (S->seq[a]);
	    l2=strlen ( S->seq[a]);
	    if ( l1!=l2)
	      {
		free_sequence (S, S->nseq);
		return 0;
	      }
	  }
	free_sequence (S, S->nseq);
	return 1;
      }
    else
      {
	return 0;
      }
    }

int format_is_fasta ( char *file)
    {
      if ( !isfile(file))return 0;
      if ( get_first_non_white_char (file)=='>')return 1;
      return 0;
    }

int format_is_pir_aln ( char *file)
    {
      if ( format_is_pir(file) && !format_is_pir_seq(file))return 1;
      else return 0;
    }

int format_is_pir_seq ( char *file)
    {
      int a, l1, l2;
      Sequence *S;


    if ( format_is_pir (file))
      {
	S=get_pir_sequence (file, NULL);
	for ( a=0; a< S->nseq; a++)
	  {
	    l1=strlen ( S->seq[a]);
	    ungap (S->seq[a]);
	    l2=strlen ( S->seq[a]);
	    if ( l1!=l2)
	      {
		free_sequence (S, S->nseq);
		return 0;
	      }
	  }
	return 1;
      }
    else
      {
	return 0;
      }
    }


int format_is_pir ( char *file)
    {
      Sequence *S;
      int pir_name=1, star_end=1, a;

      S=get_fasta_sequence (file, NULL);
      if (!S)return 0;
      else if (!S->seq[0])return 0;

      pir_name=1; star_end=1;
      for (a=0; a< S->nseq; a++)
	{
	  int l;
	  if (!is_pir_name(S->name[a]))pir_name=0;
	  l=strlen (S->seq[a]);
	  if (!l || (l && S->seq[a][l-1]!='*'))
	    star_end=0;
	}
      free_sequence(S,-1);
      if ( pir_name && star_end) return 1;
      else return 0;
    }


int is_pir_name (char *name)
{
	if (name[2] != ';')
		return 0;
  if ( strstr (name, "P1;"))return 1;
  if ( strstr (name, "F1;"))return 1;
  if ( strstr (name, "DL;"))return 1;
  if ( strstr (name, "DC;"))return 1;
  if ( strstr (name, "RL;"))return 1;
  if ( strstr (name, "RC;"))return 1;
  if ( strstr (name, "XX;"))return 1;
  return 0;
}


int format_is_conc_aln (char *file)
{
  FILE *fp;
  if ( (fp=find_token_in_file (file, NULL, "CONC_MSF_FORMAT_01"))){vfclose (fp); return 1;}
  return 0;
}
int format_is_saga ( char *file)
    {
    
    if (token_is_in_file_n(file, "SAGA", 2))return 1;
    else if (token_is_in_file_n(file, "CLUSTAL", 2))return 1;
    else if (token_is_in_file_n(file, "CLUSTALW", 2))return 1;
    else if (token_is_in_file_n(file, "ClustalW", 2))return 1;
    else if (token_is_in_file_n(file, "clustalw", 2))return 1;
    else if (token_is_in_file_n(file, "T-COFFEE_MSA", 2))return 1;
    else if (token_is_in_file_n(file, "INTERLEAVED_MSA", 2))return 1;
    else return 0;
    
    }



/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT STUFF                                              */
/*                                                                                         */
/***************************************************************************************** */
int output_format_aln ( char *format, Alignment *inA, Alignment *inEA,char *name)
        {
	Sequence_data_struc *D1=NULL;
	Sequence_data_struc *D2=NULL;
	Alignment *A=NULL;
	Alignment *EA=NULL;


	A =copy_aln (inA, NULL);
	A->CL=inA->CL;
	EA=copy_aln (inEA,NULL);
	A =expand_aln(A);
	EA=expand_number_aln(inA,EA);

	
	

	if (A && A->expanded_order  )A=reorder_aln ( A, A->expanded_order,A->nseq);
	if (EA && EA->expanded_order)EA=reorder_aln ( EA, EA->expanded_order,EA->nseq);


        D1=(Sequence_data_struc*)vcalloc ( 1, sizeof (Sequence_data_struc));
	D1->A=A;
	if (EA)
	   {
	   D2=(Sequence_data_struc*)vcalloc ( 1, sizeof (Sequence_data_struc));
	   D2->A=EA;
	   }

	main_output ( D1, NULL,D2, format, name);

	vfree(D1);
	vfree(D2);
	free_aln (A);
	free_aln (EA);
	return 1;
	}
int main_output  (Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char *out_format, char *out_file)

	{
	FILE *fp;
	int value;
	Alignment *BUF_A;
	int expanded=0;
	if ( big())
	  {
	    
	    if (D1->file[0]=='\0') return EXIT_SUCCESS;
	    
	    if (!out_file) fp=stdout;
	    else if ( strm (out_file, "stdout"))fp=stdout;
	    else if ( strm (out_file, "stderr"))fp=stderr;
	    else fp=vfopen (out_file, "w");
	    
	    vfclose (display_file_content (fp,D1->file));
	    D1->file[0]='\0';
	    return EXIT_SUCCESS;
	  }
	if ( !out_format[0])return 0;
	if ( D1 && D1->rm_gap)ungap_aln ((D1->A));

	if ( (strstr (out_format, "expanded_")))
	  {
	    if (!D1) return 1;
	    out_format+=strlen ("expanded_");
	    BUF_A=copy_aln (D1->A, NULL);
	    (D1->A)=thread_profile_files2aln ((D1->A), NULL, NULL);
	    expanded=1;
	  }

	if ( strm (out_format, "") || strm (out_format, "no") || strm (out_format, "NO"))return 0;
	else if ( strm (out_format, "sp_ascii"))
	  {
	    if (!D1) return 1;
	    sp_triplet_coffee_evaluate_output (D1->A, (D1->A)->CL, out_file);
	  }
	else if ( strm (out_format, "sp_lib"))
	  {
	    if (!D1) return 1;
	    sp_triplet_coffee_evaluate_output2 (D1->A, (D1->A)->CL, out_file);
	  }
	else if (    ( strm (out_format, "aln2lib")))
	  {
	    int a, b, c;
	    int r1,r2,s1, s2,s;
	    Constraint_list *CL;
	    FILE *fp;
	    Alignment *IN;
	    int **pos;

	    if (!D1)return 1;
	    IN=D1->A;
	    CL=(D1->A)->CL;
	    pos=aln2pos_simple(IN, IN->nseq);
	    fp=vfopen (out_file, "w");
	    fp=save_list_header (fp,CL);


	    for ( b=0; b< IN->nseq-1; b++)
	      {
		for ( c=b+1; c< IN->nseq; c++)
		  {
		    s1=IN->order[b][0];
		    s2=IN->order[c][0];
		    fprintf ( fp, "#%d %d\n", s1+1, s2+1);
		    for ( a=0; a< IN->len_aln; a++)
		      {
			r1=pos[b][a];
			r2=pos[c][a];

			if ( s1==s2 && !CL->do_self)continue;

			if ( s1< s2)s=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
			else        s=(CL->evaluate_residue_pair)( CL, s2, r2, s1, r1);

			s=(s!=UNDEFINED)?s:0;
			if ( r1>0 && r2>0)
			  {
			    fprintf (fp, "\t%5d %5d %5d \n", r1, r2, s);
			  }
		      }
		  }
	      }
	     vfclose (save_list_footer (fp, CL));
	  }
	
	else if (strm (out_format, "score") || strm (out_format, "score_aln") || strm (out_format, "score_seq") || strm (out_format, "score_pw"))
	  {
	    int i, x, y;
	    double **dm;
	    FILE *fp;

	    fp=vfopen (out_file, "w");
	    fprintf ( fp,"#SCORE_FORMAT_01\n"); 
	    
	    if ( strm (out_format, "score_aln") || strm (out_format, "score"))
	      fprintf (fp, "%d [SCORE_ALN]\n", (D1->A)->score_aln);
	    if ( strm (out_format, "score_seq") || strm (out_format, "score"))
	      {
		for (x=0; x<(D1->A)->nseq; x++)
		  fprintf (stdout, ">%-20s %d [SCORE_SEQ]\n", (D1->A)->name[x], (D1->A)->score_seq[x]);
	      }
	    if ( strm (out_format, "score_pw") || strm (out_format, "score"))
	      {
		if ((D1->A)->dm)dm=(D1->A)->dm;
		else if (D2 && D2->A && (D2->A)->dm)dm=(D2->A)->dm;
		else if (DST && DST->A && (DST->A)->dm)dm=(DST->A)->dm;
		if (dm)
		  {
		    for (x=0; x<(D1->A)->nseq; x++)
		      for (y=0; y<(D1->A)->nseq; y++)
			if (x!=y)fprintf (fp, ">%-20s %-20s %.3f [SCORE_PW]\n", (D1->A)->name[x], (D1->A)->name[y],(float)dm[x][y]);
		  }
	      }
	    vfclose (fp);
	  }
	else if      ( strncmp (out_format, "score",5)==0 || strm (out_format, "html"))
	  {
	    int a,l,n;
	    Alignment *BUF;
	    
	    if (!D1)return 1;

	    
	    if ( !DST)
	      {
		fprintf ( stderr,"\n[You Need an evaluation File: Change the output format or use +evaluate][FATAL:%s]\n", PROGRAM);
		myexit(EXIT_FAILURE);
	      }
	    if ( !strm ("html", out_format))while ( out_format[0]!='_' && out_format[0]!='\0' )out_format++;
	    
	    //l=(DST->A)->len_aln;
	    //n=(DST->A)->nseq;
	    //for (a=0; a<l; a++)HERE ("%s %s %d", out_format,(DST->A)->name[n],(DST->A)->seq_al[(DST->A)->nseq][a]);
	    
	    D1->S=aln2seq(D1->A);
	    BUF=copy_aln (DST->A, NULL);
	    
	   
	    
	    DST->A=aln2number (DST->A);
	    
	    
	    if     ( strstr ( out_format, "html" ))output_reliability_html  ( D1->A,  DST->A, out_file);
	    else if( strm ( out_format, "_ps"    ))output_reliability_ps    ( D1->A,  DST->A, out_file);
	    else if( strm ( out_format, "_pdf"   ))output_reliability_pdf   ( D1->A,  DST->A, out_file);
	    else if( strm ( out_format, "_ascii" ))output_reliability_ascii ( D1->A,  DST->A, out_file);
	    else if( strm ( out_format, "_fasta" ))output_reliability_fasta ( D1->A,  DST->A, out_file);
	    else if( strm ( out_format, "_raw"   ))output_raw_score ( D1->A,  DST->A, out_file);
	    else
	      {
		DST->A=BUF;
		main_output (DST, NULL, NULL, out_format+1, out_file);
	      }
	  }
	else if (strm(out_format, "compressed_ps"))
	  {
	    if (!D1) return 1;
	    aln2compressed_ps (D1->A, out_file);
	  }
	else if (strm(out_format, "compressed_pdf"))
	  {
	    if (!D1) return 1;
	    aln2compressed_pdf (D1->A, out_file);
	  }
	else if (strm (out_format, "sec_html") || strm (out_format, "_E_html"))
	  {
	    Alignment *ST, *A;
	    Sequence *S;

	    int a, b,c,i, ns=0;
	    char *buf;
	    if (!D1)return 1;
	    A=D1->A;


	    S=A->S;
	    ST=copy_aln (A, NULL);
	    for (a=0; a<ST->nseq; a++)
	      {
		i=name_is_in_list (ST->name[a],S->name, S->nseq, 100);
		if ( i!=-1)
		  {
		    buf=seq2E_template_string(S, i);
		    if ( buf==NULL)continue;
		    else ns++;
		    for (c=0,b=0; b<ST->len_aln; b++)
		      {
			int r1, s;
			r1=ST->seq_al[a][b];
			if ( r1!='-')
			  {
			    s=tolower (buf[c]);
			    if (s=='e')r1='0';
			    else if (s=='h')r1='9';
			    else if (s=='c')r1='5';
			    c++;
			  }
			ST->seq_al[a][b]=r1;
		      }
		  }
	      }

	    if (!ns)
	      {
		add_warning ( stderr, "Cannot output sec_html:_E_ template file (sec. struc.) is required for this output", PROGRAM);
	      }
	    output_color_html  ( A, ST, out_file);
	  }
	else if (strm (out_format, "tm_html") || strm (out_format, "_T_html"))
	  {
	    Alignment *ST, *A;
	    Sequence *S;

	    int a, b,c,i, ns=0;
	    char *buf;
	    if (!D1)return 1;
	    A=D1->A;
	    A->output_tm = 1;

	    S=A->S;
	    ST=copy_aln (A, NULL);
	    for (a=0; a<ST->nseq; a++)
	      {
		i=name_is_in_list (ST->name[a],S->name, S->nseq, 100);
		if ( i!=-1)
		  {
		    buf=seq2T_template_string(S, i);
		    if ( buf==NULL)continue;
		    else ns++;
		    for (c=0,b=0; b<ST->len_aln; b++)
		      {
			int r1, s;
			r1=ST->seq_al[a][b];
			if ( r1!='-')
			  {
			    s=tolower (buf[c]);
			    if (s=='o')r1='0';
			    else if (s=='h')r1='9';
			    else if (s=='i')r1='5';
			    c++;
			  }
			ST->seq_al[a][b]=r1;
		      }
		  }
	      }

	    if (!ns)
	      {
		add_warning ( stderr, "Cannot output tm_html:_T_ template file (trans. Memb. ) is required for this output", PROGRAM);
	      }
	    output_color_html  ( A, ST, out_file);
	  }

	else if (strm (out_format, "color_exoset"))
	  {
	    Alignment *ST, *EX, *A;
	    Constraint_list *CL;
	    int a, b, n;
	    char *buf;

	    if ( !DST->A)
	      {
		printf_exit ( EXIT_FAILURE, stderr, "\nYou must provide an obj file via the -struc_in flag [FATAL:%s]", PROGRAM);
	      }
	    EX=DST->A;
	    A=D1->A;

	    CL=declare_constraint_list ( DST->S,NULL, NULL, 0,NULL, read_matrice("pam250mt"));

	    ST=copy_aln (A, NULL);
	    buf=(char*)vcalloc ( EX->len_aln+1, sizeof (int));

	    for ( a=0; a< A->nseq; a++)
	      {
		int i;

		i=name_is_in_list (A->name[a],EX->name, EX->nseq, -1);
		if ( i==-1)continue;

		sprintf ( buf, "%s", EX->seq_al[i]);
		ungap (buf);

		for (n=0,b=0; b<A->len_aln; b++)
		  {
		    if (!is_gap(A->seq_al[a][b]))
		      {
			if ( buf[n]=='o')
			  ST->seq_al[a][b]='0';
			else if ( buf[n]=='j')
			  ST->seq_al[a][b]='1';
			else if ( buf[n]=='b')
			  ST->seq_al[a][b]='2';
			n++;
		      }
		  }
	      }
	    vfree (buf);

	    output_color_html  ( A, ST, out_file);
	    return EXIT_SUCCESS;
	  }

	else if (strm (out_format, "color_protogene"))
	  {
	    int n, a, b;
	    DST->A=copy_aln (D1->A, NULL);
	    for (n=1,a=0; a< (D1->A)->len_aln; a++, n++)
	      {
		for ( b=0; b<(D1->A)->nseq; b++)
		  {
		    if (is_gap((D1->A)->seq_al[b][a]));
		    else if ( n<=3)(DST->A)->seq_al[b][a]=2;
		    else if ( n>3)(DST->A)->seq_al[b][a]=9;
		  }

		if ( n==6)n=0;
	      }
	    output_color_html  ( D1->A,  DST->A, out_file);
	    return EXIT_SUCCESS;

	  }

	else if      ( strncmp (out_format, "color",5)==0)
	 {
	   Alignment *BUF;

	   if (!D1)return 1;

	   if ( !DST)
	     {
	       fprintf ( stderr,"\n[You Need an evaluation File: Change the output format or use +evaluate][FATAL:%s]\n", PROGRAM);
	       myexit(EXIT_FAILURE);
	     }
	   while ( out_format[0]!='_' && out_format[0]!='\0' )out_format++;

	   BUF=copy_aln (DST->A, NULL);




	   if     ( strm ( out_format, "_html"  ))output_color_html  ( D1->A,  DST->A, out_file);
	   else if( strm ( out_format, "_ps"    ))output_color_ps    ( D1->A,  DST->A, out_file);
	   else if( strm ( out_format, "_pdf"   ))output_color_pdf   ( D1->A,  DST->A, out_file);
	   else if( strm ( out_format, "_ascii"   ))output_color_ascii   ( D1->A,  DST->A, out_file);
	   else
	     {
	       DST->A=BUF;
	       return main_output (DST, NULL, NULL, out_format+1, out_file);
	     }
	   return EXIT_SUCCESS;
	 }
	else if ( strm4  ( out_format, "tc_aln","t_coffee_aln", "t_coffee", "tcoffee"))
	  {
	    if (!D1)return 1;
	    vfclose (output_aln ( D1->A, vfopen (out_file, "w")));
	  }
	else if ( strm  ( out_format, "analyse_pdb"))
	  {
	    if (!D1)return 1;
	    if ( !DST)
	      {
		fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL:%s]\n", PROGRAM);
		myexit(EXIT_FAILURE);
	      }
	    analyse_pdb ( D1->A,DST->A, "stdout");
	    (DST->A)=aln2number (DST->A);
	    output_reliability_ps    ( D1->A,  DST->A, out_file);
	  }
	else if ( strm4 ( out_format, "lower0", "lower1", "lower2", "lower3") || strm4(out_format, "lower4", "lower5", "lower6", "lower7") || strm4 (out_format,"lower8", "lower9", "align_pdb", "malign_pdb") )
	  {
	    if (!D1)return 1;
	    if ( !DST)
	      {
		fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL:%s]\n", PROGRAM);
		myexit(EXIT_FAILURE);
	      }



	    (DST->A)=aln2number (DST->A);
	    if ( strm (out_format, "align_pdb"))value=0;
	    else if (  strm (out_format, "malign_pdb"))value=5;
	    else value=atoi(out_format+5);

	    D1->A=filter_aln_upper_lower (D1->A, DST->A,0, value);
	    output_clustal_aln ( out_file, D1->A);
	  }
	else if ( strnm (out_format, "repeat", 6))
	  {
	    int size;
	    int a, b, c;
	    Alignment *CONC;

	    if ( !D1)return 1;
	    size=atoi (out_format+6);
	    print_aln (D1->A);
	    CONC=declare_aln2 ( (D1->A)->nseq, ((D1->A)->len_aln+1)*size+1);

	    for ( a=0; a< (D1->A)->nseq; a++)(D1->A)->seq_al[a][(D1->A)->len_aln]='\0';
	    for ( c=0,a=0; a< (D1->A)->nseq;c++)
	      {

		sprintf ( CONC->name[c], "%s", (D1->A)->name[a]);
		for ( b=0; b<size; b++, a++)
		  {
		    strcat (CONC->seq_al[c], (D1->A)->seq_al[a]);
		    strcat (CONC->seq_al[c], "O");
		  }
	      }
	    CONC->nseq=c;CONC->len_aln=strlen (CONC->seq_al[0]);
	    output_clustal_aln ( out_file, CONC);
	    free_aln (CONC);
	  }

	else if ( strnm (out_format, "upper", 5))
	  {
	    
	    if (!D1)return 1;
	    if ( !DST)
	      {
		fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL:%s]\n", PROGRAM);
		myexit(EXIT_FAILURE);
	      }
	    
	    
	    (DST->A)=aln2number (DST->A);
	    
	    value=atoi(out_format+5);
	    
	    D1->A=filter_aln_lower_upper (D1->A, DST->A,0, value);
	    output_clustal_aln ( out_file, D1->A);
	  }

	else if ( strstr (out_format, "tcs"))
	  {
	    int target;
	    int action;
	    char *p;
	    int value, c, s, filter;
	    Alignment *A, *EA;

	    if (!D1)return 1;
	    if ( !DST)
	      {
		fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL:%s]\n", PROGRAM);
		myexit(EXIT_FAILURE);
	      }
	    
	    if (D2 && D2->S)
	      {
		//must be DNA file
		D1->A =thread_dnaseq_on_prot_aln (D2->S, D1->A);
		DST->A=aln2X(DST->A,3);
	      }
	    
	    EA=(DST->A)=aln2number (DST->A);
	    A=D1->A;
	    
	   
	    
	    if (strstr (out_format, "residue"))target=0;
	    else target=1;

	    p=NULL;
	    if      (p=strstr (out_format, "filter")){action=0;p+=strlen ("filter");}
	    else if (p=strstr (out_format, "lower")){action=1;p+=strlen("lower");}
	    else if (p=strstr (out_format, "weighted")){action=2; p=NULL;}
	    else if (p=strstr (out_format, "replicate")){action=3; p+=strlen ("replicate");}
	    else action=1;
	    while (p && p[0] && p[0]=='_')p++;
	    if (p && p[0])sscanf (p, "%d", &filter);
	    else filter=3;
	    
	    if (action==2 || action ==3)
	      {
		Alignment *B;
		if ( check_file_exists (out_file))remove (out_file);
		
		
		B=copy_aln (A, NULL);
		B=realloc_aln (B, B->len_aln*10+1);
		B->len_aln=0;
		for (c=0; c<A->len_aln; c++) 
		  {
		    int vc=EA->seq_al[A->nseq][c];
		    int n;
		    
		    for (n=0; n<=vc; n++, B->len_aln++)
		      {
			for (s=0; s<A->nseq; s++)
			  B->seq_al[s][B->len_aln]=A->seq_al[s][c];
		      }
		  }
		if (action==2)
		  if (strstr (out_format, "fasta"))output_mfasta_aln ( out_file,B);
		  else if (strstr (out_format, "rphylip"))output_rphylip_aln ( out_file, B);
		  else output_phylip_aln ( out_file,B,"w");
		else if (action==3)
		  {
		    Alignment *C;
		    int a, b, pos, c;
		    
		    
		    C=copy_aln (A, NULL);
		    if (filter==3)filter=100;
		    for (a=0; a<filter; a++)
		      {
			for (b=0; b<A->len_aln; b++)
			  {
			    pos=(int)rand()%B->len_aln;
			    for (c=0; c<A->nseq; c++)
			      C->seq_al[c][b]=B->seq_al[c][pos];
			  }
			
			if (strstr (out_format, "fasta"))output_mfasta_aln ( out_file,C);
			else if (strstr (out_format, "rphylip"))output_rphylip_aln ( out_file, C);
			else output_phylip_aln ( out_file,C, "w");
		      }
		    free_aln(C);
		  }
		free_aln (B);
	      }
	    else
	      {
		    
		for (c=0; c<A->len_aln; c++)
		  for (s=0;s<A->nseq; s++)
		    {
		      int vr=EA->seq_al[s][c];
		      int vc=EA->seq_al[A->nseq][c];
		      int res=A->seq_al[s][c];
		      
		      if (is_gap(res))continue;
		      else if (target==0)
			{
			  if      (action==0 && vr<filter)A->seq_al[s][c]='-';
			  else if (action==1)
			    {
			      A->seq_al[s][c]=toupper (A->seq_al[s][c]);
			      if ( vr<=filter)A->seq_al[s][c]=tolower (A->seq_al[s][c]);
			    }
			}
		      else if (target ==1)
			{
			  if      (action==0 && vc<filter)A->seq_al[s][c]='-';
			  else if (action==1)
			    {
			      A->seq_al[s][c]=toupper (A->seq_al[s][c]);
			      if ( vc<=filter)A->seq_al[s][c]=tolower (A->seq_al[s][c]);
			    }
			}
		    }
		if (strstr (out_format, "fasta"))output_mfasta_aln ( out_file,D1->A);
		else if (strstr (out_format, "rphylip"))output_rphylip_aln ( out_file, D1->A);
		else output_phylip_aln ( out_file, D1->A, "w");
	      }
	    
	  }
	else if ( strm4 ( out_format, "filter0", "filter1", "filter2", "filter3"))
	  {
	    if (!D1)return 1;
	    if ( !DST)
	      {
		fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL:%s]\n", PROGRAM);
		myexit(EXIT_FAILURE);
	      }
	    
	    (DST->A)=aln2number (DST->A);
	    
	    D1->A=filter_aln (D1->A, DST->A, atoi(out_format+6));
	    output_clustal_aln ( out_file, D1->A);
	  }
	else if ( strm3 ( out_format, "rphylip_aln", "rphylip", "rphy"))
	  {
	    if (!D1)return 1;
	    if ( check_file_exists (out_file))remove (out_file);
	    output_rphylip_aln ( out_file, D1->A);
	  }
	else if ( strm3 ( out_format, "phylip_aln", "phylip", "phy"))
	  {
	    if (!D1)return 1;
	    if ( check_file_exists (out_file))remove (out_file);
	    output_phylip_aln ( out_file, D1->A, "w");
	  }
	else if ( strm ( out_format, "mocca_aln"))
	  {
	    if (!D1)return 1;
	    output_mocca_aln ( out_file, D1->A, DST->A);
	  }
	else if ( strm ( out_format, "saga_pw_sd_weights") )
	  {
	    if (!D1)return 1;
	    output_pw_weights4saga ((D1->W),(D1->W)->PW_SD, out_file);
	  }
	else if ( strm ( out_format, "saga_aln"))
	  {
		if (!D1)return 1;
		output_saga_aln (out_file, D1->A);
		}
	else if (strm2 ( out_format, "aln","clustal_tc")|| strm (out_format, "msa"))
	  {

	    if (!D1)return 1;
	    output_clustal_aln (out_file, D1->A);
	  }
	else if (strm5 ( out_format, "strict_clustal","clustal_aln", "clustalw","clustal", "clustalw_aln") || strm (out_format,"number_aln"))
	  {
	    if (!D1)return 1;
	    output_strict_clustal_aln (out_file, D1->A);
	  }
	else if ( strm ( out_format, "conc_aln"))
	      {
		if (!D1)return 1;
		output_conc_aln (out_file, D1->A);
		}
	else if ( strm2 ( out_format, "lalign_aln","lalign"))
	        {
		if (!D1)return 1;
		output_lalign (out_file, D1->A);
		}
	else if ( strm2 ( out_format, "glalign_aln","glalign"))
	        {
		if (!D1)return 1;
		output_glalign (out_file, D1->A, DST->A);
		}

	else if ( strm3 ( out_format, "fasta_aln","fasta","raw_fasta" ) || strm (out_format, "blast_aln"))
		{
		if (!D1)return 1;
		output_fasta_aln( out_file, D1->A);
		}
		else if ( strm(out_format, "xmfa"))
		{
			if (!D1)return 1;
			output_xmfa_aln( out_file, D1->A);
		}

	else if ( strstr (out_format, "overaln"))
		{

		  char *s, mode[100];
		  OveralnP *F;
		  int eb=0;
		  if (!D1) return 1;
		  F=(OveralnP*)vcalloc (1, sizeof (OveralnP));
		  ungap_aln (D1->A);
		  string_array_upper ((D1->A)->seq_al, (D1->A)->nseq);
		  if ( D2 && D2->A)
		    {
		      D1->A=mark_exon_boundaries (D1->A, D2->A);
		      eb=1;
		    }
		  else if ( (s=get_string_variable ("exon_boundaries")))
		    {
		      Sequence *S;
		      Alignment *EB;
		      EB=seq2aln(S=main_read_seq(s),NULL, 0);
		      D1->A=mark_exon_boundaries (D1->A, EB);
		      free_sequence (S, S->nseq); free_aln (EB);
		      eb=1;
		    }
		  if ( strstr (out_format, "lower")) sprintf (F->mode,"lower");
		  else if (strstr (out_format, "unalign2"))sprintf (F->mode, "unalign2");
		  else if (strstr (out_format, "unalign"))sprintf (F->mode, "unalign");
		  else sprintf (F->mode, "%s", ((s=get_string_variable ("overaln_mode")))?s:"lower");
		  if (!strm (F->mode, "lower") && !strm (F->mode, "unalign") && !strm (F->mode, "unalign2"))printf_exit (EXIT_FAILURE,stderr,"\nERROR: unknown overaln_mode in overaln output [%s] [FATAL:%s]", mode, PROGRAM);

		  if (int_variable_isset ("overaln_threshold"))F->t=get_int_variable ("overaln_threshold");
		  if (int_variable_isset ("overaln_target"))F->f=get_int_variable ("overaln_target");
		  if (int_variable_isset ("overaln_P1"))F->p1=get_int_variable ("overaln_P1");
		  if (int_variable_isset ("overaln_P2"))F->p2=get_int_variable ("overaln_P2");
		  if (int_variable_isset ("overaln_P3"))F->p3=get_int_variable ("overaln_P3");
		  if (int_variable_isset ("overaln_P4"))F->p4=get_int_variable ("overaln_P4");

		  if (eb)sprintf (F->model, "fsa2");
		  else   sprintf (F->model, "fsa1");
		  D1->A=aln2clean_pw_aln (D1->A, F);

		  //if (eb)D1->A=aln2clean_pw_aln (D1->A, mode,t, f,p1,p2,p3, "fsa2");
		  //else   D1->A=aln2clean_pw_aln (D1->A, mode,t, f,p1,p2,p3, "fsa1");

		  D1->S=aln2seq(D1->A);
		  output_clustal_aln (out_file, D1->A);
		}
	else if ( strm ( out_format, "est_prf" ))
		{
		if (!D1)return 1;
		output_est_prf( out_file, D1->A);
		}
	else if ( strm ( out_format, "clean_est_fasta_seq" ))
		{
		if (!D1)return 1;
		D1->A=clean_est(D1->A);
		output_fasta_seq(out_file, D1->A);

		}

	else if ( strm3 ( out_format, "msf_aln", "gcg", "msf"))
		{
		if (!D1)return 1;
		output_msf_aln( out_file, D1->A);
		}
	else if ( strm ( out_format, "rnalign"))
		{
		if (!D1)return 1;
		output_rnalign (out_file, D1->A, DST->S);
		}
	else if ( strm ( out_format, "tblastx_db1"))
	  {
	    seq2tblastx_db (out_file,D1->S,1);
	  }
	else if ( strm ( out_format, "tblastx_db") || strm (out_format, "tblastx_db3"))
	  {
	    seq2tblastx_db (out_file,D1->S,3);
	  }
	else if ( strm ( out_format, "tblastx_db2"))
	  {
	    seq2tblastx_db (out_file,D1->S,2);
	  }
	else if ( strm ( out_format, "fasta_seq") ||strm ( out_format, "list")||strm ( out_format, "file_list"))
	  {
	    
		if (!D1)return 1;
		output_fasta_seq (out_file,D1->A);
	  }
	else if (strm (out_format, "treelist") ||strm (out_format, "nexus") )
	  {
	    output_treelist (out_file, D1->A);
	  }
	else if (strm (out_format, "fasta_tree") )
	  {
	    if (!D1)return 1;
	    output_fasta_tree (out_file,D1->A);
	  }
	
	else if ( strm ( out_format, "gotoh_seq"))
	  {
	    if (!D1)return 1;
	    output_gotoh_seq (out_file,D1->A);
	  }
	else if ( strm (out_format, "fasta_seq1"))
	        {
		if (!D1)return 1;
		output_fasta_seq1 (out_file, D1->A);
		}
	else if ( strm (out_format, "fasta_seq2"))
	  {
	    if (!D1)return 1;
	    output_fasta_seq2 (out_file, D1->A);
	  }
	else if ( strm2 (out_format, "pir_aln", "pir"))
		{
		if (!D1)return 1;
		output_pir_aln (out_file, D1->A);
		}
	else if ( strm (out_format, "pir_seq"))
		{
		if (!D1)return 1;
		output_pir_seq (out_file, D1->A);
		}
        else if ( strm (out_format, "gor_seq"))
		{
                if (!D1)return 1;
		output_gor_seq (out_file, D1->A);
		}
	else if ( strm (out_format, "pir_seq1"))
		{
		  if (!D1)return 1;
		output_pir_seq1 (out_file, D1->A);
		}
	else if ( strm (out_format, "pw_lib_saga_aln"))
		{
		  if (!D1)return 1;
		output_pw_lib_saga_aln (out_file, D1->A);
		}
	else if ( strm (out_format, "lib"))
		{
		  if (!D1)return 1;
		output_lib (out_file, D1->A);
		}
	else if ( strm (out_format, "pdb_constraint_list"))
	        {
		  if (!D1)return 1;
		output_constraints (out_file, "pdb",D1->A);
		}
	else if ( strm  (out_format, "vienna2tc_lib"))
	  {
	    D1->CL=vienna2tc_lib (out_file, D1->S,D2->S);
	    save_contact_constraint_list (D1->CL, out_file);
	  }
	else if ( strm  (out_format, "vienna2template"))
	  {
	    vienna2template_file (out_file, D1->S,D2->S);
	    
	  }
	else if ( strm3 (out_format, "contact_list","contact_lib","c_lib"))
	  {
	    save_contact_constraint_list (D1->CL, out_file);
	  }
	else if ( strm2 (out_format, "constraint_list","tc_lib"))
	  {
	    
	    if (!D1)return 1;
	    else if (!D1->CL)output_constraints (out_file,"sim", D1->A);
	    else if (D1->CL) vfclose ( save_constraint_list ( D1->CL, 0, (D1->CL)->ne, out_file, NULL, "ascii",(D1->CL)->S));
	  }
	else if (  strm2 (out_format, "extended_lib","extended_cosmetic"))
	  {
		  if (!D1)return 1;
		  output_constraints (out_file,out_format, D1->A);
	  }
	else if ( strncmp (out_format, "extended_pair", 13)==0)
	  {
	    if (!D1)return 1;
	    output_constraints (out_file,out_format, D1->A);
	  }
	else if ( strm (out_format, "cache_id"))
	  {
	    if (!D1)return 1;
		  cache_id (D1->A);
		  output_saga_aln (out_file, D1->A);
	  }
        else if ( strm (out_format, "compress_aln"))
	  {
	    if (!D1)return 1;
	    compress_aln (D1->A);
	    output_saga_aln (out_file, D1->A);
	  }
	
	else if (strm (out_format, "n_seq") ||strm (out_format, "nseq") )
		{
		  if (!D1)return 1;
		fp=vfopen ( out_file, "w");
		fprintf ( fp, "%d\n", (D1->A)->nseq);
                vfclose (fp);
		}

	
	else if ( strm ( out_format, "tdna_fasta_seq1"))
	        {if (!D1)return 1;
		D1->A=translate_dna_aln (D1->A,0);
		output_fasta_seq1 (out_file, D1->A);
		}
	else if (strm (out_format, "exons"))
	  {
	    Alignment *A;
	    //exons come in upper case
	    //output alternates them upper/lower
	    if (!D1)return 1;
	    A=copy_aln (D1->A, NULL);
	    A->seq_al=gene2exons(A->seq_al,A->nseq);
	    output_fasta_seq (out_file,A);
	    free_aln (A);
	  }
	else if ( strm (out_format, "wexons"))
	  {
	    if (!D1)return 1;
	    output_wexons (out_file,D1->A);

	  }
	else if ( strm (out_format, "texons"))
	  {
	    Alignment *A;
	    Sequence *S;
	    //exons come in upper case
	    //output alternate amino acids in upper/lower case
	    //amino acid has the case of its first nucleotide
	    if (!D1)return 1;
	    A=copy_aln (D1->A, NULL);
	    A->seq_al=gene2exons(A->seq_al,A->nseq);
	    S=aln2seq(A);
	    output_fasta_seqS (out_file,S=translate_dna_seqS(S,1,'X'));
	  }
	else if ( strm (out_format, "sexons"))
	  {
	    Alignment *A;

	    //exons come in upper case
	    //output alternate amino acids in upper/lower case
	    //amino acid has the case of its first nucleotide
	    if (!D1)return 1;
	    A=copy_aln (D1->A, NULL);
	    output_fasta_seq ( out_file, D1->A);
	  }

	else if ( strm ( out_format, "tdna_aln"))
	        {if (!D1)return 1;
		D1->A=translate_dna_aln (D1->A,0);
		output_saga_aln ( out_file, D1->A);
		}
	else if ( strm ( out_format, "cdna_fasta_seq1"))
	        {if (!D1)return 1;
		D1->A= gene2prot(D1->A);
		output_fasta_seq1 ( out_file, D1->A);
		}
	else if ( strm ( out_format, "mutate_cdna_aln"))
	        {if (!D1)return 1;
		    D1->A= mutate_cdna_aln ( D1->A);
		    output_clustal_aln ( out_file, D1->A);
		}
	else if ( strm ( out_format, "tdna_sp_aln"))
	        { if (!D1)return 1;
	        if ( !DST)
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL:%s]\n", PROGRAM);
		   myexit(EXIT_FAILURE);
		   }
	       (DST->A)=aln2number (DST->A);
		D1->A=translate_splice_dna_aln (D1->A, DST->A);
		output_saga_aln ( out_file, D1->A);
		}
	else if (out_format && out_format[0] && (strcmp ( out_format,"rna_graph_fasta")==0))
		{
		  if (!D1)return 1;
		sprintf ( (D1->A)->seq_al[0], "%s",(DST->S)->seq[0]);
		(D1->A)->nseq=0;
		output_fasta_seq (out_file, DST->A);
		}
	else if (strm ( out_format, "freq_mat"))
	        {
		  if (!D1)return 1;
		output_freq_mat (out_file, D1->A);
		}
	else if (strm ( out_format, "maln_pval"))
	        {if (!D1)return 1;
		output_maln_pval ( out_file, D1->A);
		}
	else if ( strm ( out_format, "model_aln"))
	        {
		  if (!D1)return 1;
		output_model_aln ( out_file, D1->A);
		}
	else if (strncmp (out_format, "mult",4)==0)
	        {
		  if (!D1)return 1;
		output_mult_fasta_seq ( out_file, D1->A, atoi(out_format+4));
		}
	else if (strm (out_format, "conservation"))
	  {
	    output_conservation_statistics (out_file, D1->A);
	  }
	else if (strm (out_format, "len"))
	        {
		  if (!D1)return 1;
		  output_statistics (out_file, D1->A, "nrl");
		}
	else if ( strm (out_format, "name"))
	        {
		  if (!D1)return 1;
		  if ( D1->A)output_statistics (out_file, D1->A, "n");
		  if ( D1->T)
		    {
		      Sequence *TS;
		      TS=tree2seq(D1->T, NULL);print_array_char (vfopen(out_file, "w"), TS->name, TS->nseq, "\n");
		    }
		}
	else if ( strm (out_format, "code_name"))
	        {
		  char **nl=NULL;
		  int num, n=0;
		  Sequence *TS;
		  FILE *lfp;
		  if ( D1->A){n=(D1->A)->nseq, nl=(D1->A)->name;}
		  if ( D1->T){TS=tree2seq(D1->T, NULL);nl=TS->name;n=TS->nseq;}

		  lfp=vfopen (out_file, "w");
		  for ( num=0; num<n; num++)
		    fprintf (lfp, "\n%s C%d", nl[num], num+1);
		  fprintf (lfp, "\n");
		  vfclose (lfp);
		}
	else if ( strm ( out_format, "seq2struc"))
		  {
		    output_seq2struc (out_file, D1->A);
		  }
	else if ( strstr  ( out_format, "pavie_age_channel"))
	  {
	    output_n_pavie_age_channel ( D1->S,out_file, atoi((out_format+strlen ("pavie_age_channel"))));
	    return EXIT_SUCCESS;
	  }
	else if ( strstr ( out_format, "age_matrix"))
		  {
		    output_age_matrix (out_file, atoi((out_format+10)));
		  }
	else if ( strm ( out_format, "transitions"))
		  {
		    output_transitions (out_file, D1->A);
		  }

	else if ( strncmp (out_format, "statistics",10)==0)
	        {
		  if (!D1)return 1;

		  output_statistics (out_file, D1->A,out_format+10);
		}

	
	else if ( strm  (out_format, "newick_shuffle"))
	  {
	    no_rec_print_tree_shuffle (D1->T, stdout);
	    fprintf (stdout, ";\n");
	  }
	else if ( strm  (out_format, "newick_randomize"))
	  {
	    no_rec_print_tree_randomize (D1->T, stdout);
	    fprintf (stdout, ";\n");
	  }
	else if ( strm  (out_format, "dm"))
	  {
	    int x;
	    FILE *fpx;
	    if (D1 && D1->A && (D1->A)->Tree && ((D1->A)->Tree)->dmF_list)
	      {
		Alignment *T=(D1->A)->Tree;
		FILE *fpx=vfopen (out_file, "w");
		display_file_content (fpx,T->dmF_list[0]);
		if (int_variable_isset ("print_replicates")) 
		  for (x=1; x<T->nseq; x++)
		    {
		      display_file_content (fpx,T->dmF_list[x]);
		    }
		vfclose (fpx);
	      }
	    else
	      {
		printf_exit ( EXIT_FAILURE, stderr, "Distance Matrix Not produced\n", PROGRAM);
	      }
	  }
	else if ( strm  (out_format, "newick_dm"))
	  {
	    int x;
	    FILE *fpx;
	    if (D1 && D1->A && (D1->A)->Tree && ((D1->A)->Tree)->dmF_list)
	      {
		Alignment *T=(D1->A)->Tree;
		FILE *fpx=vfopen (out_file, "w");
		fprintf (fpx, "%s\n", T->seq_al[0]);
		display_file_content (fpx,T->dmF_list[0]);
		if (int_variable_isset ("print_replicates")) 
		  for (x=1; x<T->nseq; x++)
		    {
		      fprintf (fpx, "%s\n", T->seq_al[0]);
		      display_file_content (fpx,T->dmF_list[x]);
		    }
		vfclose (fpx);
	      }
	    else
	      {
		printf_exit ( EXIT_FAILURE, stderr, "Distance Matrix Not produced\n", PROGRAM);
	      }
	  }
	else if ( strm  (out_format, "mafftnewick"))
	  {
	    if (D2 && D2->S)
	      vfclose (tree2file (D1->T, D2->S, "mafftnewick", vfopen (out_file, "w")));
	    else if (D1 && D1->S)
	      vfclose (tree2file (D1->T, D1->S, "mafftnewick", vfopen (out_file, "w")));
	    else
	      printf_exit (EXIT_FAILURE,stderr,"ERROR -output=mafftnewick requires -in2=<seqfile> [FATAL]");
	  }
	else if ( strm  (out_format, "mafftdndmatrix"))
	  {
	    FILE *lfp1;
	    FILE *lfp2;
	    char *mafftdnd=vtmpnam (NULL);
	    char *mafftmatrix=vtmpnam (NULL);
	    char cc;
	    reset_tree_distances (D1->T, 1);
	    if (D2 && D2->S)
	      vfclose (tree2file (D1->T, D2->S, "mafftnewick", vfopen (mafftdnd, "w")));
	    else if ( D1 && D1->S)
	      vfclose (tree2file (D1->T, D1->S, "mafftnewick", vfopen (mafftdnd, "w")));
	    else
	      printf_exit (EXIT_FAILURE,stderr,"ERROR -output=mafftdndmat requires -in2=<seqfile> [FATAL]");
	    
	   
	    printf_system ("newick2mafft.rb %s > %s 2>/dev/null", mafftdnd,mafftmatrix);
	    

	    
	    lfp1=vfopen (mafftmatrix, "r");
	    lfp2=vfopen (out_file, "w");
	    while ((cc=fgetc (lfp1))!=EOF)fprintf (lfp2, "%c", cc);
	    vfclose (lfp1);vfclose (lfp2);
	  }
	else if ( strm4 (out_format, "newick_tree","newick","binary","nh"))
	        {
		  if (!D1)return 1;
		  if (!D1->T && (D1->A))
		    {		      
		      Alignment *T;
		      FILE *fpx;
		      int x;
		      
		      if ((D1->A)->Tree)
			{
			  T=(D1->A)->Tree;
			}
		      else T=(D1->A);
		      fpx=vfopen (out_file, "w");
		      fprintf (fpx, "%s\n", T->seq_al[0]);
		      if (int_variable_isset ("print_replicates"))
			{
			  for (x=1; x<T->nseq; x++)
			    fprintf (fpx, "%s\n", T->seq_al[x]);
			}
		      vfclose (fpx);
		    }
		  else
		    {
		      print_newick_tree (D1->T, out_file);
		    }
		}
	
	    
	else if ( strncmp (out_format, "sarsim", 6)==0)
	        {
		  if (!D1)return 1;
		  compare_sar_sequence (D1->S, (D2 &&D2->S)?D2->S:D1->S, atoi(out_format+6));
		  return EXIT_SUCCESS;
		}

	else if ( strncmp (out_format, "sim",3)==0)
	        {
		  if (!D1)return 1;
		  output_similarities (out_file, D1->A,out_format);
		}

	else if ( strncmp (out_format, "cov",3)==0)
	        {
		  if (!D1)return 1;
		  output_similarities (out_file, D1->A,out_format);
		}
	else if ( strm (out_format, "stockholm_aln") ||strm (out_format, "stockholm") )
	  {
	    output_stockholm_aln (out_file,D1->A, (D2)?D2->A:NULL);
	  }
	else if ( strm (out_format, "pair_sim"))
	  {
	    if ( !D2)
	      {
		fprintf ( stderr, "\n-output=pair_sim: provide aln1 via -in and aln2 via -in2 [FATAL:%s]\n", PROGRAM);
		myexit (EXIT_FAILURE);
	      }
	    output_similarities_pw (out_file, D1->A,D2->A,out_format);
	  }
	else if ( strm (out_format, "matrix") || strm (out_format, "blast_matrix"))
	  {
	    output_blast_mat (D1->M, out_file);
	  }
	else if ( strm (out_format, "header_matrix"))
	  {
	    output_header_mat( D1->M, out_file);
	  }

	else
	        {

		    fprintf ( stderr, "\n%s is an UNKNOWN OUTPUT FORMAT [FATAL:%s]\n",out_format, PROGRAM);
		    myexit (EXIT_FAILURE);

		}

	//Remove the expansion
	if ( expanded)
	  {
	    free_aln (D1->A);
	    D1->A=BUF_A;
	  }
	return 0;
	}
int is_in_format_list ( char *name)
	{
	if ( strcmp ( name, "saga_aln")==0)return 1;
	if ( strcmp ( name, "number_aln")==0)return 1;
	if ( strcmp ( name, "clustal_aln")==0)return 1;
	if ( strcmp ( name, "fasta_aln")==0)return 1;
	if ( strcmp ( name, "number_fasta")==0)return 1;
	if ( strcmp ( name, "fasta_seq")==0)return 1;
	if ( strcmp ( name, "pdb")==0)return 1;
	if ( strcmp ( name, "msf_aln")==0)return 1;
	if ( strcmp ( name, "dali_aln")==0)return 1;
	if ( strcmp ( name, "dali_seq")==0)return 1;
	if ( strcmp ( name, "barton_list_tc")==0)return 1;
	if ( strcmp ( name, "est_prf")==0)return 1;

	if ( strcmp ( name, "gotoh_aln")==0)return 1;
	if ( strcmp ( name, "amps_aln")==0)return 1;
	if ( strcmp ( name, "pir_aln")==0)return 1;
	if ( strcmp ( name, "pir_seq")==0)return 1;
	if ( strcmp ( name, "est_fasta")==0)return 1;
	if ( strcmp ( name, "amps_sd_scores")==0)return 1;
	if ( strcmp ( name, "pima_aln")==0)return 1;
	if ( strcmp ( name, "dialign_aln")==0)return 1;
	if ( strcmp ( name, "gor_seq")==0)return 1;
	if ( strcmp ( name, "gor_struc")==0)return 1;
	if ( strcmp ( name, "stockholm_aln")==0)return 1;

	return 0;
	}
int is_struc_in_format_list ( char *name)
	{
	if ( strcmp ( name, "rna_number")==0)return 1;
	if ( strcmp ( name, "fasta_seq")==0)return 1;
	return 0;
	}
char *format_name2aln_format_name (char *name)
        {
	  if ( strm (name, "gcg"))sprintf (name, "msf");
	  else if ( strm (name, "fasta"))sprintf (name, "fasta_aln");
	  return name;
	}
int is_out_format_list ( char *name)
	{
	  return main_output (NULL, NULL, NULL, name, NULL);
	}

int is_struc_out_format_list ( char *name)
	{
	  return main_output (NULL, NULL, NULL, name, NULL);
	}

/**************************************************************************************************/
/*************************************REFORMAT UTIL*************************************************/
/**************************************************************************************************/

/*************************************REFORMAT IN**************************************************/
/**************************************************************************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               READ COG FILE                                             */
/*                                                                                         */
/***************************************************************************************** */

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT WEIGHTS                                            */
/*                                                                                         */
/***************************************************************************************** */

Weights* get_amps_sd_scores ( char *fname)
	{
	FILE *fp;
	char *buf;
	char *buf2;
	int nseq;
	Weights *W;
	int a, b,e;
	int c;
	float array[20];

	buf=(char*)vcalloc ( 1001, sizeof (char));
	buf2=(char*)vcalloc ( 1001, sizeof (char));

	fp=vfopen ( fname, "r");
	set_fp_id ( fp, "Index");
	buf=fgets ( buf, 1000, fp);
	fscanf ( fp, "%s", buf2);

	nseq=0;
	while ( isalnum(buf2[0]) && !isalpha(buf2[0]))
		{
		nseq++;
		buf=fgets ( buf, 1000, fp);
		fscanf ( fp, "%s", buf2);
		}
	vfclose ( fp);

	W=declare_weights (nseq);

	fp=vfopen ( fname, "r");
	set_fp_id ( fp, "Index");
	buf=fgets ( buf, 1000, fp);
	fscanf ( fp, "%s", buf2);

	a=0;
	while ( isalnum(buf2[0]) && !isalpha(buf2[0]))
		{
		fp=set_fp_after_char (fp, '>');
		fscanf ( fp, "%s",W->seq_name[a]);
		buf=fgets ( buf, 1000, fp);
		fscanf ( fp, "%s", buf2);
		a++;
		}
	buf=fgets ( buf, 1000, fp);
	c=1;
	while ( c!=0)
		{
		for ( e=0; e< 16; e++)
			{
			c=fscanf ( fp, "%f", &array[e]);
			}
		fscanf ( fp, "\n");
		if ( c!=0)
			{

			a=(int)array[0]-1;
			b=(int)array[1]-1;
			W->PW_ID[b][a]=W->PW_ID[a][b]=array[9];
			W->PW_SD[b][a]=W->PW_SD[a][b]=array[14];
			}

		}
	vfclose ( fp);
	sprintf ( W->comments, "SD WEIGHTS GENERATED WITH THE PROGRAM AMPS IN PAIRWISE MODE");
	vfree ( buf);
	return W;
	}

Weights *read_seq_weight (char **name, int nseq, char* seq_weight)
       {
       int a, p;
       Weights *W;
       float w;

       FILE *fp;
       char line[LONG_STRING];
       char sname[MAXNAMES];


       /*Read sequence weights:
	* comment
	name1 weight1
	.....


	NOTE:
	weights must be between 0 and 1;

	sequences not in S do not get any weight
	sequences in S but not in file get a weight of 1
       */
       if ( !is_single_seq_weight_file (seq_weight))
	 {
	   fprintf ( stderr, "\nERROR: File %s is not in Format SINGLE_SEQ_WEIGHT_FORMAT_01 [FATA:%s]", seq_weight,PROGRAM);
	   myexit (EXIT_FAILURE);
	   return NULL;
	 }
       else
	 {
	   W=declare_weights(nseq);
	   for ( a=0; a< nseq; a++)
	     {
	       sprintf ( W->seq_name[a], "%s", name[a]);
	       W->SEQ_W[a]=1;
	     }
	   sprintf ( W->mode, "%s", seq_weight);
	   fp=vfopen (seq_weight, "r");


	   while ( fgets( line,LONG_STRING-1, fp))
	     {
	       if ( line[0]=='*' ||line[0]=='#' || isblanc(line));
	       else
		 {
		   if (sscanf(line, "%s %f", sname, &w)!=2)continue;
		   if ( (p=name_is_in_list ( sname, W->seq_name, nseq, MAXNAMES-1))!=-1)
		     {
		       W->SEQ_W[p]=w;
		     }
		 }
	     }
	   vfclose (fp);
	   return W;
	 }
       }


/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT MISC                                               */
/*                                                                                         */
/***************************************************************************************** */

char *** read_rename_file ( char *fname, int code)
{
  int n;
  FILE *fp;
  char ***convert=NULL;

  convert=(char***)declare_arrayN(3, sizeof (char),count_n_line_in_file(fname) +1,2,MAXNAMES+1);
  fp=vfopen (fname, "r");
  n=0;
  if ( code==CODE)      while ( fscanf ( fp, "%s %s\n", convert[n][0], convert[n][1])==2)n++;
  else if (code==DECODE)while ( fscanf ( fp, "%s %s\n", convert[n][1], convert[n][0])==2)n++;
  vfclose (fp);
  return convert;
}

void get_barton_list_tc_seq ( char *in_file)
	{
	FILE *fp, *fp_make, *fp_length, *fp_long;
	FILE *fp_small[9];

	static char *buf;
	int len_buf=10000;
	char name[100];

	char pwd[100];
	int a,c,nseq;
	int k=0;
	int *length;
	int longest=0;

	c=0;
	length=(int*)vcalloc ( 1000, sizeof(int));
	if ( buf==NULL)buf=(char*)vcalloc ( len_buf, sizeof (char));
	fp=vfopen (in_file, "r");
	fp_long=vfopen ( "barton_seq_list_large", "w");
	fp_make=vfopen ( "make_dir", "w");
	fp_length=vfopen ( "barton_length", "w");
	for ( a=0; a< 9; a++)
		{
		sprintf ( name, "barton_nseq%d",a);
		fp_small[a]=vfopen ( name, "w");
		}
	get_pwd (pwd);


	while ( c!=EOF)
		{a=0;
		while ( (c=fgetc(fp))!='#');
		while ( (c=fgetc(fp))=='#');
		ungetc ( c, fp);
		while ( (c=fgetc(fp))!='#')buf[a++]=c;
		buf[a]='\0';

		sprintf ( name, "%s", buf);

		while ( (c=fgetc(fp))=='#');

		if ( c!=EOF)
			{
			a=0;
			while ( (c=fgetc(fp))!='#' && c!=EOF)
				{
				buf[a++]=c;
				if (a==len_buf)
					{
					len_buf+=10000;
					buf=(char*)vrealloc ( buf, len_buf*sizeof (char));
					}
				}
			buf[a]='\0';
			if (c!=EOF)
				{

				nseq=process_barton_entry ( buf,name);
				length[nseq]++;
				longest=(longest<nseq)?nseq:longest;

				if ( nseq<=8) fprintf ( fp_small[nseq], "%s.pep\n", name);
				else fprintf ( fp_long, "%s.pep\n",name);
				fprintf ( fp_make, "mkdir %s\nmv %s.pep %s\nmv %s.check %s\n", name, name, name, name, name);
				k++;
				}


			}
		}

	vfclose (fp);
	vfclose (fp_long);
	for ( a=0; a< 9; a++)vfclose (fp_small[a]);
	vfclose (fp_make);
	for ( a=0; a<= longest; a++)fprintf ( fp_length, "%d: %d\n", a, length[a]);
	vfclose ( fp_length);

	}

int process_barton_entry (char *buf, char *name)
    {
    Alignment *LA;
    Sequence *LS;
    int a,c;
    static char *buf2;
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0;
    int l;
    char fname[100];
    char com_name[100];
    int rm_gap=1;

    sprintf ( fname, "%s.pep", name);
    sprintf ( com_name, "%s.check",name);

    if ( buf2==NULL)buf2=(char*)vcalloc ( 10000, sizeof (char));
    a=0;
    while (buf[a]!='\0')
	 	{
		 if (buf[a]=='>')
			{
			a=get_string_line (a,2, buf, buf2);
			while ((c=buf[a++])!='*')
				if (isalnum (c)|| c=='.' || c=='-')
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			 clen=0;
			}
		if ( buf[a]!='\0')a++;
		}


    LS=declare_sequence (  min_len_seq,  max_len_seq,  nseq);
    LS->nseq=nseq;


    for (a=0, current=0; current< nseq; current++)
    	{
    	a=get_string_line ( a, 1, buf, buf2);
    	sscanf ( buf2, ">P1;%s", LS->name[current]);
    	a=get_string_line ( a, 1, buf, buf2);
    	l=strlen ( buf2);
	buf2[l-1]='\0';
    	sprintf ( LS->seq_comment[current],"%s", buf2);

    	p=0;
    	while ( (c=buf[a++])!='*')
    		{
    		if (isalpha (c))
			LS->seq[current][p++]=tolower (c);
		else if ( isgraph(c))
			LS->seq[current][p++]=(c);
		}
    	a++;
    	}

    LA=declare_Alignment(LS);
    seq2aln ( LS, LA,rm_gap);
    output_fasta_seq (fname,LA);
    output_pir_check (com_name,LA->nseq, LA->seq_comment);
    free_Alignment ( LA);
    free_sequence ( LS, nseq);

    return nseq;
    }




Structure *read_rna_struc_number (Alignment *A,char *fname)
	{
	FILE *fp;
	int a;
	char x,y;
	float f;
	Sequence *SA;
	Structure *ST;
	int first, last;

	SA=declare_sequence ( A->len_aln, A->len_aln, 1);
	SA->len[0]=A->len[0];
	for ( a=0; a< SA->len[0]; a++)
		SA->seq[0][a]='.';
	ST=declare_rna_structure_num (SA);
	ST->S=SA;

	fp=vfopen ( fname, "r");
	fscanf ( fp, "%c\n%d\n",&x, &(ST)->tot_list);
	for ( a=0; a<(ST)->tot_list; a++)
		{
		fscanf ( fp, "%d %d %d %c %c %f\n", &(ST)->list[a][0],&(ST)->list[a][1],&(ST)->list[a][2], &x, &y, &f);
		(ST)->list[a][0]--;
		(ST)->list[a][1]--;
		(ST)->list[a][2]--;
		if ( a==0)
			{
			(ST)->stem[0][0]=(ST)->list[a][0];
			(ST)->stem[0][1]=a;
			}
		else if ( (ST)->stem[(ST)->tot_stem][0]==(ST)->list[a][0]);
		else if ( (ST)->stem[(ST)->tot_stem][0]!=(ST)->list[a][0])
			{
			(ST)->stem[(ST)->tot_stem][2]=a-1;
			(ST)->tot_stem++;
			(ST)->stem[(ST)->tot_stem][0]=(ST)->list[a][0];
			(ST)->stem[(ST)->tot_stem][1]=a;
			}

		SA->seq[0][(ST)->list[a][1]]='-';
		SA->seq[0][(ST)->list[a][2]]='-';
		}
	(ST)->stem[(ST)->tot_stem][2]=a-1;
	(ST)->tot_stem++;
	for ( a=0; a< (ST)->tot_stem; a++)
     		{

     		first=(ST)->stem[a][1];
     		last=(ST)->stem[a][2];
     		SA->seq[0][(ST)->list[first][1]]='>';
     		SA->seq[0][(ST)->list[first][2]]='<';
     		SA->seq[0][(ST)->list[last][1]]='>';
     		SA->seq[0][(ST)->list[last][2]]='<';
     		}

	return ST;
	}

Structure * declare_rna_structure_num (Sequence *SA)
	{
	Structure *ST;
	ST=( Structure*)vcalloc ( 1, sizeof ( Structure));
	ST->list=declare_int ( SA->len[0], 3);
	ST->stem=declare_int ( SA->len[0], 3);
	return ST;
	}
char ** read_lib_list (char *name, int *n)
{

  char **lines;
  char **list;
  int a, b, l;

  lines=file2lines (name);
  l=atoi (lines[0]);

  list=(char**)vcalloc (l, sizeof (char*));
  for ( n[0]=0,a=1; a<l; a++,b++)
    if ( !strstr (lines[a], "TC_LIB_LIST_FORMAT_01"))list[n[0]++]=lines[a];
  vfree (lines);
  return list;
}

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT SEQ                                                */
/*                                                                                         */
/***************************************************************************************** */
char ***read_group ( char *file)
{
  /*Format: Fasta like, the name fo the group followed with the name of the sequences
    ><Group name> <First Seq> <second seq> ....
    Groups must NOT be overlaping
    list[group_index][0]="number of sequences"
    list[group_index][1]="group name"
    list[group_index][2...N]="sequence"
  */

  FILE *fp;
  char *buf;
  char ***list;
  int a, c, l;



  l=measure_longest_line_in_file (file)+1;
  buf=(char*)vcalloc (l, sizeof (char));
  list=(char***)vcalloc ( count_n_line_in_file (file )+1, sizeof (char**));

  fp=vfopen (file, "r");

  a=0;
  while ((c=fgetc(fp))!=EOF)
    {
      buf=fgets (buf,l-1, fp);
      if ( c=='>')list[a++]=string2list (buf);
    }
  vfclose (fp);
  vfree (buf);
  return list;
}

int *pdb2atom_pos_list(char *pdb)
{
  char **ll=file2lines (pdb);
  int n, nn,li, ci;
  int *pos;
  
  nn=1;
  while (ll[nn])nn++;
  pos=(int*)vcalloc (nn+2, sizeof (int));
  
  nn=1;n=0;li=-1;
  while (ll[nn])
    {
      char **lw=string2list (ll[nn]);
      if (strm (lw[1], "ATOM"))
	{
	  ci=atoi (lw[6]);
	  if (ci!=li)
	    {
	      pos[n++]=ci;
	    }
	  li=ci;
	}	
      free_char (lw, -1);
      nn++;
    }
  pos[n]=-1;
  free_char (ll, -2);
  return pos;
}  

int seqres_equal_atom (char*fname)
{
  Sequence *S1=get_pdb_sequence_from_field(fname, "SEQRES");
  Sequence *S2=get_pdb_sequence_from_field(fname, "ATOM");
  int value;
  
  if (!S1||!S2)value=0;
  else if (!strm(S1->seq[0], S2->seq[0]))value=0;
  else value=1;
  
  free_sequence (S1, NULL);
  free_sequence (S2, NULL);
  return value;
}
  

Sequence* get_pdb_sequence   (char *fname)
{
  Sequence *S;

  if ( (S=get_pdb_sequence_from_field(fname, "SEQRES"))!=NULL);
  else if ( (S=get_pdb_sequence_from_field(fname, "ATOM"))!=NULL)
    {
      add_warning (stderr,"Read Sequence from ATOM field in %s [%s:WARNING]", fname, PROGRAM);
    }
  else
    {
      add_warning ( stderr, "failed to extract sequence from %s [%s:WARNING]", fname, PROGRAM);
      S=NULL;
    }
  return S;
}


Sequence* get_pdb_sequence_from_field   (char *fname, char *field)
{
	char *tp_name;
	char *command;
	char *pdbid;
	Sequence *S;


	command=(char*)vcalloc ( LONG_STRING, sizeof (char));
	tp_name=vtmpnam (NULL);
	sprintf ( command, "extract_from_pdb -seq_field %s -chain FIRST -infile \'%s\' -mode fasta > %s", field, check_file_exists(fname), tp_name);
// 	printf("CO: %s\n", command);
// 	char *x = vcalloc ( LONG_STRING, sizeof (char));
// 	sprintf(x, "cp %s ~/Desktop/erg.txt", tp_name);
// 	my_system(x);
	if ( getenv4debug ("DEBUG_EXTRACT_FROM_PDB"))fprintf ( stderr, "\n[DEBUG_EXTRACT_FROM_PDB:get_pdb_seq] %s\n", command);
	my_system ( command);


	S=get_fasta_sequence ( tp_name, NULL);
	if (S==NULL)return NULL;

	if ( (pdbid=get_pdb_id (fname))){sprintf ( S->name[0], "%s",pdbid);vfree (pdbid);}
	S->nseq=1;

	sprintf ( S->file[0], "%s", fname);
	S->max_len=S->min_len=S->len[0];
	if ( S->len[0]==0)
	{
		free_sequence (S, -1);
		S=NULL;
	}

	vremove ( tp_name);
	vfree ( command);

	return S;
}

char * get_pdb_file   ( char *fname)
     {
	 char *file;
	 int a, c;
	 FILE *fp;


	 a=0;
	 file=(char*)vcalloc ( sizeof (char),count_n_char_in_file ( fname)+1);
	 fp=vfopen ( fname, "r");
	 while ( (c=fgetc(fp))!=EOF)file[a++]=c;
	 file[a]='\0';
	 return file;
     }

Sequence* get_struc_gor ( char *fname)
    {
    int nseq, min_len, max_len;
    int a, c;
    int len;
    char name[STRING];


    FILE *fp;
    Sequence *S;

    min_len=max_len=-1;
    fp=vfopen ( fname, "r");
    nseq=0;
    while ( (c=fgetc(fp))!=EOF)
	    {
	    if ( c!='!');
	    else
		{
		nseq++;
		fscanf ( fp, "%s %d", name, &len);
		if (min_len==-1)min_len=max_len=len;
		else
		    {
		    min_len=(len>min_len)?min_len:len;
		    max_len=(len>max_len)?len:max_len;
		    }
		}

	    }
    vfclose (fp);

    S=declare_sequence (  min_len,  max_len+1,nseq);
    S->nseq=0;

    fp=vfopen (fname,"r");
     while ( (c=fgetc(fp))!=EOF)
	     {
	     if ( c!='!');
	     else
	        {
		fscanf ( fp, "%s %d\n",S->name[S->nseq], &(S->len[S->nseq]));

		while ( (c=fgetc(fp))!='\n');

		for ( a=0; a<S->len[S->nseq]; a++)
		    fscanf ( fp, " %*c %c %*f %*f %*f\n",&(S->seq[S->nseq][a]));

		S->seq[S->nseq][a]='\0';
		while ( (c=fgetc(fp))!='!' && c!=EOF);
		ungetc (c, fp);
		S->nseq++;
		}

	     }
    vfclose (fp);
    return S;
    }

Sequence* get_sequence_dali (char *fname)
    {
    Sequence *LS;
    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0;

    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  myexit(EXIT_FAILURE);
	 }
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (isdigit(c))
			{
			ungetc(c, fp);
			fscanf (fp, "%s",name);
			while (!isdigit(c=fgetc(fp)) && c!=EOF)
				if (isalnum (c) || c=='.' || c=='-')
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);

    LS=declare_sequence (  min_len_seq,  max_len_seq+1,nseq);
    LS->nseq=nseq;

    fp=vfopen (fname,"r");

    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (isdigit(c))
			{
			ungetc(c, fp);
			fscanf_seq_name (fp, LS->name[current]);
			p=0;
			while (!isdigit(c=fgetc(fp)) && c!=EOF)
				{
				if (isalpha (c))
				    LS->seq[current][p++]=tolower (c);
				else if ( c=='.')
				    LS->seq[current][p++]='-';
				else if ( c=='-')
				    LS->seq[current][p++]='-';
				}
			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);


    return LS;
    }

Sequence* get_dialign_sequence (char *fname)
    {
    Sequence *LS;
    FILE *fp;
    int c;

    char name[10000];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0, l=0;
    char *buf;

    buf=(char*)vcalloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  myexit(EXIT_FAILURE);
	 }
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{fscanf (fp, "%s",name);

			buf=fgets ( buf, 1000, fp);
			while ((c=fgetc(fp))!='>' && c!=EOF && c!=' ' && c!='\t')
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);

    LS=declare_sequence (  min_len_seq,  max_len_seq, nseq);
    LS->nseq=nseq;

    fp=vfopen (fname,"r");

    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='>')
			{

			fscanf_seq_name (fp, LS->name[current]);
			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==',')LS->name[current][l-1]='\0';
			buf=fgets ( buf, 1000, fp);
			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF && c!=EOF && c!=' ' && c!='\t')
				if (isalpha (c))
				    LS->seq[current][p++]=tolower (c);
				else if ( isgraph(c))
				    LS->seq[current][p++]=(c);
			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);
    return LS;
    }

Sequence* get_pima_sequence (char *fname)
    {
    Sequence *LS;

    FILE *fp;
    int c;

    char name[10000];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0, l=0, len=0;
    char *buf, *buf2;
    char prefix[1000];

    sprintf (  prefix, "%s",fname);

    buf=strstr(prefix, "-");
    buf[0]='\0';
    len=strlen (prefix);



    buf=(char*)vcalloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  myexit(EXIT_FAILURE);
	 }
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{
			  fscanf_seq_name (fp,name);
			  if ( strlen(name)>=len && strncmp ( name, prefix, len)==0)
				{
				  c=fgetc(fp);
				}
			  else
				{

				buf=fgets ( buf, 1000, fp);
				while ((c=fgetc(fp))!='>' && c!=EOF)
					if (isalnum (c)|| is_gap(c))
						clen++;
				 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 	 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 	nseq++;
				clen=0;
				}
			}
		else
		    	c=fgetc (fp);
		}
    vfclose (fp);

    LS=declare_sequence (  min_len_seq,  max_len_seq, nseq);
    LS->nseq=nseq;

    fp=vfopen (fname,"r");

    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='>')
			{
			fscanf_seq_name (fp,LS->name[current]);
			if ( strlen(LS->name[current])>=len && strncmp ( LS->name[current], prefix, len)==0)
				c=fgetc (fp);
			else
				{
				buf2=strstr (LS->name[current], ".");
				if ( buf2!=NULL) buf2[0]='\0';

				l=strlen ( LS->name[current]);
				if ( LS->name[current][l-1]==','||LS->name[current][l-1]==',')LS->name[current][l-1]='\0';
				buf=fgets ( buf, 1000, fp);
				p=0;
				while ((c=fgetc(fp))!='>' && c!=EOF)
					if (isalpha (c))
					    LS->seq[current][p++]=tolower (c);
					else if ( isgraph(c))
					    LS->seq[current][p++]=(c);
				LS->seq[current][p]='\0';
				LS->len[current]=strlen ( LS->seq[current]);
				current++;
				}
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);
    return LS;
    }

Sequence* perl_reformat2fasta (char *perl_command, char *fname)
{
	char command[1000];
	char *file;

	file=vtmpnam (NULL);

	check_program_is_installed ( perl_command,"", perl_command,EMAIL,IS_FATAL);
	sprintf ( command, "%s %s > %s", perl_command, fname, file);
	my_system ( command);
	return get_fasta_sequence (file, NULL);
}


// int fscanf_seq_name ( FILE *fp, char *sname)
// {
// 	static char *name;
// 	int r;
// 	if ( !sname) return 0;
//
// 	if ( !name)name=vcalloc ( 10000, sizeof (char));
// 	fscanf (fp, "%s", name);
// 	r=strlen (name);
// 	if ( r>MAXNAMES)
// 		add_warning (stderr, "Seq Name Too long: [%s]. Truncated to %d", name, MAXNAMES);
// 	name[MAXNAMES]='\0';
// 	sprintf ( sname, "%s", name);
// 	return r;
// }


void check_seq_name (char *sname)
{

	char *tmp = strchr(sname, ' ');
	if (tmp != NULL)
		*tmp='\0';
	int x = strlen (sname);
	if ( x>MAXNAMES)
		add_warning (stderr, "Seq Name Too long: [%s]. Truncated to %d", sname, MAXNAMES);
// 	if (sname[x-1] == '\n')
// 		sname[MAXNAMES]='\0';
}




Sequence *get_file_list ( char *fname)
{

  char ***list;
  char *tmp;
  int a;
  FILE *fp;

  tmp=vtmpnam (NULL);
  list=file2list (fname, "\n");
  fp=vfopen (tmp, "w");
  a=0;
  while (list[a] && !isspace(list[a][1][0]))
    {

      fprintf ( fp, ">%s\n", list[a][1]);
      a++;
    }
  vfclose (fp);
  free_arrayN((void ***)list, 3);
  return get_fasta_sequence (tmp, NULL);
}

Sequence *get_fasta_sequence_num (char *file, char *comment_out)  
{
  return get_fasta_sequence(file, comment_out);
}
Sequence *get_fasta_sequence_raw (char *file, char *comment_out)  
{
  return get_fasta_sequence(file, comment_out);
}
Sequence *reload_seq(Sequence *A)
{
  int a;
  FILE *fp;
  char *tmp=vtmpnam (NULL);
  dump_seq(A, tmp);
  free_sequence (A,-1);
  return get_fasta_sequence(tmp, NULL);
}
Sequence *get_fasta_sequence (char *file, char *comment_out)
{
  long *map;
  char *s;
  int i, nseq,a;
  Sequence *A;
  int maxlen=0;

  
  map=fasta2map(file);
  
  nseq=read_array_size_new (map)-1;
  
  A=declare_sequence(0,0,nseq);
  A->nseq=nseq;
  
  A->nseq=nseq;
  for(i=0;i<nseq; i++)
    {
      s=file2record_it(file,i, map);
      A->seq[i]=csprintf (A->seq[i], "%s", seq2clean(FastaRecord2seq(s)));
      A->name[i]=csprintf (A->name[i], "%s", FastaRecord2name(s));
      A->seq_comment[i]=csprintf (A->seq_comment[i], "%s", FastaRecord2comment(s));
      A->len[i]=strlen (A->seq[i]);
      A->max_len=( A->max_len>A->len[i])?A->len[i]:A->max_len;
    }
  file2record_it(NULL,0, NULL);
  vfree (map);
  
  return A;
}


Sequence* get_pir_sequence (char *fname, char *comment_out)
    {
    Sequence *LS;

    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0, l=0;
    char *buf;

    buf=(char*)vcalloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  myexit(EXIT_FAILURE);
	 }
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{
			if ( (c=fgetc(fp))=='P')while ( (c=fgetc(fp))!=';');
			else ungetc ( c, fp);
			fscanf_seq_name (fp,name);

			buf=fgets ( buf, 1000, fp);
			while ((c=fgetc(fp))!='>' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);



    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq);
    LS->nseq=nseq;

    fp=vfopen (fname,"r");

    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='>')
			{
			if ( (c=fgetc(fp))=='P')while ( (c=fgetc(fp))!=';');
			else ungetc ( c, fp);

			fscanf_seq_name (fp,LS->name[current]);

			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==',')LS->name[current][l-1]='\0';
			LS->name[current]=translate_name ( LS->name[current]);
			buf=fgets ( buf, 1000, fp);

			LS->seq_comment[current]=fgets ( LS->seq_comment[current],COMMENT_SIZE-1, fp);
			LS->seq_comment[current][strlen(LS->seq_comment[current])-1]='\0';
			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF)
				if (isalpha (c))
				    LS->seq[current][p++]=tolower (c);
				else if ( !isspace(c) && c!='*')
				    LS->seq[current][p++]=(c);
			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);
    if (comment_out!=NULL) output_pir_check ( comment_out,LS->nseq, LS->seq_comment);
    return LS;
    }

Sequence* get_gor_sequence (char *fname, char *comment_out)
    {
    Sequence *LS;

    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=99999;
    int nseq=0;
    char *buf;

    buf=(char*)vcalloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  myexit(EXIT_FAILURE);
	 }
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='!')
			{
			fscanf_seq_name (fp,name);

			buf=fgets ( buf, 1000, fp);
			while ((c=fgetc(fp))!='!' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);

    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq);
    LS->nseq=nseq;

    fp=vfopen (fname,"r");

    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='!')
			{


			fscanf_seq_name (fp,LS->name[current]);
			LS->name[current]=translate_name ( LS->name[current]);
			buf=fgets ( buf, 1000, fp);

			p=0;
			while ((c=fgetc(fp))!='!' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
				    LS->seq[current][p++]=tolower (c);

			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);

    return LS;
    }



int fscanf_seq_name ( FILE *fp, char *sname)
{
	static char *name;
	int r;
	if ( !sname) return 0;

	if ( !name)name=(char*)vcalloc ( 10000, sizeof (char));
	fscanf (fp, "%s", name);
	r=strlen (name);
	if ( r>MAXNAMES)
		add_warning (stderr, "Seq Name Too long: [%s]. Truncated to %d", name, MAXNAMES);
	name[MAXNAMES]='\0';
	sprintf ( sname, "%s", name);
	return r;
}


/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT ALN                                                 */
/*                                                                                         */
/***************************************************************************************** */
Alignment* undump_msa ( Alignment *A, char *tmp)
{
  return quick_read_fasta_aln (A, tmp);
}
char* dump_msa (Alignment *A, char *tmp)
{
  FILE*fp;
  int a;

  if (!A)return NULL;
  if (!tmp)tmp=vtmpnam(NULL);
  if (!(fp=vfopen (tmp, "w")))return NULL;
  
  for (a=0; a<A->nseq; a++)
    fprintf ( fp, ">%s %s\n%s\n", A->name[a], A->aln_comment[a], A->seq_al[a]);
  vfclose (fp);
  return tmp;
}
char* dump_seq (Sequence *A, char *tmp)
{
  FILE*fp;
  int a;

  if (!A)return NULL;
  if (!tmp)tmp=vtmpnam(NULL);
  if (!(fp=vfopen (tmp, "w")))return NULL;
  
  for (a=0; a<A->nseq; a++)
    fprintf ( fp, ">%s %s\n%s\n", A->name[a], A->seq_comment[a], A->seq[a]);
  vfclose (fp);
  return tmp;
}    

void undump_msa_old ( Alignment *A, char *tmp)
{
  FILE *fp;
  int m;
  char *buf;
  int index;

  if ( !A || !tmp || !check_file_exists (tmp))return;
  m=measure_longest_line_in_file (tmp );
  A=realloc_aln2 ( A,A->max_n_seq,m+1);

  buf=(char*)vcalloc (m+1, sizeof (char));
  fp=vfopen (tmp, "r");
  while (fscanf (fp, "%d %s\n", &index, buf)==2)
    {
      sprintf ( A->seq_al[index], "%s", buf);
    }
  vfclose (fp);
  vfree (buf);
}
void dump_msa ( char *file,Alignment *A, int nseq, int *lseq)
{
  FILE *fp;
  int a;
  fp=vfopen (file, "w");
  for (a=0; a<nseq; a++)
    fprintf ( fp, "%d %s\n", lseq[a], A->seq_al[lseq[a]]);
  vfclose (fp);
}

void read_stockholm_aln (char *file_name, Alignment *A)
{
  char *tmp_name;
  Sequence *S;

  
  tmp_name=vtmpnam (NULL);
  if (printf_system ( "clustalw_aln2fasta_aln.pl %s > %s",file_name, tmp_name)!=EXIT_SUCCESS)
    {
      printf_exit ( EXIT_FAILURE, stderr, "Could Not Read File %s [FATAL:%s]\n", file_name, PROGRAM);
    }
  else
    {
      int a;
      S=get_fasta_sequence ( tmp_name,NULL);
      for (a=0; a<S->nseq; a++)
	{
	  if (strstr (S->name[a], "_stockholm"))
	    {
	      substitute ( S->name[a], "_stockholmspace_", " ");
	      substitute ( S->name[a], "_stockholmhasch_", "#");
	    }
	}
      A=seq2aln (S, A, 0);
    }
  return;
}
Alignment* read_blast_aln ( char *file_name, Alignment *A)
{
  char *tmp_name;

  int type;


  if ( !(type=is_blast_file (file_name)))
    {
      myexit (EXIT_FAILURE);
    }
  tmp_name=vtmpnam ( NULL);
  if (type==BLAST_TXT)
    {
      printf_system("cat %s | blast_aln2fasta_aln.pl | fasta_aln2fasta_aln_unique_name.pl >%s", file_name, tmp_name);
    }
  else if (type==BLAST_XML)
    {

      printf_system("blast_xml2fasta_aln.pl %s >%s", file_name, tmp_name);
    }

  main_read_aln (tmp_name, A);
  return A;
}


void read_number_aln ( char *file_name, Alignment *A)
   {
    FILE *fp, *fp2;
    int * ptr_aln;
    int a,b,d;
    int c;
    char *buf=NULL;

    int tot=0;
    int flag=0;
    char *fname;
    int n_comment=0;

    int nseq=0;
    int max_len=0;


    fp=vfopen ( file_name, "r");

    fname=vtmpnam(NULL);
    fp2=vfopen ( fname, "w");
    while ( (c=fgetc(fp))!=EOF)
	{
	    fprintf ( fp2, "%c", c);
	}
    vfclose (fp);
    vfclose (fp2);


    /*1 Count The number of sequences*/
    fp=vfopen ( fname, "r");
    buf=vfgets ( buf,fp);
    if ( !isblanc (buf));
    while ( isblanc (buf))
    	{
	  buf=vfgets ( buf, fp);
    	}
    while (!isblanc (buf))
    	{
    	buf=vfgets ( buf,fp);
	}
    while ( !isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);
    	buf=vfgets ( buf,fp);
    	}

    if ( c!='\n')ungetc(c,fp);

    while ( isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);
    	a=0;
    	while ( isgraph ((c=fgetc(fp))));
        nseq++;
    	buf=vfgets ( buf, fp);
    	}
    vfclose (fp);

    /*DONE*/
    /*2 get_max_len*/
    max_len=count_n_char_in_file(fname)/nseq;
    A=realloc_alignment2( A, nseq+1, max_len+1);

    /*DONE*/


    fp=vfopen ( fname, "r");
    buf=vfgets ( buf, fp);
    if ( !isblanc (buf))sprintf (A->aln_comment[n_comment++], "%s", buf);
    while ( isblanc (buf))
    	{
        buf=vfgets ( buf,fp);
    	}
    while (!isblanc (buf))
    	{
    	buf=vfgets ( buf, fp);
	sprintf ( A->aln_comment[n_comment++], "%s", buf);

	}
    while ( !isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);
    	buf=vfgets ( buf, fp);

    	}

    if ( c!='\n')ungetc(c,fp);

    while ( isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);

	fscanf_seq_name (fp, A->name[A->nseq]);

	if ( name_is_in_list (A->name[A->nseq], A->name, A->nseq, 100)!=-1)
	  {
	    fprintf ( stderr, "(read_number_aln): Sequence %s Duplicated in File %s ", A->name[A->nseq], A->file[A->nseq]);
	    if (!getenv("ALLOW_DUPLICATE"))
	      {
		fprintf ( stderr, " [FATAL:%s]\n", PROGRAM);
		myexit (EXIT_FAILURE);
	      }
	  }
    	A->nseq++;
    	buf=vfgets ( buf,fp);
    	}

    vfclose (fp);



    if ((fp=vfopen ( fname, "r"))==NULL)
	printf ( "\nCOULDN'T READ %s", fname);

    ptr_aln=(int*)vcalloc ( A->nseq, sizeof(int));
    while ( flag==0)
	{
	while (  (c=fgetc(fp))!='\n');
	if ( (c=fgetc(fp))=='\n')
	    flag=1;
	}
    while ( !isalnum(c=fgetc(fp)));
    ungetc ( c, fp);
    while ( c!=EOF)
	{
	tot=0;
	while(tot< A->nseq && c!=EOF)
	    {
	     b=0;
	     while ( !isgraph (c=fgetc(fp)) && c!=EOF);
	     if ( c!=EOF)ungetc(c, fp);
	     while ( isgraph((buf[b++]=fgetc(fp))));
	     buf[b-1]='\0';
	     for ( a=-1,d=0; d< A->nseq; d++)
		if ( strcmp (A->name[d], buf)==0)
		    {a=d;
		     tot++;
		    }

	     if ( a==-1) while ( (c=fgetc(fp))!='\n' && c!=EOF);
	     else
	       {
		 while ( (c=fgetc(fp))!='\n')
		   {
		     if ( isgraph(c) || is_gap(c))
		       {if ( isalpha(c))
			 c=(A->residue_case==KEEP_CASE)?c:tolower(c);

		       if (!isspace(c))A->seq_al[a][ptr_aln[a]++]=c;
		       }
		   }
	       }
	     }
	 while ( !isalnum(c=getc(fp)) && c!=EOF);
	 if ( c!=EOF)
	    ungetc (c, fp);
	 }

    vfclose (fp);


    for ( a=0; a< A->nseq; a++)
	{A->seq_al[a][ptr_aln[a]]='\0';
	 A->order[a][0]=a;
	 A->order[a][1]=0;
	}

    A->len_aln= strlen(A->seq_al[0]);

    vfree (buf);
    vfree(ptr_aln);
    vremove (fname);

    }
void read_amps_aln ( char *in_file, Alignment *A)
	{
	FILE *fp;
	int a, b, c, cont=1;
	A->nseq=get_amps_seq_name ( A->name, in_file);

	fp=vfopen ( in_file, "r");
	fp=set_fp_id(fp, "1*");
	while ( (c=fgetc(fp))!='\n');
	b=0;
	while ( cont==1)
		{
		c=fgetc ( fp);
		c=fgetc(fp);
		if ( c=='*')
			{
			cont=0;
			for ( a=0; a<A->nseq; a++)
				A->seq_al[a][b]='\0';
			A->len_aln=b;
			}

		else
		    	{
		    	ungetc (c, fp);
		  	for ( a=0; a< A->nseq; a++)
		  		{
		  		c=fgetc(fp);
		  		if ( c==' ')A->seq_al[a][b]='-';
		  		else
		  			{
		  			A->seq_al[a][b]=c;
		  			A->len[a]++;
		  			}
		  		}
		  	while ((c=fgetc(fp))!='\n');
		  	b++;
		  	}
		}
	}






int get_amps_seq_name ( char **name, char* fname)
	{
	FILE *fp;
	int nseq=0;

	fp=vfopen ( fname, "r");
	fp=set_fp_id ( fp, "Index");
	while ( (fgetc(fp))!='\n');
	while ( isspace(fgetc(fp)))
		{fscanf (fp, "%*d >%s", name[nseq++]);
		 while ( (fgetc(fp))!='\n');
		}
	vfclose ( fp);
	return nseq;
	}
Alignment * read_gotoh_aln ( char *fname, Alignment *A)
   {
    FILE *fp;
    int * ptr_aln;
    int a,b,d,e;


    char *buf;
    char buf2[VERY_LONG_STRING+1];
    char buf3[VERY_LONG_STRING+1];
    char buf4[VERY_LONG_STRING+1];

    int tot=0;

    int l;
    int nseq, max_len;


    if ( !check_file_exists (fname))return NULL;
    fp=vfopen ( fname, "r");

/*1 GET THE NUMBER OF SEQUENCES*/
    nseq=0;
    buf=(char*)vcalloc ( VERY_LONG_STRING+1, sizeof (char));
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while (!isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( !isblanc ( buf) && buf!=NULL)
    	{
    	a=-1;
     	d=sscanf ( buf, "%d %s %s %s", &a, buf2, A->name[A->nseq],buf3);
    	if ( a!=-1)
    		{
		if ( name_is_in_list (A->name[A->nseq], A->name, A->nseq, 100)!=-1)
		  {
		    fprintf ( stderr, "\nWARNING (get_amps_seq_name): Sequence %s Duplicated in File %s ", A->name[A->nseq], A->file[A->nseq]);
		    if (!getenv("ALLOW_DUPLICATE"))
		      {
			fprintf ( stderr, " [FATAL:%s]\n", PROGRAM);
			myexit (EXIT_FAILURE);
		      }
		  }
    		nseq++;
    		fgets(buf, VERY_LONG_STRING, fp);
    		}
    	else ( buf=NULL);
    	}
    vfclose (fp);
/*2 Get the MAX Len and Reallocate*/
    max_len=count_n_char_in_file(fname)/nseq;
    A=realloc_aln2( A, nseq+1, max_len+1);
/*3 Get The Sequences Names*/
    A->nseq=0;
    fp=vfopen ( fname, "r");
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while (!isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( !isblanc ( buf) && buf!=NULL)
    	{
    	a=-1;
     	d=sscanf ( buf, "%d %s %s %s", &a, buf2, A->name[A->nseq],buf3);
    	if ( a!=-1)
    		{
    		if ( d==4)sprintf (A->name[A->nseq],"%s", buf3);
    		A->nseq++;
    		fgets(buf, VERY_LONG_STRING, fp);
    		}
    	else ( buf=NULL);
    	}
    vfclose (fp);

/*READ THE ALN*/
    fp=vfopen ( fname, "r");

    buf=(char*)vcalloc ( VERY_LONG_STRING+1, sizeof (char));;
    ptr_aln=(int*)vcalloc ( A->nseq, sizeof(int));

    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while (!isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));


    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));

    while (buf!=NULL)
	{
	tot=0;
	while(tot< A->nseq)
	    {

	    e=sscanf (buf, "%d %s %s %s", &e, buf2, buf3, buf4);
	    if ( e==4)sprintf( buf3, "%s", buf4);


	    for ( d=0; d< A->nseq; d++)
		{

		if ( strcmp (A->name[d], buf3)==0)
		    {a=d;
		     tot++;
		    }
		}
	     l=strlen (buf2);
	     if ( buf2[l-1]=='|')l--;
	     buf2[l]='\0';

	     for (b=0; b<l; b++)
	     	{
	     	if ( isgraph (buf2[b]))
	     	 	A->seq_al[a][ptr_aln[a]++]=(A->residue_case==KEEP_CASE)?buf2[b]:tolower (buf2[b]);
	     	 }
	     buf=fgets(buf, VERY_LONG_STRING, fp);
	     }
	 if ( buf!=NULL)
	 	{
	 	buf=fgets(buf, VERY_LONG_STRING, fp);
	 	while ( isblanc (buf) && buf!=NULL)
	 		{
	 		buf=fgets ( buf, VERY_LONG_STRING, fp);
	 		}
	 	}

	 }

    vfclose (fp);


    for ( a=0; a< A->nseq; a++)
	{A->seq_al[a][ptr_aln[a]]='\0';
	}

    A->len_aln= strlen(A->seq_al[0]);



    for ( a=0; a< A->nseq; a++)
    	{
    	for ( b=0; b< A->len_aln; b++)
    		A->len[a]+=1-is_gap(A->seq_al[a][b]);
    	}
    for ( a=0, b=0; a< A->len_aln; a++)
    	{
    	if ( !is_gap(A->seq_al[0][a]) &&!is_gap(A->seq_al[1][a]))b++;
    	}
    return A;
    }





void read_msf_aln ( char *fname, Alignment *A)
   {
    char command[1000];
    char *tmp_name;
    Sequence *S;

    tmp_name=vtmpnam(NULL);
    sprintf ( command, "msf_aln2fasta_aln.pl %s > %s", fname, tmp_name);

    if ( my_system (command)!=EXIT_SUCCESS)
      {
	fprintf ( stderr, "\nERROR: file %s does not have a legal msf format [FATAL:%s]", fname,PROGRAM);
	myexit (EXIT_FAILURE);
      }

    S=get_fasta_sequence ( tmp_name,NULL);
    A=seq2aln (S, A, 0);
    vremove (tmp_name);
    return;
    }

/**************************************************************************************************/
/*************************************REFORMAT OUT*************************************************/
/**************************************************************************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT MATRICES                                           */
/*                                                                                         */
/***************************************************************************************** */



int output_freq_mat ( char *outfile, Alignment *A)
    { /*
	function documentation: start

	int output_freq_mat ( char *outfile, Aligmnent *A)

	This function counts the number of residues in each column of an alignment (Prot)
	It outputs these values in the following format

	A | 0 0 0 1 0
	B | 1 0 0 0 1
	- | 0 1 1 0 0

	This format can be piped into:
	The routine used for computing the p-value  gmat-inf-gc-v2c

	function documentation: end
      */

    int a, b;
    int **freq_mat;
    FILE *fp;


    freq_mat=aln2count_mat (A);

    fp=vfopen ( outfile, "w");
    for ( b=0; b< 26; b++)
      {
	fprintf (fp, "%c |", 'A'+b);
	for ( a=0; a< A->len_aln; a++)fprintf (fp,"%d ", freq_mat[b][a]);
	fprintf (fp, "\n");
      }
    fprintf (fp, "- |");
    for ( a=0; a< A->len_aln; a++)fprintf (fp,"%d ", freq_mat[26][a]);

    free_int (freq_mat, -1);
    vfclose ( fp);
    return 1;
    }
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT P-Values                                           */
/*                                                                                         */
/***************************************************************************************** */
float output_maln_pval ( char *outfile, Alignment *A)
    {
      /*
	function documentation: start
	float output_maln_pval ( char *outfile, Aligmnent *A)

	This function outputs the p-value of a multiple alignmnet as described
	in Hertz, Stormo, Bioinformatics, 15-7/8, 563/577
	    ftp beagle.colorado.edu /pub/cosensus
	Locally
	    packages/consensus/gmat-inf-gc-v2c


	The routine used for computing the p-value is the program gmat-inf-gc-v2c
	function documentation: end
      */


    char *mat;
    char *result;
    FILE *fp;
    float value;
    char command[LONG_STRING];
    char string[STRING];
    mat=vtmpnam (NULL);
    result=vtmpnam (NULL);

    output_freq_mat (mat,A);
    sprintf ( command, "more %s | gmat-inf-gc-v2c -A abcdefghijklmnopqrstuvwxyz> %s",mat, result);
    my_system ( command);

    if ( !check_file_exists(result))return 0;
    fp=find_token_in_file ( result, NULL, "ln(p-value):");

    fscanf ( fp, "%s",string);
    value=atof ( string);
    vfclose ( fp);

    vremove ( mat);
    vremove ( result);

    fp=vfopen ( outfile, "w");
    fprintf ( fp, "%.6f\n", value);
    vfclose ( fp);

    return value;
    }


/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT WEIGHTS                                            */
/*                                                                                         */
/***************************************************************************************** */
int output_seq_weights ( Weights *W, char *wfile)
        {
	FILE*fp;
	int a;

	if ( W==NULL)return 0;

	fp=vfopen (wfile, "w");
	if ( fp==NULL)return 0;


	for ( a=0; a< W->nseq; a++)
		{

		  fprintf ( fp, "%s %.2f\n", W->seq_name[a],W->SEQ_W[a]);
		}
	vfclose ( fp);
	return 1;
	}
void output_pw_weights4saga ( Weights *W, float **w_list, char *wfile)
	{
	FILE*fp;
	int a, b;
	fp=vfopen (wfile, "w");

	fprintf ( fp, "%s\n$\n", W->comments);
	for ( a=0; a< W->nseq-1; a++)
		{
		for (b=a+1; b< W->nseq; b++)
			{
			fprintf ( fp, "%s %s %f\n", W->seq_name[a], W->seq_name[b],w_list[a][b]);
			}
		}
	fprintf ( fp, "$\n");
	vfclose ( fp);
	}

FILE * display_weights (Weights *W, FILE *fp)
{
  int a;
  int max_len;

  if ( W==NULL || strm (W->mode, "no_seq_weight"))
    {
      fprintf ( fp, "\n\nUN-WEIGHTED MODE: EVERY SEQUENCE WEIGHTS 1\n");
      return fp;
    }

  fprintf ( fp, "\n\nWEIGHTED MODE:%s\n\n", (W)->mode);
  if (W->nseq>MAX_NSEQ_4_DISPLAY)return fp;
  for ( a=0, max_len=0; a< W->nseq; a++)max_len=MAX(max_len, strlen (W->seq_name[a]));
  for ( a=0; a< (W->nseq); a++)
    {
      fprintf ( fp, "\t%*s %.2f\n", max_len,(W)->seq_name[a],W->SEQ_W[a]);
    }
  fprintf ( fp, "\n");
  return fp;
}

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT SEQ                                                */
/*                                                                                         */
/***************************************************************************************** */
int ** input_similarities (char *file, Alignment *A, char *mode)
{
  int a, b, i, n;
  int **sim;
  float score;
  char name[1000];
  FILE *fp=NULL;
  char *buf1=NULL, *buf2=NULL;
  int new_aln=0;



  if ( !check_file_exists (file) || !is_distance_matrix_file (file) ||!is_similarity_matrix_file (file) )
    {
      return NULL;
    }

  if ( A)
    {
      fp=vfopen (file, "r");
      while ((buf2=vfgets (buf1,fp))!=NULL )
	{
	  if (strstr (buf2, "SEQ_INDEX"))
	    {
	      buf1=buf2;
	      sscanf (buf1, "# SEQ_INDEX %s %d",name, &i);
	      if ( !strm (A->name[i], name))
		{
		  return NULL;
		}
	    }
	}
      vfclose (fp);
    }
  else
    {

      A=similarities_file2aln(file);
      new_aln=1;
    }

  sim=declare_int ( A->nseq, A->nseq);
  for ( a=0; a<A->nseq; a++)sim[a][a]=100;


  fp=find_token_in_file (file, NULL, "PW_SEQ_DISTANCES");
  fp=find_token_in_file (file, fp, "BOT");
  while ((buf2=vfgets (buf1,fp))!=NULL )
    {
      if ( !(strstr (buf2, "BOT\t") || strstr (buf2, "TOP\t")))continue;
      buf1=buf2;
      n=sscanf (buf1, "%*s %d %d %f", &a, &b, &score);
      if ( n!=3)
	{
	  free_int (sim, -1);
	  return NULL;
	}
      else sim[a][b]=sim[b][a]=(int)score;
    }
  vfclose (fp);
  vfree (buf1);
  if (new_aln)free_aln(A);
  return sim;
}

Alignment * similarities_file2aln ( char *file)
{
  int nseq=0, i;
  FILE *fp;
  char name[1000];
  Alignment *A;


  fp=vfopen (file, "r");
  while ((fp=find_token_in_file (file,fp, "SEQ_INDEX")))nseq++;
  A=declare_aln2 (nseq+1, 10);

  while ((fp=find_token_in_file (file,fp, "SEQ_INDEX")))
    {
      fscanf (fp, "%s %d", name,&i);
      sprintf ( A->name[i], "%s", name);
    }
  A->nseq=nseq;

  return A;
}

void output_similarities (char *file, Alignment *A, char *mode)
{
  float s;
  float *tot;
  double bigtot=0;
  double bigtotsq=0;
  int n, max;
  FILE *fp;
  int a, b;
  char *p;
  int **M=NULL;
  for (max=0, a=0; a< A->nseq; a++)max=MAX(max,(strlen (A->name[a])));


  tot=(float*)vcalloc ( A->nseq, sizeof (float));
  fp=vfopen (file, "w");
  if (!strstr (mode, "avg")&& !strstr (mode, "std"))
    {
      fprintf (fp, "# TC_SIMILARITY_MATRIX_FORMAT_01\n");
      for ( a=0; a<A->nseq; a++)
	fprintf ( fp, "# SEQ_INDEX %s %d\n",A->name[a],a);
      fprintf ( fp, "# PW_SEQ_DISTANCES \n");
    }
  for (n=0,a=0;a< A->nseq-1; a++)
    {
      for ( b=a+1; b<A->nseq; b++, n++)
	{
	   if (strstr (mode, "_sarmat2"))
	    {
	      s=get_sar_sim (A->seq_al[a], A->seq_al[b]);
	    }
	  else if (strstr (mode, "_sar"))
	    {
	      s=get_sar_sim (A->seq_al[a], A->seq_al[b]);
	    }
	  else if ( (p=strstr (mode, "_memory_")))
	    {
	      int **sim;
	      sscanf ( p, "_memory_%ld", (long int*)&sim);
	      s=sim[a][b];
	    }
	  else if ( strstr (mode, "_idscore") || strstr ( mode, "_covscore"))
	    {
	      static Sequence *S;
	      if (a==0 && b==1)
		{
		  free_sequence (S, -1);
		  if ( strstr (mode, "idscoreDNA"))
		    M=read_matrice ("idmat");
		  else
		    M=read_matrice("blosum62mt");

	          S=aln2seq(A);
		}
	      if ( strstr (mode, "_idscore"))s=idscore_pairseq(S->seq[a], S->seq[b], -10,-1, M, "sim");
	      else 	      s=idscore_pairseq(S->seq[a], S->seq[b], -10,-1, M, "cov");
	    }
	  else if ( strstr (mode, "cov"))
	    {
	      s=get_seq_sim ( A->seq_al[a], A->seq_al[b],GAP_LIST, "cov");
	    }
	  else
	    {
	      s=get_seq_fsim2 (A->seq_al[a], A->seq_al[b],GAP_LIST, mode);
	    }
	   if ( !strstr (mode, "avg")&& !strstr (mode, "std"))
	     {
	       fprintf (fp, "BOT\t %4d %4d\t %5.2f %*s\t %*s\t %5.2f\n", a,b,s,max,A->name[a], max, A->name[b], s);
	       fprintf (fp, "TOP\t %4d %4d\t %5.2f %*s\t %*s\t %5.2f\n", b,a,s,max,A->name[b], max, A->name[a], s);
	     }
	  tot[a]+=s;
	  tot[b]+=s;
	  bigtot+=(double)s;
	  bigtotsq+=(double)(s*s);
	}
    }
  for ( a=0; a< A->nseq; a++)
    {
      if ( !strstr (mode, "avg") && !strstr (mode, "std"))fprintf (fp, "AVG\t %d\t %*s\t %*s\t %5.2f\n", a,max,A->name[a], max, "*", tot[a]/(A->nseq-1));

    }
  vfree (tot);free_int (M, -1);
  if ( strstr (mode, "avg"))
    {
      fprintf (fp, "%5.2f\n",(float)(bigtot/(double)n));
    }
  else if (strstr (mode, "std"))
    {
      fprintf (fp, "%5.2f\n",(float)sqrt(((bigtotsq-(bigtot * bigtot)/(double)n)/(double)(n-1))));
    }
  else
    {
      fprintf (fp, "TOT\t %*s\t %*s\t %5.2f\n", max,"TOT", max, "*", (float)(bigtot/(double)n));
      fprintf (fp, "AVG\t %*s\t %*s\t %5.2f\n", max,"AVG", max, "*", (float)(bigtot/(double)n));
      fprintf (fp, "VAR\t %*s\t %*s\t %5.2f\n", max,"VAR", max, "*", (float)(bigtotsq-(bigtot * bigtot)/(double)n)/(double)(n-1));
      fprintf (fp, "STD\t %*s\t %*s\t %5.2f\n", max,"STD", max, "*", (float)sqrt(((bigtotsq-(bigtot * bigtot)/(double)n)/(double)(n-1))));
    }
  vfclose (fp);
}

void output_similarities_pw (char *file, Alignment *A, Alignment *B,char *mode)
{
  float s;
  float *tot;
  float bigtot=0;
  int n, max;
  FILE *fp;
  int a, b;

  int **M=NULL;
  Sequence *SA, *SB;

  if ( strstr (mode, "idscoreDNA"))
    M=read_matrice ("idmat");
  else
    M=read_matrice("blosum62mt");

  SA=aln2seq(A);
  SB=aln2seq(B);

  for (max=0, a=0; a< A->nseq; a++)max=MAX(max,(strlen (A->name[a])));
  for (a=0; a< B->nseq; a++)max=MAX(max,(strlen (B->name[a])));


  tot=(float*)vcalloc ( A->nseq, sizeof (float));
  fp=vfopen (file, "w");
  fprintf (fp, "# TC_SIMILARITY_MATRIX_FORMAT_01\n");
  for ( a=0; a<A->nseq; a++)
    fprintf ( fp, "# SEQ_INDEX %s %d\n",A->name[a],a);
  fprintf ( fp, "# PW_SEQ_DISTANCES \n");
  for (n=0,a=0;a< A->nseq; a++)
    {
      for ( b=0; b<B->nseq; b++, n++)
	{
	  s=idscore_pairseq(SA->seq[a], SB->seq[b], -10,-1, M, "sim");
	  fprintf (fp, "BOT\t %4d %4d\t %5.2f %*s\t %*s\t %5.2f\n", a,b,s,max,A->name[a], max, B->name[b], s);
	  fprintf (fp, "TOP\t %4d %4d\t %5.2f %*s\t %*s\t %5.2f\n", b,a,s,max,B->name[b], max, A->name[a], s);
	  tot[a]+=s;
	  tot[b]+=s;
	  bigtot+=s;
	}
    }

  for ( a=0; a< A->nseq; a++)
    {
      fprintf (fp, "AVG\t %d\t %*s\t %*s\t %5.2f\n", a,max,A->name[a], max, "*", tot[a]/(A->nseq-1));
    }
  vfree (tot);free_int (M, -1);
  fprintf (fp, "TOT\t %*s\t %*s\t %5.2f\n", max,"TOT", max, "*", bigtot/n);
  vfclose (fp);
}
void output_conservation_statistics ( char *file, Alignment *A)
{
  int a, b, c,c1, c2;
  double **tot;
  char aa[1000];
  int naa;

  sprintf (aa, "%s", BLAST_AA_ALPHABET);
  naa=strlen (aa);

  tot=declare_double (256, 256);


  for ( a=0; a<A->nseq; a+=2)
    {
      b=a+1;
      for ( c=0; c<A->len_aln; c++)
	{
	  c1=tolower (A->seq_al[a][c]);
	  c2=tolower (A->seq_al[b][c]);
	  if ( !is_gap(c1) && !is_gap(c2))
	    {
	      tot[c1][c2]++;
	      tot[c2][c1]++;
	      tot[c1][0]++;
	      tot[c2][0]++;
	      tot[0][0]++;
	    }
	}
    }

  fprintf ( stdout, "# BLAST_MATRIX FORMAT\n#ALPHABET=%s\n",aa);
  for (a=0; a<naa; a++)fprintf ( stdout, "%3c ", toupper(aa[a]));
  fprintf ( stdout, "\n");
  for (a=0; a< naa; a++)
    {
      fprintf (stdout, "%c", toupper(aa[a]));
      for ( b=0; b< naa; b++)
	{
	  float f1, f2, f3, r, v;
	  c1=tolower(aa[a]);c2=tolower(aa[b]);
	  f1=(float)((tot[c1][c2]*2)/tot[0][0]);
	  f2=(float)((tot[c1][0])/tot[0][0]);
	  f3=(float)((tot[c2][0])/tot[0][0]);
	  r=(float)(f2==0 || f3==0)?0:(f1/(f2*f3));
	  v=(r==0)?0:((float)10*log((double)r));
	  fprintf (stdout, " %5d",(int)v);
	}
      fprintf ( stdout, "\n");
    }
}
void output_statistics (char *file, Alignment *A, char *mode)
    {
      FILE *fp;
      int a, b, c, d=0, n;
      int maxname=0;


      if (!mode || !mode[0])
	mode="hnrglNL";
      else if ( mode[0]=='_')
	mode++;
      for ( a=0; a<A->nseq; a++)maxname=MAX(strlen(A->name[a]), maxname);
      maxname++;


      fp=vfopen (file, "w");

      if (mode[0]=='h')
	{
	  b=0;
	  while ((c=mode[b++])!='\0')
	    {
	      if ( c=='n') fprintf (fp, "%-*s ",maxname,"name");
	      if ( c=='l') fprintf (fp, "%-*s ",5,"nres");
	      if ( c=='g') fprintf (fp, "%-*s ",5,"ngap");
	      if ( c=='t') fprintf (fp, "%-*s ",5,"len");
	    }
	  if (is_in_set ( c, "nlgt"))	  fprintf (fp, "\n");
	  mode++;
	}
      b=0;
      while ((c=mode[b++])!='\0')
	{
	  if ( c=='n')break;
	  if ( c=='N'){d=1;fprintf (fp, "NSEQ %d ", A->nseq);}
	  if ( c=='L'){d=1;fprintf (fp, "LEN  %d ", A->len_aln);}
	}
      if ( d) fprintf (fp, "\n");

      for (a=0; a<A->nseq; a++)
	{
	  b=0;
	  d=0;
	  while ((c=mode[b++])!='\0')
	    {
	      if (is_in_set ( c, "nlgt"))d=1;

	      if (c=='n'){d=1;fprintf ( fp, "%-*s ", maxname,A->name[a]);}
	      if (c=='l')
		{
		  for (n=0,d=0; d<A->len_aln; d++)n+=!is_gap(A->seq_al[a][d]);
		  fprintf ( fp, "%-5d ",n);
		}
	      if (c=='g')
		{
		  for (n=0,d=0; d<A->len_aln; d++)n+=((is_gap(A->seq_al[a][d]) && !is_gap(A->seq_al[a][d+1]))||(is_gap(A->seq_al[a][d])&& A->seq_al[a][d+1]=='\0')) ;
		  fprintf ( fp, "%-5d ",n);
		}
	      if (c=='t')
		{
		  fprintf ( fp, "%-5d ",(int)strlen (A->seq_al[a]));
		}
	       if (c=='N' && d)
		{
		 fprintf ( fp, "%-5d ",A->nseq);
		}
	      if (c=='L'&& d)
		{
		 fprintf ( fp, "%-5d ",A->len_aln);
		}
	    }
	  if (d)fprintf ( fp, "\n");
	}
      vfclose (fp);
    }

int output_age_matrix ( char *outfile, int val)
{
  int **mat;
  int a, b;
  char alp[]="abcdefghij-";
  int naa;

  mat=declare_int ( 256, 256);
  naa=strlen (alp);
  for ( a=0; a<naa; a++)
    for ( b=0; b<naa; b++)
      {
	if (is_gap(alp[a]) ||is_gap(alp[b] ))mat[(int)alp[a]][(int)alp[b]]=((val==0)?1:val)*-1;
	else mat[(int)alp[a]][(int)alp[b]]=(FABS((a-b))*-1)*((val==0)?1:val);

      }
  output_mat ( mat,outfile, alp, 0);
  free_arrayN((void**)mat, 2);
  return 1;
}




int output_transitions(char *outfile, Alignment *A)
{
  double table[256][256];
  double symbols[256];
  double tot, l, freq, expected, log_odd;
  int a, b;
  char *s;
  char *alp;
  int naa=0;
  int **mat;
  float **fmat;

  FILE *fp;

  for ( a=0; a< 256; a++)
    for (b=0; b<256; b++)
      {
	symbols[b]=0;
	table[a][b]=0;
      }
  alp=(char*)vcalloc ( 256, sizeof (char));
  mat=declare_int ( 256,256);
  fmat=declare_float ( 256,256);

  for (tot=0,a=0; a< A->nseq; a++)
    {
      ungap (A->seq_al[a]);
      lower_string (A->seq_al[a]);
      s=A->seq_al[a];
      l=strlen (s);
      if ( s[0]=='\0') continue;
      symbols[(int)s[0]]++;
      for ( b=1; b< l; b++)
	{
	  symbols[(int)s[b]]++;
	  table[(int)s[b-1]][(int)s[b]]++;
	  tot++;
	}
    }
  for (naa=0, a=0; a< 256; a++)
    {
      if (symbols[a])alp[naa++]=a;
    }


  for ( a=0; a< 256; a++)
    for (b=0; b<256; b++)
      {
	if (symbols[a]&& symbols[b] && table[a][b] && tot>0)
	  {
	    freq=(table[a][b])/tot;
	    expected=(symbols[a]*symbols[b])/(tot*tot);
	    log_odd=log (freq/expected);
	    mat[a-'A'][b-'A']=log_odd*10;
	    fmat[a-'A'][b-'A']=log_odd;
	  }
	else if ( symbols[a]&& symbols[b])
	  {
	    mat[a-'A'][b-'A']=-999;
	    fmat[a-'A'][b-'A']=-999;
	  }
      }
  output_mat ( mat,outfile, alp, 'A');

  fp=vfopen (outfile, "a");
  for ( a=0; a<256; a++)
    if ( symbols[a])
      {
	fprintf (fp, "# %c tot: %6d freq: %7.5f\n", a, (int)symbols[a],(float)symbols[a]/tot);
      }

  for ( a=0; a< 256; a++)
    for (b=0; b<256; b++)
      {
	if (symbols[a]&& symbols[b])
	  {
	    freq=(table[a][b])/tot;
	    fprintf (fp, "# %c%c tot: %6d freq: %7.5f log_odd: %9.3f\n", a, b, (int)table[a][b],(float)freq,fmat[a-'A'][b-'A']);
	  }
      }
  vfclose (fp);
  vfree(alp);
  free_arrayN ((void **)mat, 2);
  free_arrayN ((void **)fmat, 2);

  return 1;
}



void output_est_prf   (char *fname, Alignment *A)
        {
	int a;
	FILE *fp;

	if ( !A->P)
	  {
	    fprintf ( stderr, "\nFormat output_est_prf Impossible: No profile\n");
	    myexit(EXIT_FAILURE);
	  }


	fp=vfopen ( fname, "w");
	fprintf ( fp, "Consensus Sequence\nReconstructed with %s (%s,%s)\n",PROGRAM,AUTHOR,DATE);
	fprintf ( fp, "%4c %4c %4c %4c %15s    Consensus\n",  'A','G','C','T', "Internal Gaps");

	for ( a=0; a< A->len_aln; a++)
	  {
	    fprintf (fp, "%4d %4d %4d %4d %15d %c\n", (A->P)->count[0][a],(A->P)->count[1][a],(A->P)->count[2][a], (A->P)->count[3][a], (A->P)->count[4][a],A->seq_al[0][a]);
	  }
	return;
	}


void output_gotoh_seq (char *fname, Alignment*A )
	{
	int a;
	FILE *fp;

	fp=vfopen ( fname, "w");
	fprintf ( fp, "%d %d\n",A->nseq, A->max_len);
	for ( a=0; a< A->nseq; a++)
		{
		ungap ( A->seq_al[a]);
		fprintf ( fp, ">%s\n", A->name[a]);
		fp=output_string_wrap ( 50,A->seq_al[a] , fp);
		fprintf ( fp, "//\n");
		}

	vfclose (fp);
	}

void output_mult_fasta_seq (char *fname, Alignment*A, int n )
	{
	int a;
	FILE *fp;

	fp=vfopen (fname, "w");
	ungap(A->seq_al[0]);
	for (a=0; a<n; a++)
	  {
	    fprintf (fp, ">%s_%d\n%s\n", A->name[0],a+1, A->seq_al[0]);
	  }
	vfclose (fp);
	}
int output_wexons (char *name, Alignment *A)
{
  int a,b,c;
  FILE *fp;
  int **w;
  fp=vfopen (name, "w");
  if (!fp) return 0;
  if (!A) {vfclose(fp);return 0;}
  w=A->score_res;
  if (!w) {vfclose (fp);return 0; }

  for (a=0; a<A->nseq; a++)
    {
      fprintf (fp, ">%s\n", A->name[a]);
      for (c=0,b=0; b<A->len_aln; b++)
	{
	  int r;
	  r=A->seq_al[a][b];
	  if (!is_gap(r) && r!='f' && r!='F')
	    {
	      fprintf (fp, " %c %d ", r,w[a][c++]);
	    }
	  else if (!is_gap(r))fprintf (fp,"%c -1 ",r);
	}
      fprintf (fp, "\n");
    }
  return 1;
}
char * output_fasta_seqX (char *name, char *mode, Sequence *S, Alignment *A, int i)
{
  FILE *fp;

  if (!name)name=vtmpnam (NULL);
  fp=vfopen (name, mode);
  if ( (S && S->nseq<=i) || (A && S->nseq<=i) || (!A && !S))
    {
      fprintf ( stderr, "\nERROR in function reformat:output_fasta_seqX[FATAL:%s]", PROGRAM);
      myexit (EXIT_FAILURE);
    }

  else if ( S)
    fprintf ( fp, ">%s %s\n%s\n", S->name[i], S->seq_comment[i], S->seq[i]);
  else if ( A)
    {
      ungap (A->seq_al[i]);
      fprintf ( fp, ">%s %s\n%s\n", A->name[i], A->seq_comment[i], A->seq_al[i]);
    }
  vfclose (fp);
  return name;
}
void output_fasta_seq2 (char *fname, Alignment*A )
	{
	char seq_name[VERY_LONG_STRING];
	int a,b;
	FILE *fp;
	char *extension;

	for (a=0; a<A->nseq; a++) ungap (A->seq_al[a]);

	for ( a=0 ; a< A->nseq-1; a++)
	  for ( b=a+1; b< A->nseq; b++)
		{
		  sprintf (seq_name, "%s__%s.fa", A->name[a], A->name[b]);
		  fp=vfopen (seq_name, "w");
		  fprintf (fp, ">%s\n%s\n>%s\n%s\n", A->name[a], A->seq_al[a], A->name[b],A->seq_al[b]);
		  vfclose (fp);
		  fprintf (stdout,">%s\n",seq_name); 
		}
	exit (0);
	}

void output_fasta_seq1 (char *fname, Alignment*A )
	{
	char seq_name[VERY_LONG_STRING];
	int a;
	FILE *fp;
	char *extension;

	for ( a=0; a< A->nseq; a++)
		{
		if ( strncmp( fname, "name",4)==0)
		  {
		    if ( (fname+4)[0]!='\0')extension=fname+5;
		    else
		      extension=NULL;

		     sprintf ( seq_name,"%s.%s", A->name[a],(extension==NULL)?"seq":extension);
		  }
		else
		   sprintf ( seq_name,"%s.seq",A->name[a]);

		ungap ( A->seq_al[a]);
		fp=vfopen (seq_name, "w");
		fprintf (fp, ">%s %s\n", A->name[a], A->seq_comment[a]);
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "\n");
		vfclose (fp);
		}
	}
void output_pir_check (char *fname,int nseq, char **comment )
	{
	int a;
	FILE *fp;

	if ( fname==NULL)return;
	fp=vfopen ( fname, "w");

	for ( a=0; a< nseq; a++)fprintf (fp, "%s\n", comment[a]);
	vfclose (fp);
	}
void output_fasta_seqS (char *fname, Sequence *S)
{
  Alignment *A;
  A=seq2aln (S,NULL,RM_GAP);
  output_fasta_seq (fname, A);
  free_aln (A);
}
void output_fasta_simple (char *fname, Sequence *S)
{
  FILE *fp;
  int a;

  if (!S)return;
  if (!fname)return;

  if (!(fp=vfopen (fname, "w")))return;
  for (a=0; a<S->nseq; a++)fprintf ( fp, ">%s\n%s\n", S->name[a], S->seq[a]);
  vfclose (fp);
}

void output_fasta_seq (char *fname, Alignment*A)
{
  main_output_fasta_seq (fname, A, HEADER);
}

void main_output_fasta_seq (char *fname, Alignment*A,int header )
	{
	int a;
	FILE *fp;

	fp=vfopen ( fname, "w");

	for ( a=0; a< A->nseq; a++)
		{
		ungap(A->seq_al[a]);
		fprintf ( fp, ">%s", A->name[a]);
		if (header==HEADER && A->seq_comment[a][0] && !isblanc(A->seq_comment[a]))fprintf (fp," %s\n",A->seq_comment[a]);
		else fprintf ( fp, "\n");
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "\n");
		}
	vfclose (fp);
	}
void output_gor_seq (char *fname, Alignment*A )
	{
	int a;
	FILE *fp;

	fp=vfopen ( fname, "w");

	for ( a=0; a< A->nseq; a++)
		{
		ungap(A->seq_al[a]);
		fprintf ( fp, "!%s %d \n", A->name[a], (int)strlen(A->seq_al[a]));
		upper_string ( A->seq_al[a]);
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "@\n");
		}
	vfclose (fp);
	}
void output_pir_seq (char *fname, Alignment*A )
	{
	int a;
	for ( a=0; a< A->nseq; a++)ungap(A->seq_al[a]);
	output_pir_aln (fname, A);
	}
void output_pir_seq1 (char *fname, Alignment*A )
	{
	char seq_name[VERY_LONG_STRING];
	int a;
	FILE *fp;
	char type[20];


	for ( a=0; a< A->nseq; a++)
		{
		if      ( strm ( get_string_type (A->seq_al[a]),"DNA") || strm ( get_string_type (A->seq_al[a]),"RNA"))sprintf(type, "DL");
		else if ( strm ( get_string_type (A->seq_al[a]),"PROTEIN"))sprintf(type, "P1");
		sprintf ( seq_name,"%s;%s_%s.seq",type, fname,A->name[a]);
		ungap ( A->seq_al[a]);
		fp=vfopen (seq_name, "w");
		fprintf (fp, ">%s\n\n", A->name[a]);
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "\n*\n");
		vfclose (fp);
		}
	}
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT ALN                                                */
/*                                                                                         */
/***************************************************************************************** */
void output_mocca_aln (char *outfile, Alignment *A, Alignment *S)
    {
    FILE *fp;
    int **score;
    char **new_name_order;
    int a, maxl;

    score=declare_int (S->nseq, 2);
    new_name_order=declare_char ( S->nseq,MAXNAMES+1);
    for ( a=0; a<A->nseq; a++)
      {
	score[a][0]=a;
	score[a][1]=S->score_seq[a];
      }
    sort_int_inv (score+1,2,1,0,S->nseq-2);
    for ( a=0; a<A->nseq; a++)
      {
	sprintf ( new_name_order[a], "%s", A->name[score[a][0]]);
      }
    A=reorder_aln (A, new_name_order, A->nseq);

    fp=vfopen (outfile, "w");
    fprintf ( fp, "MOCCA,(%s,%s, C. Notredame)\nSCORE %d\nNSEQ  %d\nLEN   %d\n",VERSION,DATE, A->score_aln, A->nseq, A->len_aln);

    maxl=return_maxlen ( new_name_order, A->nseq);


    for (a=0; a< A->nseq; a++)
      {
	fprintf (fp, "%-*s: %3d\n", maxl, A->name[a], score[a][1]);
      }

    fprintf ( fp, "\n");

    fp=output_Alignment_without_header ( A, fp);
    vfclose (fp);
    free_int  (score, -1);
    free_char (new_name_order, -1);
    return ;
    }

void print_sub_aln ( Alignment *B, int *ns, int **ls)
{
  Alignment *X;
  int a, b;


  X=copy_aln (B, NULL);
  X->nseq=0;
  X->len_aln=strlen ( B->seq_al[ls[0][0]]);


  for (a=0; a< 2; a++)
    for ( b=0; b<ns[a]; b++, X->nseq++)
      {
	sprintf ( X->seq_al[X->nseq], "%s", B->seq_al[ls[a][b]]);
	sprintf ( X->name[X->nseq], "%s", B->name[ls[a][b]]);
      }
  X->name[X->nseq][0]='\0';

  print_aln (X);
  free_aln (X);
}
void print_aln ( Alignment *B)
    {

    while(B)
      {
	output_Alignment_without_header ( B, stderr);
	B=B->A;
      }
    }


FILE * output_aln ( Alignment *B, FILE *fp){return output_Alignment(B, fp);}
FILE * output_Alignment ( Alignment *B, FILE *fp)
    {
      fprintf ( fp, "%s, %s (%s) [%s] [MODE: %s]\n%s\nCPU   %d sec\nSCORE %d\nNSEQ  %d\nLEN   %d\n",PROGRAM,VERSION,DATE,retrieve_mode(),URL,AUTHOR,  (B->cpu+get_time())/1000, B->score_aln, B->nseq, B->len_aln);

      return output_Alignment_without_header ( B, fp);
    }

FILE * output_Alignment_without_header ( Alignment *B, FILE *fp)
    {
    int a,b, c;
    int max_len=0;
    int line;
    int *n_residues;
    char s;


    if (fp==NULL)return fp;
    for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
    max_len=MAX(max_len+2, 16);
    line=get_msa_line_length (0, 0);
    n_residues=(int*)vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a<B->nseq; a++)n_residues[a]=(B->output_res_num==2)?B->order[a][1]:0;




  fprintf ( fp, "\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<=B->nseq; b++)
	     {
	       fprintf (fp,"%-*s",max_len,B->name[b]);
	       if (B->output_res_num)fprintf (fp, " %4d ", n_residues[b]+1);
	       for (c=a;c<a+line && c<B->len_aln;c++)
		 {
		   if (b==B->nseq){n_residues[b]++;s=analyse_aln_column ( B, c);}
		   else
		     {n_residues[b]+=!is_gap(B->seq_al[b][c]);
		       s=GET_CASE(B->residue_case, B->seq_al[b][c]);
		     }

		   fprintf (fp,"%c",s );
			        }
	       if (B->output_res_num)fprintf (fp, " %4d", n_residues[b]);
	       fprintf (fp,"\n");
	     }

	     fprintf (fp,"\n");
	   }

     fprintf (fp,"\n\n");
     vfree (n_residues);

     return fp;
    }
FILE * output_aln_score ( Alignment *B, FILE *fp){return output_Alignment_score(B, fp);}
FILE * output_Alignment_score ( Alignment *B, FILE *fp)
    {
    int a, b, c;
    static int max_len=0;
    static int line;
    int ch;

    if (fp==NULL)return fp;
    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	max_len+=4;

	}
   line=get_msa_line_length(0, 0);
   sprintf (B->name[B->nseq], "CONS");
   fprintf ( fp, "T_COFFEE ALIGNMENT\nCPU TIME:%d sec.\n", (B->cpu+get_time())/1000);
   fprintf ( fp, "SCORE=%d\n", B->score_aln);
   for ( a=0;a<B->nseq; a++)fprintf ( fp, "%s: %d\n", B->name[a], B->score_seq[a]);
   fprintf ( fp, "\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	      {
	      fprintf (fp,"%-*s",max_len,B->name[b]);
	      for (c=a;c<a+line && c<B->len_aln;c++)
		{
		ch=B->seq_al[b][c];
		if (ch==NO_COLOR_RESIDUE)fprintf (fp,"-");
		else if ( ch==NO_COLOR_GAP)fprintf (fp,"*");
		else if ( ch<10 && ch>=0)fprintf (fp,"%d",ch);
		else if ( ch>10)fprintf (fp,"#");
		else if ( ch<0)fprintf  (fp,".");
		else fprintf (fp,"9");
		}
	      fprintf (fp,"\n");
	      }
	    fprintf (fp,"\n");
	    fprintf (fp,"%-*s",max_len,B->name[b]);
	    for (c=a;c<a+line && c<B->len_aln;c++)
	      {
	      ch=B->seq_al[b][c];
	      if (ch==NO_COLOR_RESIDUE)fprintf (fp,"-");
	      else if ( ch==NO_COLOR_GAP)fprintf ( fp, "*");
	      else if ( ch<10 && ch>=0)fprintf (fp,"%d",ch);
	      else if ( ch>10)fprintf (fp,"#");
	      else if ( ch<0)fprintf (fp,".");
	      else fprintf (fp,"9");
	      }
	    fprintf (fp,"\n\n\n");
	   }
    fprintf (fp,"\n\n");
    return fp;
    }
FILE * output_aln_with_res_number ( Alignment *B, FILE *fp){return  output_Alignment_with_res_number(B, fp);}
FILE * output_Alignment_with_res_number ( Alignment *B, FILE *fp)
    {
    int a, b, c;
    static int max_len=0;
    static int line;
    int**order;

    if (fp==NULL)return fp;
    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	max_len+=4;
	line=60;
	}
    order=copy_int ( B->order,declare_int ( B->nseq, 2));

   fprintf ( fp, "T_COFFEE ALIGNMENT\nCPU TIME:%d sec.\n", (B->cpu+get_time())/1000);
   fprintf ( fp, "\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {
	      fprintf (fp,"%-*s %3d %4d ",max_len,B->name[b], order[b][0], order[b][1] );
	      for (c=a;c<a+line && c<B->len_aln;c++)
		{
		order[b][1]+=1-is_gap(B->seq_al[b][c]);
		fprintf (fp,"%c",toupper(B->seq_al[b][c]) );
		}
	      fprintf (fp," %4d\n", order[b][1] );
	      }
	    fprintf (fp,"\n");
	    }
    fprintf (fp,"\n\n");

    free_int (order, -1);
    return fp;
    }

void output_constraints ( char *fname, char *mode,Alignment *A)
	{
	FILE *fp;
	Constraint_list *CL;
	char *buf;
	char **name_list;

	if ( !A->CL || strm ( mode, "pdb"))
	   {
	       if (!A->S)
	          {
		      A->S=aln2seq(A);
		  }

	       CL=declare_constraint_list ( A->S, NULL, NULL, 0, NULL, NULL);
	       CL=aln2constraint_list (A,CL, mode);
	       fp=save_constraint_list ( CL, 0, CL->ne,fname, NULL, "lib",A->S);
	       vfclose (fp);
	       free_constraint_list (CL);
	       return;
	   }
	else if ( strncmp ( mode, "extended_pair", 13)==0)
	  {
	    buf=duplicate_string (mode+14);

	    name_list=(char**)vcalloc(2, sizeof(char*));
	    name_list[0]=strtok (buf,"_");
	    name_list[1]=strtok (NULL,"_");
	    mode[13]='\0';


	    CL=A->CL;
	    fp=save_sub_list_header (vfopen(fname, "w"),2, name_list,CL);
	    fp=save_extended_constraint_list_pair (CL, "pair",name_list[0],name_list[1],fp);
	    fp=save_list_footer (fp, CL);
	    vfree (buf);
	  }
	else if ( strm2 (mode, "extended_lib","extended_cosmetic"))
	  {
	    CL=A->CL;
	    fp=save_extended_constraint_list ( CL,mode+9, vfopen(fname, "w"));
	  }
	else
	   {
	       CL=(Constraint_list *)A->CL;
	       fp=save_constraint_list ( CL, 0, CL->ne,fname, NULL, "lib",A->S);
	   }
	vfclose ( fp);

	if ( (Constraint_list *)A->CL !=CL)free_constraint_list (CL);

	return;

	}
void output_model_aln (char *fname, Alignment*A )
        {
	  FILE *fp;
	  int a;
	  Dp_Model *M;
	  Dp_Result *R;
	  char *string;

	  if ( A->Dp_result==NULL)
	    {
	      fprintf ( stderr, "\nWARNING Could Not Output Model %s [%s]", fname, PROGRAM);
	    }
	  R=A->Dp_result;
	  M=R->Dp_model;

	  fp=vfopen ( fname, "w");
	  for (a=0; a<M->nstate; a++)
	    {
	      if (M->model_comments[a][0])fprintf ( fp, "#STATE %c: %s\n", 'a'+a, M->model_comments[a]);
	    }
	  string=(char*)vcalloc ( R->len+1, sizeof (char));
	  for (a=0; a<R->len; a++)string[a]=R->traceback[a]+'a';
	  fprintf ( fp, ">%s\n",fname);
	  fp=output_string_wrap ( 50,string, fp);
	  vfree(string);
	  fprintf ( fp, "\n");

	  vfclose (fp);
	  return;
	}
char * output_fasta_sub_aln (char *fname, Alignment*A, int ns, int *ls  )
{
  int a,s;
  FILE *fp;
  if (fname==NULL)fname=vtmpnam (NULL);
  fp=vfopen (fname, "w");
  for (a=0; a<ns; a++)
    {
      s=ls[a];
      fprintf (fp, ">%s %s\n%s\n", A->name[s],A->seq_comment[s],A->seq_al[s]);
    }
  vfclose (fp);
  return fname;
}
char * output_fasta_sub_aln2 (char *fname, Alignment*A, int *ns, int **ls  )
{
  int a,g,s;
  FILE *fp;
  if (fname==NULL)fname=vtmpnam (NULL);
  fp=vfopen (fname, "w");
  for ( g=0; g<2; g++)
    for (a=0; a<ns[g]; a++)
      {
	s=ls[g][a];
	fprintf (fp, ">%s %s\n%s\n", A->name[s],A->seq_comment[s],A->seq_al[s]);
      }
  vfclose (fp);
  return fname;
}

int output_suchard_aln (char *out_file, Alignment *A)
{
  int a, b, c, d;
  FILE *fp;

  A=back_translate_dna_aln (A);

  for ( c=0,a=0; a<A->len_aln; a++, c++)
	     {
	       if (c==3)c=0;
	       for (b=0; b<A->nseq; b++)
		 {
		 if (c==2)
		   {
		     A->seq_al[b][a]='-';
		   }
		 }
	     }
  A=ungap_aln_n (A, 1);
  fp=vfopen (out_file, "w");
  for ( a=0; a< A->nseq; a++)
    {
      for (b=0; b< A->len_aln; b++)
	{
	  c=tolower(A->seq_al[a][b]);
	  if ( c=='a')d=1;
	  else if ( c=='g')d=2;
	  else if ( c=='c')d=3;
	  else if ( c=='t')d=4;
	  else if ( c=='u')d=5;
	  else d=6;

	  fprintf ( fp, "%d", d);
	}
      fprintf ( fp, "\n");
    }
  vfclose (fp);
  myexit (EXIT_SUCCESS);
}
void output_mfasta_aln (char *fname, Alignment *A )
{
	FILE *fp;
	int a,b,c,l;
	int line=0;
	static int rep;
	
	line=get_msa_line_length (line, A->len_aln+1);
	if (check_file_exists(fname)) 
	  {
	    fp= vfopen ( fname, "a");//make it possible to output replicates
	    fprintf (fp, "#Replicate %d\n", rep);
	    rep ++;
	  }
	else
	  {
	    fp= vfopen ( fname, "w");
	    rep=1;
	  }
	for ( a=0; a< A->nseq; a++)
	{
		fprintf ( fp, ">%s", A->name[a]);
		if ( A->seq_comment[a][0] && !isblanc (A->seq_comment[a]))
			fprintf ( fp, " %s", A->seq_comment[a]);
		fprintf ( fp, "\n");
		l=strlen (A->seq_al[a]);
		for (b=0, c=0;b<l; b++, c++)
		  {
		    if (line>0 && c==line){fprintf (fp, "\n");c=0;}
		    fprintf (fp, "%c", GET_CASE(A->residue_case, A->seq_al[a][b]));
		  }
		fprintf ( fp, "\n");
	}
	vfclose (fp);
}
void output_fasta_aln (char *fname, Alignment *A )
{
	FILE *fp;
	int a,b,c,l;
	int line=0;

	line=get_msa_line_length (line, A->len_aln+1);
	fp=vfopen ( fname, "w");

	for ( a=0; a< A->nseq; a++)
	{
		fprintf ( fp, ">%s", A->name[a]);
		if ( A->seq_comment[a][0] && !isblanc (A->seq_comment[a]))
			fprintf ( fp, " %s", A->seq_comment[a]);
		fprintf ( fp, "\n");
		l=strlen (A->seq_al[a]);
		for (b=0, c=0;b<l; b++, c++)
		  {
		    if (line>0 && c==line){fprintf (fp, "\n");c=0;}
		    fprintf (fp, "%c", GET_CASE(A->residue_case, A->seq_al[a][b]));
		  }
		fprintf ( fp, "\n");
	}
	vfclose (fp);
}



/**
 * Prints the alignment in XMFA format_is_conc_aln
 * \param fname The file to print the output to.
 * \param A	The alignment to print.
*/
void
output_xmfa_aln(char *fname, Alignment *A)
{
	if (A->S->genome_co == NULL)
	{
		fprintf(stderr, "WARNING: Genomic information not available. Producing fasta file instead!\n");
		output_fasta_aln (fname, A );
		return;
	}
	FILE *fp;

	int line=0;
	line=get_msa_line_length (line, A->len_aln+1);
	fp=vfopen ( fname, "w");
	Genomic_info *genome_co = A->S->genome_co;
	unsigned int aln_index = 0;
	unsigned int seq_index;
	unsigned int num_seqs = A->S->nseq;
	char **names = A->S->name;
	for (; aln_index< A->nseq; ++aln_index)
	{
		seq_index = name_is_in_list(A->name[aln_index], names, num_seqs, 0);
		fprintf ( fp, ">%i:%i-%i %c %s %s", aln_index+1, genome_co[seq_index].start+1, genome_co[seq_index].end+1, genome_co[seq_index].strand, genome_co[seq_index].seg_name, A->name[aln_index]);
		if ( A->seq_comment[aln_index][0] && !isblanc (A->seq_comment[aln_index]))
			fprintf ( fp, " %s", A->seq_comment[aln_index]);

		fprintf ( fp, "\n");
		fp=output_string_wrap ( line,A->seq_al[aln_index] , fp);
		fprintf ( fp, "\n");
	}
	fprintf(fp,"=\n");
	vfclose (fp);

}





void output_pir_aln (char *fname, Alignment*A )
	{
	int a;
	FILE *fp;
	char type[20];





	fp=vfopen ( fname, "w");
	for ( a=0; a< A->nseq; a++)
		{
		if      ( strm ( get_string_type (A->seq_al[a]),"DNA") || strm ( get_string_type (A->seq_al[a]),"RNA"))sprintf(type, "DL");
		else if ( strm ( get_string_type (A->seq_al[a]),"PROTEIN"))sprintf(type, "P1");
		fprintf ( fp, ">%s;%s\n%s\n",type, A->name[a], A->seq_comment[a]);
		fp=output_string_wrap ( 50,A->seq_al[a] , fp);
		fprintf ( fp, "\n*\n");
		}

	vfclose (fp);
	}

int landscape_msa;
int  set_landscape_msa (int len)
{
  if ( len==0)landscape_msa=-1;
  else
    {
      landscape_msa=len;
    }
  return landscape_msa;
}
int get_msa_line_length (int line, int aln_len)
{
  if (landscape_msa==-1) return aln_len;
  else if ( landscape_msa)return landscape_msa;
  else if (line) return line;
  else
    {
      return (getenv ("ALN_LINE_LENGTH"))?atoi(getenv("ALN_LINE_LENGTH")):ALN_LINE_LENGTH;
    }
}

void output_msf_aln (char *fname,Alignment *B)
        {
	int a, b, c;
	char *seq;
	int *all_checks;
	int i,j;
	long grand_checksum;
	FILE *fp;
	int max_len;
	int line=0;
	int block=10;
	int c_block;
	char aa;


	line=get_msa_line_length (line, B->len_aln+1);


	for ( max_len=0,a=0; a< B->nseq; a++)max_len= MAX(strlen ( B->name[a]),max_len);


	max_len+=5;

	fp=vfopen (fname, "w");

	seq =(char*)vcalloc(B->len_aln,  sizeof(char));
	all_checks =(int*)vcalloc(B->nseq, sizeof(int));
	for ( i=0; i< B->nseq; i++)
	  {
	    for ( j=0; j<B->len_aln; j++)
	      {
		if ( is_gap(B->seq_al[i][j]))seq[j]='.';
		else seq[j]=B->seq_al[i][j]=toupper(B->seq_al[i][j]);

	      }
	    all_checks[i] = SeqGCGCheckSum(seq, (int)B->len_aln);
	  }
	grand_checksum = 0;
	for(i=0; i<B->nseq; i++) grand_checksum += all_checks[i];
	grand_checksum = grand_checksum % 10000;
	fprintf(fp,"PileUp\n\n");
	B=get_aln_type(B);
	fprintf(fp,"\n\n   MSF:%5d  Type: ",B->len_aln);
	if(strm ( (B->S)->type, "DNA") || strm ( (B->S)->type, "RNA"))
		fprintf(fp,"N");
	else
		fprintf(fp,"P");
	fprintf(fp,"    Check:%6ld   .. \n\n", (long)grand_checksum);
	for (i=0; i< B->nseq; i++)
	  {
	    fprintf ( fp, " Name: %s oo  Len:%5d  Check:%6ld  Weight:  %.3f\n", B->name[i], B->len_aln,(long)all_checks[i],(B->S)->W?((B->S)->W)->SEQ_W[i]:1.00);
	  }
	fprintf(fp,"\n//\n\n");

	for (a=0; a<B->len_aln; a+=line)
	   {
	     fprintf ( fp,"\n\n");
	     for (b=0; b<B->nseq; b++)
	       {
		 fprintf (fp,"%-*s ",max_len,B->name[b]);
		 for (c_block=0,c=a;c<a+line && c<B->len_aln;c++)
		   {
		     if ( c_block==block)
			    {
			      fprintf (fp, " ");
			      c_block=0;
			    }
			c_block++;
		     aa=(is_gap(B->seq_al[b][c]))?'.': toupper(B->seq_al[b][c]);
		     fprintf (fp,"%c",aa );
		   }
		 if ( c_block==block)
			    {
			      fprintf (fp, " ");
			      c_block=0;
			    }
		 fprintf (fp,"\n");

	       }
	   }
    	fprintf ( fp,"\n");
	vfclose ( fp);


	vfree(seq);
	vfree(all_checks);


	return;
}
int SeqGCGCheckSum(char *seq, int len)
{
	int  i;
        long check;

        for( i=0, check=0; i< len; i++,seq++)
                check += ((i % 57)+1) * toupper(*seq);

        return(check % 10000);
}
void old_output_msf_aln (char *fname,Alignment *B)
	{
	FILE *fp;
	static int *put_seq;
	int a, b, c;
	int line=0;
	char aa;
	char *buf;
	int max_len;
    	int seq_max_len;

	line=get_msa_line_length (line, B->len_aln+1);


    	for ( max_len=0,a=0; a< B->nseq; a++)max_len= MAX(strlen ( B->name[a]),max_len);
	for ( seq_max_len=0,a=0; a< B->nseq; a++)seq_max_len= MAX(strlen ( B->seq_al[a]),max_len);


	buf=(char*)vcalloc(seq_max_len+1, sizeof (int));

	if ( put_seq==NULL)
		put_seq=(int*)vcalloc ( B->nseq, sizeof (int));
	put_seq[0]=1;


	for ( b=1; b< B->nseq; b++)
		{
		sprintf ( buf, "%s", B->seq_al[b]);
		ungap(buf);
		put_seq[b]=( strlen (buf)>0)?1:0;
		}

	fp=vfopen ( fname, "w");
	fprintf ( fp, "MSF: %d Type P Check: 5083 ..\n", B->len_aln);
	for ( a=0; a< B->nseq; a++)
		{
		if ( put_seq[a]==1)
			fprintf ( fp,"Name: %s\n",B->name[a]);
		}
	 fprintf ( fp, "//\n");
	 for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {
	     if ( put_seq[b]==1)
	     	{
	     	fprintf (fp,"%-*s ",max_len,B->name[b]);
	        for (c=a;c<a+line && c<B->len_aln;c++)
			{



			aa=(B->seq_al[b][c]=='-')?'.': toupper(B->seq_al[b][c]);
			fprintf (fp,"%c",aa );
			}
	      	fprintf (fp,"\n");
	      	}
	      }
	    fprintf (fp,"\n");
	    }
    	fprintf ( fp,"\n\n");
	vfclose ( fp);

	vfree (buf);
	vfree(put_seq);
	}

void output_saga_aln ( char *name, Alignment *B)
    {
    int a, b, c;
    FILE *fp;



    int max_len;
    int line=0;

    line=get_msa_line_length (line, B->len_aln+1);



    for ( max_len=0,a=0; a< B->nseq; a++)max_len= (strlen ( B->name[a])>max_len)?(strlen ( B->name[a])):max_len;




    fp= vfopen ( name, "w");

    fprintf (fp, "\nSAGA FORMAT\nalignement  %s nseq=%d len=%d\n", name, B->nseq, B->len_aln);

    fprintf (fp, "\n\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {fprintf (fp,"%-*s ",max_len,B->name[b]);
	      for (c=a;c<a+line && c<B->len_aln;c++)
		{
		  fprintf (fp,"%c",(B->seq_al[b][c]) );
		}
	      fprintf (fp,"\n");
	      }
	    fprintf (fp,"\n");
	    }
    fprintf (fp,"\n\n");
    vfclose ( fp);
    }
void output_compact_aln ( char *name, Alignment *B)
    {
    int a, b, c;
    FILE *fp;
    int do_print=0;


    int max_len;
    int line=0;

    line=get_msa_line_length (line, B->len_aln+1);


    for ( max_len=0,a=0; a< B->nseq; a++)max_len= (strlen ( B->name[a])>max_len)?(strlen ( B->name[a])):max_len;




    fp= vfopen ( name, "w");

    fprintf (fp, "\nSAGA FORMAT\nalignement  %s nseq=%d len=%d", name, B->nseq, B->len_aln);
    fprintf (fp, "\n\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {

	     for ( do_print=0, c=a;c<a+line && c<B->len_aln;c++)
	       do_print+=1-is_gap(B->seq_al[b][c]);
	     if ( do_print>0)
	           {
		     fprintf (fp,"%-*s ",max_len,B->name[b]);



		     for (c=a;c<a+line && c<B->len_aln;c++)
		       {
			 if ( is_gap(B->seq_al[b][c])&& B->seq_al[b][c]!='-' )fprintf (fp,"%c", '-');
			 else fprintf (fp,"%c",(B->seq_al[b][c]) );
		       }
		     fprintf (fp,"\n");
		   }
	      }
	    fprintf (fp,"\n");
	    }
    fprintf (fp,"\n\n");
    vfclose ( fp);
    }

void output_clustal_aln ( char *name, Alignment *B)
{
  return output_generic_clustal_aln (name, B, "tc_clustal");
}
void output_strict_clustal_aln ( char *name, Alignment *B)
{
  return output_generic_clustal_aln (name, B, "strict_clustal");
}

void output_generic_clustal_aln ( char *name, Alignment *B, char *mode)
    {
    int a, b, c;
    FILE *fp;
    int max_len=0;
    int line=0;
    int *n_residues;

    if ( getenv ("SEP_4_TCOFFEE"))
      {
	while ( line<B->len_aln && B->seq_al[0][line]!='o' && B->seq_al[0][line]!='O')line++;
	if ( B->seq_al[0][line]=='O' || B->seq_al[0][line]=='o')line++;
      }
    else
      {
	while ( line<B->len_aln)line++;
      }

    if ( line==B->len_aln)line=get_msa_line_length (0, B->len_aln+1);

    n_residues=(int*)vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    n_residues[a]=B->order[a][1];
	    }
    max_len=MAX(max_len+2, 16);


    fp= vfopen ( name, "w");

    if ( strm (mode, "strict_clustal"))
      fprintf ( fp, "CLUSTAL W (1.83) multiple sequence alignment");
    else
      {
	fprintf (fp, "CLUSTAL FORMAT for %s %s [%s] [MODE: %s ], CPU=%.2f sec, SCORE=%d, Nseq=%d, Len=%d ", PROGRAM, VERSION,URL, retrieve_mode (),(float)(B->cpu+get_time())/1000, B->score_aln, B->nseq, B->len_aln);
	if (B->ibit>0)
	  {
	    float ibit=(float)log ((double)B->ibit)/log ((double)2);
	    float nibit=(float)log(ibit/(B->len_aln*B->nseq));
	    fprintf ( fp, " Ties: %.1f bits (%d alternative)\n",ibit, B->ibit-1);

	  }
      }
    fprintf (fp, "\n\n");


    if ( B->len_aln==0)
      {
	for (b=0; b<=B->nseq; b++)
	  fprintf (fp,"%-*s -\n",max_len, B->name[b]);
      }

    else
      {
	for (a=0; a<B->len_aln; a+=line)
	  {for (b=0; b<=B->nseq; b++)
	    {
	      if (b!=B->nseq)
		{
		  fprintf (fp,"%-*s",max_len, B->name[b]);
		  for (c=a;c<a+line && c<B->len_aln;c++)
		    {
		      if ( is_gap(B->seq_al[b][c]))fprintf (fp,"%c", '-');
		      else
			{
			  n_residues[b]++;
			  fprintf (fp, "%c", GET_CASE(B->residue_case, B->seq_al[b][c]));

			}

		    }
		  if (B->output_res_num)fprintf (fp, " %d", n_residues[b]);
		  fprintf (fp,"\n");
		}
	      else if ( b==B->nseq)
		{
		  fprintf (fp,"%-*s",max_len," ");
		  for (c=a;c<a+line && c<B->len_aln;c++)
		    {
		      fprintf ( fp, "%c", analyse_aln_column (B, c));
		    }
		  fprintf (fp,"\n");
		}
	    }
	  fprintf (fp,"\n");
	  }
      }
    fprintf (fp,"\n\n");
    vfree (n_residues);
    vfclose ( fp);
    }
FILE * output_generic_interleaved_aln (FILE *fp, Alignment *B, int line, char gap, char *mode)
    {
    int a, b, c;
    int max_len=0;
    int *n_residues;


    n_residues=(int*)vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    n_residues[a]=B->order[a][1];
	    }
    max_len=MAX(max_len+2, 16);




    if ( B->len_aln==0)
      {
	for (b=0; b<=B->nseq; b++)
	  fprintf (fp,"%-*s -\n",max_len, B->name[b]);
      }

    else
      {
	for (a=0; a<B->len_aln; a+=line)
	  {for (b=0; b<=B->nseq; b++)
	    {
	      if (b!=B->nseq)
		{
		  fprintf (fp,"%-*s",max_len, B->name[b]);
		  for (c=a;c<a+line && c<B->len_aln;c++)
		    {
		      if ( is_gap(B->seq_al[b][c]))fprintf (fp,"%c", gap);
		      else
			{
			  n_residues[b]++;
			  fprintf (fp, "%c", GET_CASE(B->residue_case, B->seq_al[b][c]));

			}

		    }
		  if (B->output_res_num)fprintf (fp, " %d", n_residues[b]);
		  fprintf (fp,"\n");
		}
	    }
	  fprintf (fp,"\n");
	  }
      }
    vfree (n_residues);
    return fp;
    }

static int mnl;
void output_rphylip_aln ( char *name, Alignment *B)
{
  mnl=-1;
  output_phylip_aln (name, B, "w");
  
}

void output_phylip_aln ( char *name, Alignment *B, char *mode)
    {
      int a, b, c, d;
      FILE *fp;
      
      int *print_name;
      static int line=0;
      
      if (!line)
	{
	  if (!getenv("ALN_LINE_LENGTH"))cputenv ("ALN_LINE_LENGTH=60");
	  line=get_msa_line_length(0, 0);
	}
      
      if (mnl==-1)
	{
	  for (a=0; a<B->nseq; a++)
	    {
	      int nl=(B->name && B->name[a])?strlen (B->name[a]):0;
	      if (nl>mnl)mnl=nl;
	    }
	}
      
      


      print_name=(int*)vcalloc ( B->nseq, sizeof (int));
      if (check_file_exists(name)) fp= vfopen ( name, "a");//make it possible to output replicates
      else fp= vfopen ( name, mode);
      
      fprintf (fp, "%5d  %d\n", B->nseq, B->len_aln);
      for (a=0; a<B->len_aln; a+=line)
	{for (b=0; b<B->nseq; b++)
	    {
	      if (!mnl)
		{
		  if ( print_name[b]==0)
		    {
		      
		      fprintf (fp,"%-10.10s ",B->name[b]);
		      print_name[b]=1;
		    }
		  else
		    {
		      fprintf (fp, "%10.10s ", " ");
		    }
		}
	      else
		{
		   if ( print_name[b]==0)
		    {
		      
		      fprintf (fp,"%-*.*s ",mnl,mnl,B->name[b]);
		      print_name[b]=1;
		    }
		  else
		    {
		      fprintf (fp, "%*.*s ", mnl, mnl," ");
		    }
		}
	      
	      for (d=0,c=a;c<a+line && c<B->len_aln;c++, d++)
		{
		  if ( d==10)
		    {
		      fprintf ( fp, " ");
		      d=0;
		    }
		  if ( is_gap(B->seq_al[b][c])&& B->seq_al[b][c]!='-' )fprintf (fp,"%c", '-');
		  else fprintf (fp,"%c",(B->seq_al[b][c]) );
		}
	      fprintf (fp,"\n");
	    }
	  fprintf (fp,"\n");
	}
      fprintf (fp,"\n\n");
      
      vfclose ( fp);
      mnl=0;
    }

void output_rnalign (char *out_file, Alignment *A, Sequence *STRUC)
    {
    int a, b;
    FILE *fp;
    char bank_file[100];
    char pep_file[100];
    char *buf;

    sprintf ( bank_file, "%s.mss", out_file);
    sprintf ( pep_file, "%s.one_rna", out_file);


    buf=(char*)vcalloc ( strlen ( A->seq_al[0]+1), sizeof (char));

    for ( b=0,a=0; a< strlen(A->seq_al[0]); a++)
    	{
    	if ( is_gap(A->seq_al[0][a]))
    		buf[a]='.';
    	else
    		buf[a]=STRUC->seq[0][b++];
    	}
    buf[a]='\0';

    fp=vfopen ( bank_file, "w");

    fprintf ( fp, "ST\n");
    fp=output_string_wrap ( 50, buf, fp);
    fprintf ( fp, "\n\n");

    for ( a=0; a<A->nseq-1; a++)
    	{
    	fprintf ( fp, "AS %s\n ", A->name[a]);
    	fp=output_string_wrap ( 50, A->seq_al[a], fp);
	fprintf ( fp, "\n\n");
	}
    vfclose ( fp);
    fp=vfopen ( pep_file, "w");
    fprintf ( fp, ">%s\n", A->name[A->nseq-1]);
    fp=output_string_wrap ( 50, A->seq_al[A->nseq-1], fp);
    fprintf ( fp, "\n");
    vfclose (fp);
    }

void output_lib (char *pw_lib_saga_aln_name, Alignment *A )
    {
    Alignment *B;
    char fname[VERY_LONG_STRING];
    int a,b;

    B=declare_Alignment (NULL);

    B->nseq=2;

    for ( a=0; a< A->nseq-1; a++)
    	{
    	for ( b=a+1; b<A->nseq; b++)
    		{
    		sprintf ( B->seq_al[0], "%s", A->seq_al[a]);
    		sprintf ( B->name[0], "%s", A->name[a]);
    		sprintf(B->name[1], "%s", A->name[b]);
    		sprintf ( B->seq_al[1], "%s",A->seq_al[b]);
    		B->nseq=2;
    		sprintf ( fname, "%s_%s_%s.lib",pw_lib_saga_aln_name, A->name[a], A->name[b]);

    		B->len_aln=strlen ( B->seq_al[0]);
    		ungap_aln (B);
    		output_clustal_aln (fname,B);
        	}
        }
    }
void output_pw_lib_saga_aln (char *pw_lib_saga_aln_name, Alignment *A )
    {
    Alignment *B;
    char fname[VERY_LONG_STRING];
    int a,b;

    B=declare_Alignment (NULL);

    B->nseq=2;

    for ( a=0; a< A->nseq-1; a++)
    	{
    	for ( b=a+1; b<A->nseq; b++)
    		{
    		sprintf ( B->seq_al[0], "%s", A->seq_al[a]);
    		sprintf ( B->name[0], "%s", A->name[a]);
    		sprintf(B->name[1], "%s", A->name[b]);
    		sprintf ( B->seq_al[1], "%s",A->seq_al[b]);
    		B->nseq=2;
    		sprintf ( fname, "%s_%s_%s.pw_lib_saga_aln",pw_lib_saga_aln_name, A->name[a], A->name[b]);

    		B->len_aln=strlen ( B->seq_al[0]);
    		ungap_aln (B);
    		output_clustal_aln (fname,B);
        	}
        }
    }
void output_lalign_header( char *name, Alignment *A)
    {
    FILE *fp;

    fp=vfopen ( name, "w");
    fprintf ( fp, " Lalign mode: best local alignments between two sequences\n");
    fprintf ( fp, " %s(%s) [%s]\n\n", VERSION, DATE, URL);
    fprintf ( fp, " Comparison of:\n(A) %s\t%s\t-%d aa\n", (A->S)->file[A->order[0][0]],(A->S)->name[A->order[0][0]], (A->S)->len[A->order[0][0]]);
    fprintf ( fp, "(B) %s\t%s\t-%d aa\n", (A->S)->file[A->order[1][0]],(A->S)->name[A->order[1][0]], (A->S)->len[A->order[1][0]]);


    vfclose ( fp);
    return;
    }
void output_stockholm_aln (char *file, Alignment *A, Alignment *ST)
{
  FILE *fp;
  int a,b;

  for (a=0; a<A->nseq; a++)
    for (b=0; b<A->len_aln; b++)
      if (A->seq_al[a][b]==STOCKHOLM_CHAR)A->seq_al[a][b]='.';

  fp=vfopen (file, "w");
  fprintf ( fp, "# STOCKHOLM 1.0\n\n");
  output_generic_interleaved_aln (fp,A, 50, '.', NULL);
  fprintf ( fp, "//\n");
  vfclose (fp);
}

void output_glalign ( char *name, Alignment *B, Alignment *S)
{
  int a, b, g, s;
  int naln=0;
  FILE *fp;
  int **nr;
  B=B->A;
  if ( B==NULL){return;}

  fp=vfopen (name, "w");
  fprintf (fp, "Format: GLALIGN_01 [Generated with %s ]\n", PROGRAM);
  fprintf (fp, "#Each Line corresponds to a column\n");
  fprintf (fp, "#First column coresponds to first genome\n");
  fprintf (fp, "#Last Column gives the column reliability on a 0-9 scale\n");
  fprintf (fp, "#[-1] Indicates that the reliability was not evaluated\n");

  fprintf (fp, "Genome List\n");
  for ( a=0; a< B->nseq; a++)
    fprintf (fp, "\tGenome %s\n", B->name[a]);
  fprintf (fp, "Alignment List\n");
  while (B)
    {
      fprintf (fp, "Alignment %d Len %d Score %d\n", ++naln, B->len_aln, S->score_aln);
      nr=duplicate_int (B->order, -1, -1);
      for ( a=0; a< B->len_aln; a++)
	{
	  fprintf ( fp, "\t");
	  for ( b=0; b< B->nseq; b++)
	    {
	      g=is_gap (B->seq_al[b][a]);
	      nr[b][1]+=1-g;

	      if (g)fprintf (fp, "---- ");
	      else fprintf ( fp, "%4d ",nr[b][1]);
	    }
	  s=((S)?S->seq_al[S->nseq][a]:-1);
	  if (s==NO_COLOR_RESIDUE)s=-1;
	  fprintf ( fp,"[ %d ]",s);
	  fprintf ( fp, "\n");

	}
      free_int (nr, -1);
      B=B->A;
      S=S->A;
    }
  vfclose ( fp);
}
Alignment *input_conc_aln ( char *name, Alignment *IN)
{
  FILE *fp;
  char *string, *p, *file;
  Alignment *F=NULL,*A=NULL, *B=NULL;

  file=vtmpnam (NULL);

  string=file2string(name);
  string=substitute ( string, "@", "!Protected!");
  string=substitute ( string, TC_REC_SEPARATOR, "@");
  strtok (string,"@");


  while ( (p=strtok (NULL,"@"))!=NULL)
    {
      char *buf;
      HERE ("--- %s", p);
      if ( p[0]=='#')continue;
      buf=(char*)vcalloc ( strlen (p)+1, sizeof (char));
      sprintf (buf,"%s", p);
      buf=substitute (buf,"!protected!", "@");

      fp=vfopen (file, "w");
      fprintf ( fp, "%s",buf);
      vfclose (fp);
      vfree (buf);

      if ( is_aln (file))
	{
	  B=main_read_aln (file,NULL);

	  if ( !A)
	    {
	      if (IN){copy_aln (B, IN);F=A=IN;}
	      else F=A=B;
	    }
	  else
	    {
	      A->A=B;
	      A=A->A;
	    }
	}
    }

  vfree (string);
  return F;
}

void output_conc_aln ( char *name, Alignment *B)
{
  FILE *fp;
  int a;

  fp=vfopen (name, "w");
  fprintf (fp, "# CONC_MSF_FORMAT_01\n");
  while (B)
    {
      fprintf (fp, "%s\n", TC_REC_SEPARATOR);
      for ( a=0; a< B->nseq; a++)
	{
	  fprintf ( fp, ">%s\n%s\n", B->name[a], B->seq_al[a]);
	}
      B=B->A;

    }
  vfclose (fp);
}

void output_lalign ( char *name, Alignment *B)
{
  static int output_header;

  B=B->A;
  if ( B==NULL){output_header=0;return;}
    else if ( output_header==0)
      {
	output_lalign_header(name, B);
	output_header=1;
      }
  while (B)
    {
      output_lalign_aln   ( name, B);
      B=B->A;
    }
}
void output_lalign_aln   ( char *name, Alignment *B)
    {
    int a, b, c,d=0, s=0;
    char col;

    float tot=0;
    float id=0;

    FILE *fp;
    int max_len=0;
    int line;
    int *n_residues;
    int res;


    n_residues=(int*)vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    n_residues[a]=B->order[a][1];
	    }
    max_len=MAX(max_len+2, 16);
    line=60;



    fp= vfopen ( name, "a");

    for (a=0; a< B->len_aln; a++)
      {
        if ( !is_gap(B->seq_al[0][a]) && !is_gap(B->seq_al[1][a]))
	     {
	       tot++;
	       id+=(B->seq_al[0][a]==B->seq_al[1][a]);
	     }
      }

    id=(id*100)/tot;
    fprintf (fp, " %.1f%% identity in %d aa overlap; score: %d\n\n", id,(int)tot, B->score_aln);


    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<5; b++)
	     {
	         if ( b==0 || b==4)
		   {
		     if ( b==0)s=0;
		     if ( b==4)s=1;
		     fprintf (fp,"%-*s",max_len," ");
		     for (d=0,c=a;c<a+line && c<B->len_aln;c++)
		       {
			 res=!is_gap ( B->seq_al[s][c]);
			 n_residues[s]+=res;
			 if ( (n_residues[s]%10)==0 && res && (c-a+4)<line){fprintf (fp, "%-4d", n_residues[s]);d=-3;}
			 else
			   {
			     if ( d==0)fprintf (fp, " ");
			     else d++;
			   }
		       }
		     fprintf (fp,"\n");
		   }
		 else if (b==1 || b==3)
		    {
		      if ( b==1)s=0;
		      if ( b==3)s=1;
			fprintf (fp,"%-*s",max_len, B->name[s]);
			for (c=a;c<a+line && c<B->len_aln;c++)
			    {
				if ( is_gap(B->seq_al[s][c]))fprintf (fp,"%c", '-');
				else
				    {
					fprintf (fp, "%c", GET_CASE(B->residue_case, B->seq_al[s][c]));
				    }
			    }
			fprintf (fp,"\n");
		    }
		 else if ( b==2)
		    {
		    fprintf (fp,"%-*s",max_len," ");
		    for (c=a;c<a+line && c<B->len_aln;c++)
			    {
			    col=analyse_aln_column (B, c);
			    if ( col=='*')col=':';
			    else if ( col==':')col='.';
			    else if ( col=='.')col=' ';
			    fprintf ( fp, "%c", col);
			    }
		    fprintf (fp,"\n");
		    }
	     }
	   fprintf (fp,"\n");
	   }

    fprintf (fp,"\n\n----------\n\n");
    vfree (n_residues);
    vfclose ( fp);
    }


/****************************************************************************************************/
/*************************************UTIL *********************************************************/
/**************************************************************************************************/


/****************************************************************************************************/
/***************************                                    *************************************/
/***************************             PROCESSING 		*************************************/
/***************************                                    *************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              THREADING                                                  */
/***************************************************************************************** */

char *thread_aa_seq_on_dna_seq( char *s)
     {
	 int l, b, c;
	 char *array;


	 l=strlen ( s);
	 array=(char*)vcalloc ( l*3 +1, sizeof (char));
	 for ( b=0, c=0; b< l; b++, c+=3)
	     {
		 array[c]=s[b];
		 array[c+1]='o';
		 array[c+2]='o';
	     }
	 array[c]='\0';
	 return array;
     }

Alignment *thread_dnaseq_on_prot_aln (Sequence *S, Alignment *A)
        {
	    Alignment *B=NULL;
	    int a, b, c, n, la, ls, ln, m;

	    B=copy_aln ( A, B);
	    B=realloc_aln2 ( B, B->nseq, B->len_aln*3 +1);

	    for ( n=0,a=0; a< A->nseq; a++)
	        {
		for ( m=0,b=0; b< S->nseq; b++)
		    {
		    if (strm (A->name[a], S->name[b]) )
		       {
			   m=1;
			   n++;
			   ungap ( S->seq[b]);
			   B->seq_al[a][0]='\0';
			   for (la=0, ls=0, ln=0; la< A->len_aln; la++)
			       {
				 if (is_gap(A->seq_al[a][la]))
				   for (c=0; c< 3; c++)B->seq_al[a][ls++]='-';
				 else if (isupper(A->seq_al[a][la]))
				   for (c=0; c< 3; c++)B->seq_al[a][ls++]=toupper((S->seq[b][ln++]));
				 else 
				   for (c=0; c< 3; c++)B->seq_al[a][ls++]=tolower((S->seq[b][ln++]));
			       }
			   if (ln!=S->len[b])
			     {
			       printf_exit (EXIT_FAILURE, stderr, "thread_dna_seq_on_prot_aln: %s DNA/protein sequences differ [FATAL]", S->name[b]);
			     }
			   B->seq_al[a][ls]='\0';
		       }
		    }
		if ( m==0)
		       {
		       for (la=0, ls=0, ln=0; la< A->len_aln; la++)
			       {

				   B->seq_al[a][ls++]=A->seq_al[a][la];
				   B->seq_al[a][ls++]='-';
				   B->seq_al[a][ls++]='-';
			       }
		       }
		}

	    B->len_aln=strlen ( B->seq_al[0]);
	    return B;
	}
void thread_seq_struc2aln ( Alignment *A, Sequence *ST)
	{
	int a, b, c,d;
	int len, cons;

	for ( a=0; a< A->nseq; a++)
		for ( b=0; b< ST->nseq; b++)
			{
			if ( strcmp ( A->name[a], ST->name[b])==0)
				{
				ungap (ST->seq[b]);
				len=strlen(A->seq_al[a]);
				for ( c=0, d=0; c<len; c++)
					{
					if ( !is_gap(A->seq_al[a][c]))A->seq_al[a][c]=ST->seq[b][d++];
					}
				}
			}
	
	cons=name_is_in_list ("Cons", ST->name, ST->nseq, 100);
	if ( cons!=-1 && A->len_aln==strlen ( ST->seq[cons]))
	  {
	    sprintf (A->name[A->nseq], "Cons");
	    sprintf (A->seq_al[A->nseq],"%s", ST->seq[cons]);
	    A->nseq++;
	  }
	}
void cache_id ( Alignment *A)
	{
	int a, b,n;
	char r1, r2, r3;

	for ( a=0; a< A->len_aln; a++)
		{
		for ( b=0, n=0; b< A->nseq; b++)if ( !is_gap(A->seq_al[b][a]))n++;
		for ( b=0; b< A->nseq; b++)
			if ( !is_gap(A->seq_al[b][a]) && n==A->nseq)A->seq_al[b][a]='h';
			else if( !is_gap(A->seq_al[b][a]))A->seq_al[b][a]='x';
		}
	for ( a=0; a< A->nseq; a++)
		{
		for ( b=1; b< A->len_aln-1; b++)
			{
			r1=A->seq_al[a][b-1];
			r2=A->seq_al[a][b];
			r3=A->seq_al[a][b+1];
			if (r2=='h')
				{
				if ( (r1=='h' || r1=='b') && (r3=='h' || r3=='b'))A->seq_al[a][b]='h';
				else A->seq_al[a][b]='b';
				}
			}
		for ( b=1; b< A->len_aln-1; b++)if ( A->seq_al[a][b]=='b')A->seq_al[a][b]='x';
		}

	}


/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               PROCESING OF EST                                          */
/*                                                                                         */
/***************************************************************************************** */
int process_est_sequence ( Sequence *S, int *cluster_list)
	{
	char **inverted_seq;
	int T=20;
	int a, b;
	int V1, V2;
	int **sens;
	int **a_sens;
	int **best;
	int *solution;
	char buf [VERY_LONG_STRING];
	int n_clusters=0;
	int n;

	sens=declare_int ( S->nseq,S->nseq);
	a_sens=declare_int ( S->nseq,S->nseq);
	best=declare_int ( S->nseq,S->nseq);


	inverted_seq=(char**)vcalloc ( S->nseq, sizeof (char*));
	for ( a=0; a<S->nseq; a++)
		inverted_seq[a]=invert_seq ( S->seq[a]);

	for ( a=0; a< S->nseq-1; a++)
		{

		for ( b=a+1; b<S->nseq; b++)
			             {

			             V1=sens[a][b]=sens[b][a]=get_best_match ( S->seq[a], S->seq[b]);
			             V2=a_sens[a][b]=a_sens[b][a]=get_best_match ( S->seq[a],inverted_seq[b]);
				     best[a][b]=best[b][a]=(V1>V2)?V1:V2;
				     }
		}
	solution=SHC ( S->nseq, a_sens, sens);


	for ( a=0; a<S->nseq; a++)cluster_list[a]=-1;
	for ( a=0; a<S->nseq; a++)
		{
		n=search_for_cluster (a, n_clusters, cluster_list, T, S->nseq, best);
		if ( n>0)n_clusters++;
		}
	fprintf ( stderr, "\nTHERE %s %d Independant Cluster(s) in your sequences",(n_clusters>1)?"are":"is",(n_clusters));
	for (a=0; a<n_clusters; a++)
		{
		fprintf (stderr, "\n");
		for ( b=0; b<S->nseq; b++)
			{
			if ( cluster_list[b]==a)fprintf ( stderr, "%s ", S->name[b]);
			}
		}

	for ( a=0; a<S->nseq; a++)
		{
		if ( solution[a]==-1)
			{
			S->seq[a]=inverted_seq[a];
			sprintf ( buf, "i_%s", S->name[a]);
			sprintf ( S->name[a], "%s", buf);
			}
		}
	return n_clusters;
	}

int search_for_cluster ( int seq, int cluster_number, int *cluster_list, int T, int nseq, int **S)
	{
	int n=0,a;

	if (cluster_list[seq]==-1)
		{
		cluster_list[seq]=cluster_number;
		n++;
		}
	for ( a=0; a<nseq; a++)
		if ( cluster_list[a]==-1)
			{

			if (S[seq][a]>T)
				{
				n++;
				cluster_list[a]=cluster_number;
				n+=search_for_cluster ( a, cluster_number, cluster_list, T, nseq, S);
				}
			}
	return n;
	}

int * SHC ( int nseq, int **NST, int **ST)
	{
	int a;
	int mut;
	int score, new_score;
	int N_IT=VERY_LONG_STRING;
	int *sol;
	int count;

	sol=(int*)vcalloc ( nseq, sizeof (int));
	for ( a=0; a<nseq; a++)
		sol[a]=(addrand ((unsigned long)100)>49)?1:-1;

	score=evaluate_sol (sol, nseq, ST, NST);
	fprintf ( stderr, "\nI_Score=%d\n", score);
	N_IT=N_IT*nseq;

	for ( count=0,a=0; a< N_IT && score<VERY_LONG_STRING; a++, count++)
		{
		mut=mutate_sol ( sol,nseq);
		new_score=evaluate_sol (sol, nseq, ST, NST);
		if ( new_score>score)
			{
			score=new_score;
			}
		else if ( (addrand ((unsigned long)VERY_LONG_STRING))>score)
			{
			score=new_score;
			}
		else
			sol[mut]=sol[mut]*-1;
		if ( count==VERY_LONG_STRING)
			{
			count=0;
			fprintf ( stderr, "\nScore=%d", score);
			}
		}
	fprintf ( stderr, "\nScore=%d\n", score);
	return sol;
	}

int mutate_sol (int *sol, int nseq)
	{
	int n;
	n=addrand ((unsigned long)nseq);
	sol[n]=sol[n]*-1;
	return n;
	}
int evaluate_sol ( int *sol, int nseq, int **ST, int **NST)
	{
	static int max_score;
	int a, b, score=0;

	if ( max_score==0)
		{
		for ( a=0; a<nseq-1; a++)
			for ( b=a+1; b<nseq; b++)
				{
				max_score+=(ST[a][b]>NST[a][b])?ST[a][b]:NST[a][b];
				}
		}

	for ( a=0; a<nseq-1; a++)
		for (b=a+1; b<nseq; b++)
			if ( (sol[a]*sol[b])<0)score+=NST[a][b];
			else score+=ST[a][b];
	return (score*VERY_LONG_STRING)/max_score;
	}


char * invert_seq ( char *seq)
	{
	int a, b;

	char *nseq;
	int l;


	l=strlen ( seq);
	for ( a=0; a<l; a++)
		seq[a]=tolower ( seq[a]);
	nseq=(char*)vcalloc ( l+1, sizeof (char));

	for ( a=0, b=l-1; a<l; a++, b--)
		{
		if (seq[b]=='n')nseq[a]='n';
		else if (seq[b]=='g')nseq[a]='c';
		else if (seq[b]=='c')nseq[a]='g';
		else if (seq[b]=='a')nseq[a]='t';
		else if (seq[b]=='t')nseq[a]='a';
		}

	nseq[l]='\0';
	return nseq;
	}


int get_best_match ( char *seq1, char *seq2)
	{
	static int **m;
	static int ml;
	int a, b;
	int **mdiag;
	int n_mdiag=0;
	int best;
	int l1, l2;


	l1=strlen ( seq1);
	l2=strlen (seq2);
	if ( m==NULL)
		{
		ml=(l1>l2)?l1:l2;
		m=declare_int (ml, ml);
		}
	else if ( (ml<l1) || (ml<l2))
		{
		free_int (m, ml);
		ml=(l1>l2)?l1:l2;
		m=declare_int (ml, ml);
		}

	for ( a=0; a<l1; a++)
		{
		for ( b=0; b<l2; b++)
			m[a][b]=((seq1[a]==seq2[b])|| seq1[a]=='n' ||seq2[b]=='n')?1:0;
		}
	mdiag= extract_m_diag_streches ( m, l1, l2,seq1, seq2, &n_mdiag);

	for ( best=0,a=0; a<n_mdiag; a++)
		best=(mdiag[a][0]>best)?mdiag[a][0]:best;

	return best;
	}

int** extract_m_diag_streches ( int ** m, int l1, int l2,char *seq1, char *seq2, int *n_mdiag)
	{

	int b, x, y, s1, s2;
	static int **mdiag;
	int in;
	static int max_diag=VERY_LONG_STRING;

	 /*
	 diag[0]=len;
	 diag[1]=x_start;
	 diag[2]=y_start;
	 diag[3]=x_end;
	 diag[4]=y_end;
	 */

	if ( mdiag==NULL)
		mdiag=declare_int ( max_diag, 5);

	for ( s1=l1-1, s2=0;s2<l2;)
		{
		for ( in=0,x=s1, y=s2; x<l1 && y<l2; x++, y++)
			{
			if (m[x][y]>0)
				{
				if (in==1)
					mdiag[n_mdiag[0]][0]++;
				else
					{
					mdiag[n_mdiag[0]][0]=1;
					mdiag[n_mdiag[0]][1]=x;
					mdiag[n_mdiag[0]][2]=y;
					in=1;
					}
				}
			else
				if (in==1)
					{
					in=0;
					mdiag[n_mdiag[0]][3]=x-1;
					mdiag[n_mdiag[0]][4]=y-1;
					if ( !is_strech ( "ta", seq1, seq2,mdiag[n_mdiag[0]][0], mdiag[n_mdiag[0]][1],mdiag[n_mdiag[0]][2]))n_mdiag[0]++;
					}
			if (n_mdiag[0]==(max_diag-1))
				{mdiag=(int**)vrealloc (mdiag, (max_diag+VERY_LONG_STRING)*sizeof (int*));
				for ( b=max_diag; b<max_diag+VERY_LONG_STRING; b++)mdiag[b]=(int*)vcalloc ( 5, sizeof (int));
				max_diag+=VERY_LONG_STRING;
				}
			}
		s2+= (s1==0)?1:0;
		s1-= (s1==0)?0:1;
		if (in==1)
			{
			in=0;
			mdiag[n_mdiag[0]][3]=x-1;
			mdiag[n_mdiag[0]][4]=y-1;
			if ( !is_strech ( "ta", seq1, seq2,mdiag[n_mdiag[0]][0], mdiag[n_mdiag[0]][1],mdiag[n_mdiag[0]][2]))n_mdiag[0]++;
			}
		}

	return mdiag;
	}
int is_strech ( char *AA, char *seq1, char *seq2, int len, int x, int y)
	{
	int n, i, j, c,a,nr;
	int T=70;

	n=strlen ( AA);
	for ( a=0; a<n; a++)
		{
		for (nr=0, i=x, j=y, c=0; c<len; c++, i++, j++)
			if ((seq1[i]==AA[a]) && (seq2[j]==AA[a]))nr++;
		if ( ((nr*100)/len)>T)return 1;
		}
	return 0;
	}


/************************************************************************************/
/*                                                                                  */
/*                                      STRUC                                       */
/*                                                                                  */
/*                                                                                  */
/************************************************************************************/

char * oneletaa2threeletaa(char aa);
float aa2property   (char aa, char *mode);

int output_seq2struc(char *outfile, Alignment *A)
{
  FILE *fp1, *fp2;
  int a,c, l;
  float v, h, x, y, z, dx, dy, dz;
  char *s;
  char *tmpfile1, *tmpfile2;
  char command[1000];

  tmpfile1=vtmpnam(NULL);
  tmpfile2=vtmpnam(NULL);

  ungap (A->seq_al[0]);
  s=A->seq_al[0];l=strlen (s);
  fp1=vfopen (tmpfile1, "w");

  x=y=z=0;
  for ( a=0; a< l; a++)
    {
      h=aa2property ( s[a], "doolittle"   );
      v=aa2property (s[a], "volume");
      /*14.398907: peptide bond length*/
      dx=(float)sqrt ((double)(14.398907/(((h*h)/(v*v))+1)));
      dy=dx*(h/v);
      dz=0;


      x+=dx;
      y+=dy;
      z+=dz;
      fprintf (fp1, "ATOM%7d   CA %s A%4d%12.3f%8.3f%8.3f  1.00   5.30\n",a+1, oneletaa2threeletaa(s[a]),a+1, x, y, z);
    }
  vfclose (fp1);
  sprintf ( command, "extract_from_pdb -infile %s -force > %s", tmpfile1,  tmpfile2);
  my_system  (command);
  fp1=vfopen (tmpfile2, "r");
  fp2=vfopen (outfile, "w");

  while ( (c=fgetc(fp1))!=EOF)fprintf (fp2, "%c", c);
  vfclose (fp1);
  vfclose (fp2);

  return 0;
}

char * oneletaa2threeletaa(char aa)
  {
    aa=tolower (aa);
    if ( aa=='a')return "ALA";
    else if ( aa=='r') return "ARG";
    else if ( aa=='n') return "ASN";
    else if ( aa=='d') return "ASP";
    else if ( aa=='c') return "CYS";
    else if ( aa=='q') return "GLN";
    else if ( aa=='e') return "GLU";
    else if ( aa=='g') return "GLY";
    else if ( aa=='h') return "HIS";
    else if ( aa=='i') return "ILE";
    else if ( aa=='l') return "LEU";
    else if ( aa=='k') return "LYS";
    else if ( aa=='m') return "MET";
    else if ( aa=='f') return "PHE";
    else if ( aa=='p') return "PRO";
    else if ( aa=='s') return "SER";
    else if ( aa=='t') return "THR";
    else if ( aa=='w') return "TRP";
    else if ( aa=='y') return "TYR";
    else if ( aa=='v') return "VAL";
    else
      {
	fprintf ( stderr, "\nERROR: %c is not an amino acid [FATAL::aa2hydropathy::%s]", aa, PROGRAM);
	myexit (EXIT_FAILURE);
	return NULL;
      }
    return NULL;
  }

float aa2property   (char aa, char *mode)
  {
    if ( mode==NULL || strm (mode, "doolittle"))
	 {
	   aa=tolower (aa);
	   if ( aa=='i')return 4.5;
	   else if ( aa=='v') return 4.2;
	   else if ( aa=='l') return 3.8;
	   else if ( aa=='f') return 2.8;
	   else if ( aa=='c') return 2.5;
	   else if ( aa=='m') return 1.9;
	   else if ( aa=='a') return 1.8;
	   else if ( aa=='g') return -0.4;
	   else if ( aa=='t') return -0.7;
	   else if ( aa=='w') return -0.9;
	   else if ( aa=='s') return -0.8;
	   else if ( aa=='y') return -1.3;
	   else if ( aa=='p') return -1.6;
	   else if ( aa=='h') return -3.2;
	   else if ( aa=='e') return -3.5;
	   else if ( aa=='q') return -3.5;
	   else if ( aa=='d') return -3.5;
	   else if ( aa=='n') return -3.5;
	   else if ( aa=='k') return -3.9;
	   else if ( aa=='r') return -4.5;
	   else
	     {
	       fprintf ( stderr, "\nERROR: %c is not an amino acid [FATAL::aa2hydropathy::%s]", aa, PROGRAM);
	       myexit (EXIT_FAILURE);
	     }
	 }
    else if (strm (mode, "volume"))
	 {
	   aa=tolower (aa);
	   if ( aa=='a')return 0.915;
	   else if ( aa=='r') return 2.02;
	   else if ( aa=='n') return 1.35;
	   else if ( aa=='d') return 1.24;
	   else if ( aa=='c') return 1.18;
	   else if ( aa=='q') return 1.61;
	   else if ( aa=='e') return 1.55;
	   else if ( aa=='g') return 0.66;
	   else if ( aa=='h') return 1.67;
	   else if ( aa=='i') return 1.69;
	   else if ( aa=='l') return 1.68;
	   else if ( aa=='k') return 1.71;
	   else if ( aa=='m') return 1.70;
	   else if ( aa=='f') return 2.03;
	   else if ( aa=='p') return 1.29;
	   else if ( aa=='s') return 0.99;
	   else if ( aa=='t') return 1.22;
	   else if ( aa=='w') return 2.37;
	   else if ( aa=='y') return 2.03;
	   else if ( aa=='v') return 1.41;
	   else
	     {
	       fprintf ( stderr, "\nERROR: %c is not an amino acid [FATAL::aa2hydropathy::%s]", aa, PROGRAM);
	       myexit (EXIT_FAILURE);
	     }
	 }

    else
      {
	fprintf ( stderr, "\nERROR: %s is an unknown mode [FATAL::aa2hydropathy::%s]", mode  , PROGRAM);
	myexit (EXIT_FAILURE);
      }
  return 0;
  }





/************************************************************************************/
/*                                                                                  */
/*                                      DNA                                         */
/*                                                                                  */
/*                                                                                  */
/************************************************************************************/

Alignment *code_dna_aln (Alignment *A)
       {
	 int a, b,l,r;

	 for ( a=0; a< A->nseq; a++)
	   {
	     for (l=0, b=0; b< A->len_aln; b++)
	       {
		 r=A->seq_al[a][b];
		 if ( r=='-')l++;
		 else if ( r=='~')continue;
		 else if ( r=='.')l++;
		 else if ( !islower(r))A->seq_al[a][b]='4';
		 else
		   {
		     A->seq_al[a][b]=(l+3)%3+'0';
		     l++;
		   }
	       }
	   }
	 return A;
       }


Alignment *back_translate_dna_aln (Alignment *A)
       {
	 /*Given a set of aligned sequences
	   starts from left to right
	   1 aa->3 nuc
	   ambiguities are randomly resolved.
	   returns the corresponding amino acid alignment
	 */
	  int a;
	  char *seq    ;

	 ungap_aln(A);
	 A=realloc_aln (A, 10000);
	 seq=(char*)vcalloc ( 10000, sizeof (char));


	 for ( a=0; a< A->nseq; a++)
	   {
	   seq=back_translate_dna_seq (A->seq_al[a], seq, RANDOM);
	   sprintf ( A->seq_al[a], "%s", seq);
	   }
	 A->len_aln=A->len_aln*3;
	 compress_aln (A);
	 vfree (seq);
	 return A;
       }
char * back_translate_dna_seq ( char *in_seq,char *out_seq, int mode)
       {
	 int a,len;

	 len=strlen(in_seq);

	 if (out_seq==NULL)out_seq=(char*)vcalloc ( len*3+1, sizeof (char));

	 out_seq[0]='\0';
	 for (a=0; a<len; a++)
	   {
	   strcat (out_seq,  back_translate_dna_codon (in_seq[a],mode));
	   }

	 return out_seq;
       }

static Sequence *rna_seq2dna_seq (Sequence *S);
static Sequence *dna_seq2rna_seq (Sequence *S);

/**
 * Transforms RNA to DNA or the other way around.
 *
 * Will invoke ::rna_seq2dna_seq or ::dna_seq2rna_seq or give an error, if mode is unknown.
 * \param[in,out] S Points to a Sequence object.
 * \param mode Should be either \b rna2dna or \b dna2rna
 */
Sequence * transform_sequence ( Sequence *S, char *mode)
{
  if ( strm (mode, "rna2dna"))
    return rna_seq2dna_seq (S);
  else if ( strm (mode, "dna2rna"))
    return dna_seq2rna_seq (S);
  else
    printf_exit (EXIT_FAILURE, stderr, "Unknown -transform mode: %s [FATAL:%s]\n", mode,PROGRAM);
  return NULL;}

Sequence *rna_seq2dna_seq (Sequence *S)
{
  int a, b;

  if ( !strm(S->type, "DNA") && !strm (S->type, "RNA")) printf_exit (EXIT_FAILURE, stderr, "Sequences should be *RNA* type [FATAL:%s]\n", PROGRAM);
  for ( a=0; a<S->nseq; a++)
    {
      for (b=0; b<strlen (S->seq[a]); b++)
	{
	  if ( S->seq[a][b]=='u') S->seq[a][b]='t';
	  if ( S->seq[a][b]=='U') S->seq[a][b]='T';
	}
    }
  return S;
}
Sequence *dna_seq2rna_seq (Sequence *S)
{
  int a, b;

  if ( !strm(S->type, "DNA") && !strm (S->type, "RNA")) printf_exit (EXIT_FAILURE, stderr, "Sequences should be *DNA* type (type=%s) [FATAL:%s]\n", PROGRAM, S->type);
  for ( a=0; a<S->nseq; a++)
    for (b=0; b<S->len[a]; b++)
      {
	if ( S->seq[a][b]=='t') S->seq[a][b]='u';
	if ( S->seq[a][b]=='T') S->seq[a][b]='U';
      }
  return S;
}



int get_longest_frame (char *seq, int mode);
Sequence  *translate_dna_seqS     (Sequence *S, int frame, int stop)
{
  //frame: 1->3
  char *s;
  int a, b,c,l;

  for (a=0; a<S->nseq; a++)
    {
      s=S->seq[a];
      ungap(s);
      l=strlen (s);
      for (b=(frame-1); b<l; b+=3)
	{
	  s[b]=translate_dna_codon (s+b,stop);
	  for (c=b+1; c<b+3 && c<l; c++)s[c]='-';
	}
    }
  return S;
}
Alignment *translate_dna_aln (Alignment *A, int frame)
       {
	 /*Given a set of aligned sequences
	   starts from left to right
	   3 nuc->1 aa
	   2nuc+1gap, 1nuc+2gap->3 gaps
	   1 stop-> 3gaps
	   returns the corresponding amino acid alignment
	 */


	 int a, b,r;


	 if (frame==3 || frame ==4)
	   {

	     for (a=0; a< A->nseq; a++)
	       {
		 char *d, *buf, f;
		 d=A->seq_al[a];
		 f=get_longest_frame (d,frame);
		 buf=(char*)vcalloc ( strlen (d)+1, sizeof (char));
		 if ( f<3)
		   {
		     sprintf (buf, "%s", d+f);
		     sprintf (d, "%s", buf);
		     sprintf (A->seq_comment[a], " frame: %d", f);
		   }
		 else if ( f>=3)
		   {
		     f-=3;
		     sprintf ( buf, "%s", d);
		     buf=complement_string (buf);
		     sprintf (d, "%s",buf+f);
		     sprintf (A->seq_comment[a], " frame: %d Reverse Complement", f);
		   }
		 vfree (buf);
	       }
	   }
	 else
	   {

	     for ( a=0; a< A->nseq; a++)
	       for (b=0; b< frame; b++)
		 A->seq_al[a][b]='-';
	     ungap_aln(A);
	   }

	 for ( b=0; b< A->nseq; b++)
	   {
	     for ( a=0; a< A->len_aln;)
	       {
		 
		 r=translate_dna_codon (A->seq_al[b]+a, 'z');
		 
		 if (is_gap(r))
		   {
		     A->seq_al[b][a++]='-';
		     A->seq_al[b][a++]='-';
		     A->seq_al[b][a++]='-';
		   }
		 else if ( r=='x' || r=='X')
		   {
		     A->seq_al[b][a++]=(r=='x')?'o':'O';
		     A->seq_al[b][a++]='-';
		     A->seq_al[b][a++]='-';
		   }
		 else if ( r=='z' || r=='Z')
		   {
		     A->seq_al[b][a++]=(r=='z')?'x':'X';
		     A->seq_al[b][a++]='-';
		     A->seq_al[b][a++]='-';
		   }
		 else
		   {
		     A->seq_al[b][a++]=r;
		     A->seq_al[b][a++]='-';
		     A->seq_al[b][a++]='-';
		   }
	       }
	   }
	 compress_aln (A);

	 return A;
       }
int seq2blastdb (char *out, Sequence *S)
{

  output_fasta_simple (out, S);

  if ( strm (S->type, "DNA") || strm (S->type, "RNA") )
	  printf_system ("makeblastdb -in %s -dbtype nucl -logfile /dev/null", out);
  else {
	  fprintf(stderr, "WARNING: makeblastdb is used with -dbtype prot\n");
	  printf_system ("makeblastdb -in %s -dbtype prot -logfile /dev/null", out);
  }
  return 1;
}
int seq2tblastx_db (char *out,Sequence *S, int strand)
{
  //strand : same values as in ncbi blastall
  //1: direct
  //2:revers
  //3: both
  int a, b,d, l;
  char *prot, *pprot;
  int min_exon_len=5;
  FILE *fp;
  fp=vfopen (out, "w");
  for (a=0; a<S->nseq; a++)
    {
      for (b=-3; b<=3; b++)
	{
	  int f;
	  int dl;

	  dl=strlen (S->seq[a]);
	  if (b==0)continue;
	  else if ( strand==1 && b<0)continue;//only direct strand
	  else if ( strand==2 && b>0)continue;//only reverse strand
	  else if (b<0)
	    {
	      S->seq[a]=complement_string (S->seq[a]);
	      f=b*-1;
	    }
	  else
	    f=b;
	  prot=translate_dna_seq (S->seq[a], f, 'X', NULL);
	  upper_string (prot);

	  l=strlen (prot);
	  pprot=prot;
	  for (pprot=prot,d=0; d<=l; d++)
	    {
	      if (prot[d]=='\0')
		{
		  prot[d]='\0';
		  if ( strlen (pprot)>min_exon_len)
		    {
		      int start, end;
		      end=(d)*3+(f-1);
		      start=(end-(strlen (pprot))*3)+1;
		      fprintf (fp, ">%s__%c__%d__%d__%d\n%s\n", S->name[a],(b>0)?'d':'r',start,end,dl, pprot);
		      pprot=prot+d+1;
		    }
		}
	    }
	  vfree (prot);
	  if (b<0) S->seq[a]=complement_string (S->seq[a]);
	}
    }
  vfclose (fp);
  return EXIT_SUCCESS;
}

int get_longest_frame (char *in_seq, int mode)
{
  char *prot, *seq;
  int a;
  int max_l=0, l;
  int best_frame=0;
  int nf;

  seq=(char*)vcalloc (strlen (in_seq)+1, sizeof (char));
  prot=(char*)vcalloc (strlen (in_seq)+1, sizeof (char));
  sprintf ( seq, "%s", in_seq);

  if ( mode == 3)nf=3;
  else if ( mode == 4) nf=6;

  for (a=0; a<nf; a++)
    {
      int f;
      if (a==3)seq=complement_string (seq);
      f=(a>=3)?a-3:a;
      prot=translate_dna_seq ( seq,f,'\0', prot);
      l=strlen (prot);
      if (l>=max_l){max_l=l;best_frame=a;}
    }
  vfree (seq);
  vfree (prot);
  return best_frame;
}

Alignment *clean_gdna_aln (Alignment *A)
       {
	   int a, b, c, r1, r2,s, p, n, tn;
	   int *col;
	   static int **mat;
	   Alignment *T=NULL;
	   int **score;
	   char *buffer;


	   /*Viterbi Parameters*/
	   int AL=0;        /*Allowed Transition*/
	   int F=-1000000; /*Forbiden Transition*/
	   int SPLICE_PENALTY=100;
	   int ORF1=0, ORF2=1, ORF3=2, NC=3;

	   int state, pstate, best_e, best_pstate_p,best_state_p, best_pstate_v, best_state_v, v;
	   int nstate=4;
	   int **transitions;
	   int e;
	   int **v_tab_p;
	   int **v_tab;
	   int * is_dna;

	   best_state_p=best_state_v=best_pstate_p=best_pstate_v=best_e=0;
	   buffer=(char*)vcalloc ( 100000, sizeof (char));
	   is_dna=(int*)vcalloc ( A->nseq, sizeof (int));
	   score=declare_int ( A->nseq+1, A->len_aln);


	   if ( !mat)mat=read_matrice("pam250mt");
	   T=copy_aln (A, T);
	   col=(int*)vcalloc ( A->nseq, sizeof (int));

	   for (a=0; a<= A->len_aln; a++)
	       for ( b=0; b< A->nseq; b++){A->seq_al[b][a]=tolower(A->seq_al[b][a]); A->seq_al[b][a]=(A->seq_al[b][a]=='t')?'u':A->seq_al[b][a];}

	   for ( a=0; a< A->nseq; a++)
	       {
		   sprintf ( buffer, "%s", A->seq_al[a]);
		   ungap (buffer);
		   is_dna[a]=strm ( get_string_type (buffer), "DNA");
	       }


           for (a=0; a< A->len_aln-2; a++)
	       {
	       for (b=0; b< A->nseq; b++)
		       {
		       if (is_dna[b])col[b]=translate_dna_codon (A->seq_al[b]+a, 'x');
		       else col[b]=tolower ( A->seq_al[b][a]);
		       }

	       for (n=0,tn=0,b=0; b< A->nseq; b++)
		   for ( c=b; c< A->nseq; c++   )
		       {
			   r1=col[b];
			   r2=col[c];

			   if (r1=='x' || r2=='x'){score[A->nseq][a]=F;break;}
			   else if (r1=='-' && r2=='-');
			   else if (r1=='-' || r2=='-');
			   else
			       {

				   if ( is_dna[b] && is_dna[c])score[A->nseq][a]+= mat[r1-'A'][r2-'A'];
				   else score[A->nseq][a]+=mat[r1-'A'][r2-'A']* (A->nseq*A->nseq);
			       }
			   n+=( !is_gap(r1) && !is_gap(r2));
			   score[A->nseq][a]=(((tn!=0)?score[A->nseq][a]/tn:0));
		       }

	       }

	   /*initialisation*/

	   transitions=declare_int ( nstate, nstate);
	   v_tab=declare_int ( A->len_aln+2, nstate       );
	   v_tab_p=declare_int ( A->len_aln+2, nstate       );

	   for (a=0; a<nstate;a++)
	       for (b=0; b<nstate;b++)
	             {transitions[a][b]=F;}

	   transitions[ORF1][ORF2]=AL;
	   transitions[ORF2][ORF3]=AL;
	   transitions[ORF3][ORF1]=AL;

	   transitions[ORF3][NC]  =AL-SPLICE_PENALTY;
	   transitions[NC][ORF1]  =AL-SPLICE_PENALTY;


	   for ( s=0; s<A->nseq; s++)
	       {
	       for ( p=0; p<=A->len_aln; p++){for (state=0; state< nstate; state++)v_tab_p[p][state]=-1; }
	       for (p=1+2; p<= A->len_aln; p++)
	           {

		   for (state=0; state< nstate; state++)
		       {

			   if ( state==NC){e=-best_e;}
			   else
			      {
				  e=score[A->nseq][(p-1)-state];
				  if ( state==0)best_e=e;
				  else best_e=MAX(e, best_e);
			      }

			   for ( pstate=0; pstate<nstate; pstate++)
		               {
				   v=e+transitions[pstate][state]+v_tab[p-1][pstate];
				   if (pstate==0 ||(v>best_pstate_v) )
				      {
				       best_pstate_v=v;
				       best_pstate_p=pstate;
				      }
			       }

			   v_tab[p][state]=best_pstate_v;
			   v_tab_p[p][state]=best_pstate_p;
			   if (state==0 ||best_pstate_v>best_state_v )
			      {
			       best_state_p=state;
			       best_state_v=best_pstate_v;
			      }
		       }

		   }



	       for (p=0; p< A->len_aln; p++)T->seq_al[s][p]='.';
	       for (p=A->len_aln; p>0; p--)
	           {

		       if ( best_state_p==0)T->seq_al[s][p-1]=translate_dna_codon (A->seq_al[s]+(p-1), 'x');
		       else if ( best_state_p==1 || best_state_p==2)T->seq_al[s][p-1]='-';



		       best_state_p=v_tab_p[p][best_state_p];

		   }
	       }



	   vfree (col);
	   return T;
       }

Alignment *clean_cdna_aln (Alignment *A)
       {
	 /*Given an alignmnet of nucleotides
	   Returns the same alignmnent whith non coding nucleotides replaced with dots

	   at each position, the emission probability is the sum of pair of the substitution of amino-acids
	 */

	   int a, b, c,s, p;
	   static int **mat;
	   int   *emission;
	   float em1, em2;
	   char *buffer;
	   Alignment *B=NULL;




	   /*Viterbi Parameters*/
	   int AL=0;        /*Allowed Transition*/
	   int F=-1000000; /*Forbiden Transition*/
	   int PENALTY=30;
	   int NC, C1,C2, C3, START, END;
	   int nstate=0;
	   int state=0,best_state=0, score=0, best_score=0;
	   int p_state;
	   int e=0;
	   int **score_tab;
	   int **state_tab;

	   int **transitions;
	   int n;
	   int r1, r2, r3;

	   NC=nstate++;
	   C1=nstate++;
	   C2=nstate++;
	   C3=nstate++;
	   START=nstate++;
	   END=nstate++;


	   B=copy_aln (A, B);
	   buffer=(char*)vcalloc ( 100000, sizeof (char));
	   emission=(int*)vcalloc (A->len_aln, sizeof (int));

	   if ( !mat)
	     {
	       mat=read_matrice("pam250mt");
	     }

	   /*Computation of the emission proba for the coding state*/


	   for (a=0; a< A->len_aln; a++)
	     {

	       /*First component: % occupancy of the column*/
	       em1=0;
	       for ( b=0; b< A->nseq; b++) em1+=!is_gap(translate_dna_codon (A->seq_al[b]+a, '-'));
	       em1=em1/(float)A->nseq;

	       /*Second Component: % similarity within column*/
	       em2=0;
	       for (n=0,b=0; b< A->nseq-1; b++)
		 {
		   r1=translate_dna_codon (A->seq_al[b]+a, '-');

		   for (c=b+1; c<A->nseq; c++)
		     {
		       r2=translate_dna_codon (A->seq_al[c]+a, '-');
		       if (is_gap(r2) || is_gap(r1));
		       else
			 {
			   n++;
			   em2+=((mat[r1-'A'][r2-'A'])>1)?1:0;
			 }
		     }
		 }
	       em2=em2/(float)((n==0)?1:n);


	       emission[a]=(em1*100);

	     }



	   /*initialisation*/

	   transitions=declare_int ( nstate, nstate);
	   score_tab=declare_int ( A->len_aln+2, nstate       );
	   state_tab=declare_int ( A->len_aln+2, nstate       );

	   for (a=0; a<nstate;a++)
	       for (b=0; b<nstate;b++)
	             {transitions[a][b]=F;}


	   transitions[START][C1]=AL;
	   transitions[START][NC]=AL;
	   transitions[C3][END]=AL;
	   transitions[NC][END]=AL;
	   transitions[C1 ][C2 ]=AL;
	   transitions[C2 ][C3 ]=AL;
	   transitions[C3 ][C1 ]=AL;
	   transitions[C3 ][NC ]=AL-PENALTY;
	   transitions[NC ][C1 ]=AL-PENALTY;
	   transitions[NC][NC]=AL-PENALTY;



	   for ( s=0; s< A->nseq; s++)
	     {
	     for ( p=0; p<=A->len_aln; p++){for (state=0; state< nstate; state++){score_tab[p][state]=F;state_tab[p][state]=-1;} }
	     score_tab[0][START]=0;

	     for (p=1; p<= A->len_aln; p++)
	       {
		 for (state=0; state< nstate; state++)
		   {
		     if ( state==START || state==END)continue;
		     else if      ( state==NC)  e=-10;
		     else if ( state==C1)
		       {
			 e=emission[p-1];
		       }
		     else if ( state ==C2)
		       {
			 if ( p-2<0)e=F;
			 else e=emission[p-2];
		       }
		     else if ( state==C3)
		       {
			 if ( p-3<0)e=F;
			 else e=emission[p-3];
		       }

		     for (p_state=0; p_state<nstate; p_state++)
		       {

			 if (e==F)score=F;
			 else
			   {
			     score=(score_tab[p-1][p_state]==F)?F:(e+transitions[p_state][state]+score_tab[p-1][p_state]);
			   }

			 if(p_state==0 || score>best_score){ best_score=score;best_state=p_state;}

		       }

		     score_tab[p][state]=best_score;
		     state_tab[p][state]=best_state;

		   }
	       }

	     best_score=best_state=UNDEFINED;
	     for (state=0; state<nstate; state++)
	       {
		 if (state==START || state==END)continue;
		 e=transitions[state][END];
		 if (e==F || score_tab[p-1][state]==F)continue;

		 if (best_score==UNDEFINED || score_tab[p-1][state]>best_score)
		   {
		     best_score=score_tab[p-1][state]+e;
		     best_state=state;
		   }

	       }

	     for (p=A->len_aln; p>0;)
	       {
		 B->seq_al[s][p-1]=best_state+'0';
		 best_state=state_tab[p][best_state];
		 p--;
	       }
	     }

	   for ( a=0; a< A->nseq; a++)
	     for ( b=0; b< A->len_aln;)
	       {
		 s=B->seq_al[a][b];
		 if ( s==C1+'0')
		   {
		     r1=A->seq_al[a][b];
		     r2=A->seq_al[a][b+1];
		     r3=A->seq_al[a][b+2];


		     if ( is_gap(r1) ||is_gap(r2) ||  is_gap(r3))
		       {
			 A->seq_al[a][b]=(is_gap(r1))?'~':'.';
			 A->seq_al[a][b+1]=(is_gap(r2))?'~':'.';
			 A->seq_al[a][b+2]=(is_gap(r3))?'~':'.';
		       }
		     b+=3;
		   }
		 else if ( s==NC+'0')
		   {
		     A->seq_al[a][b]=(is_gap(A->seq_al[a][b]))?'~':'.';
		     b++;
		   }
		 else
		   {
		     fprintf (stderr, "\nPROBLEM: [%d %d]->%d", a, b, s-'0');
		   }
	       }


	   free_aln (B);
	   free_int (transitions, -1);
	   free_int (score_tab, -1);
	   free_int (state_tab, -1);
	   vfree (emission);
	   vfree (buffer);

	   return A;
       }




Alignment *translate_splice_dna_aln (Alignment *A, Alignment *ST)
       {
	   int a, b, c, r1, r2,s, p, n, tn;
	   int *col;
	   static int **mat;
	   Alignment *T=NULL;
	   int **score;

	   /*Viterbi Parameters*/
	   int AL=0;        /*Allowed Transition*/
	   int F=-1000000; /*Forbiden Transition*/
	   int ORF1=0, ORF2=1, ORF3=2,SPL1=3, SPL2=4, SPL3=5, SPL4=6, NC=7;
	   int SPLICE_PENALTY;
	   int frame1, frame2, frame3, best_frame;
	   int nstate=8;
	   char r;



	   int state=0, pstate=0, best_pstate_p=0,best_state_p=0, best_pstate_v=0, best_state_v=0, v=0;

	   int **transitions;
	   int e=0;
	   int **v_tab_p;
	   int **v_tab;

	   score=declare_int ( A->nseq+1, A->len_aln);


	   if ( !mat)mat=read_matrice("pam250mt");
	   T=copy_aln (A, T);
	   col=(int*)vcalloc ( A->nseq, sizeof (int));

	   for (a=0; a<= A->len_aln; a++)
	       for ( b=0; b< A->nseq; b++){A->seq_al[b][a]=tolower(A->seq_al[b][a]); A->seq_al[b][a]=(A->seq_al[b][a]=='t')?'u':A->seq_al[b][a];}




	   for (a=0; a< A->len_aln-2; a++)
	       {
	       for (b=0; b< A->nseq; b++)
		       {
		       col[b]=translate_dna_codon (A->seq_al[b]+a, 'x');
		       }

	       for (n=0,tn=0,b=0; b< A->nseq-1; b++)
		   for ( c=b+1; c< A->nseq; c++, tn++   )
		       {
			   r1=col[b];
			   r2=col[c];

			   if (r1=='x' || r2=='x')score[A->nseq][a]=F;
			   else if (r1=='-' && r2=='-');
			   else if (r1=='-' || r2=='-');
			   else
			       {
				   score[A->nseq][a]+= mat[r1-'A'][r2-'A'];

			       }
			   n+=( !is_gap(r1) && !is_gap(r2));
		       }
 	       score[A->nseq][a]=(((tn!=0)?score[A->nseq][a]/tn:0));

	       }

	   /*initialisation*/

	   transitions=declare_int ( nstate, nstate);
	   v_tab=declare_int ( A->len_aln+2, nstate*nstate);
	   v_tab_p=declare_int ( A->len_aln+2, nstate*nstate);

	   for (a=0; a<nstate;a++)
	     for (b=0; b<nstate;b++)
	       {transitions[a][b]=F;}

	   SPLICE_PENALTY=-1000;

	   transitions[ORF1][ORF2]    =AL;
	   transitions[ORF1][SPL1]    =AL-SPLICE_PENALTY;

	   transitions[ORF2][ORF3]    =AL;
	   transitions[ORF2][SPL1]    =AL-SPLICE_PENALTY;

	   transitions[ORF3][ORF1]    =AL;
	   transitions[ORF3][SPL1]    =AL-SPLICE_PENALTY;

	   transitions[ORF3][ORF1]    =AL;
	   transitions[ORF3][SPL1]    =AL-SPLICE_PENALTY;

	   transitions[ORF3][NC]=AL-100;
	   transitions[NC][ORF1]=AL-100;


	   transitions[SPL1][SPL2]=AL;
	   transitions[SPL2][NC  ]=AL-SPLICE_PENALTY;
	   transitions[NC  ][NC  ]=AL;
	   transitions[NC  ][SPL3]=AL-SPLICE_PENALTY;
	   transitions[SPL3][SPL4]=AL;
	   transitions[SPL4][ORF1]=AL;
	   transitions[SPL4][ORF2]=AL;
	   transitions[SPL4][ORF3]=AL;


	   for ( s=0; s<A->nseq; s++)
	       {
	       for ( p=0; p<=A->len_aln; p++){for (state=0; state< nstate; state++)v_tab_p[p][state]=-1; }
	       for (p=1+2; p<= A->len_aln; p++)
	           {
		    frame1=score[A->nseq][(p-1)];
		    frame2=score[A->nseq][(p-1)-1];
		    frame3=score[A->nseq][(p-1)-2];
		    best_frame=best_int (3, 1, &a, frame1, frame2, frame3);
		    for (state=0; state< nstate; state++)
		       {
			 r=tolower (A->seq_al[s][p-1]);
			 r=(r=='u')?'t':r;

			 if      (state==ORF1)e=frame1;
			 else if (state==ORF2)e=frame2;
			 else if (state==ORF3)e=frame3;
			 else if (state==SPL1)e=(r=='g')?best_frame:F;
			 else if (state==SPL2)e=(r=='t')?best_frame:F;
			 else if (state==SPL3)e=(r=='a')?best_frame:F;
			 else if (state==SPL4)e=(r=='g')?best_frame:F;
			 else if (state==NC)e=-best_frame;
			 for ( pstate=0; pstate<nstate; pstate++)
		               {
				   v=e+transitions[pstate][state]+v_tab[p-1][pstate];
				   if (pstate==0 ||(v>best_pstate_v) ){best_pstate_v=v;best_pstate_p=pstate;}
			       }

			   v_tab[p][state]=best_pstate_v;
			   v_tab_p[p][state]=best_pstate_p;
			   if (state==0 ||best_pstate_v>best_state_v ){best_state_p=state; best_state_v=best_pstate_v;}
		       }
		   }



	       for (p=0; p< A->len_aln; p++)T->seq_al[s][p]='.';
	       for (p=A->len_aln; p>0; p--)
	           {
		       if ( best_state_p==0)T->seq_al[s][p-1]=toupper(translate_dna_codon (A->seq_al[s]+(p-1), 'x'));
		       else if ( best_state_p>=SPL1  && best_state_p<=SPL4)T->seq_al[s][p-1]='-';
		       best_state_p=v_tab_p[p][best_state_p];
		   }
	       }



	   vfree (col);
	   return T;
       }

Alignment * mutate_cdna_aln ( Alignment *A)
{
    int a, b, c, n;
    int n1, n2, r1, r2;
    int **pos, ps;
    int neutral_substitution=50;
    int random_substitution=0;
    int random_deletion=0;
    int amino_acid_deletion=0;
    int amino_acid_substitution=0;
    char nuc_list[]="agct";
    char *new_codon;

    neutral_substitution=atoi(get_env_variable ("NEUTRAL_SUBSTITUTION",IS_FATAL));
    random_substitution =atoi(get_env_variable ("RANDOM_SUBSTITUTION", IS_FATAL));
    random_deletion     =atoi(get_env_variable ("RANDOM_DELETION", IS_FATAL));
    amino_acid_deletion =atoi(get_env_variable ("AMINO_ACID_DELETION", IS_FATAL));
    amino_acid_substitution =atoi(get_env_variable ("AMINO_ACID_SUBSTITUTION", IS_FATAL));


    if (A->S)free_sequence ( A->S, (A->S)->nseq);
    A->S=aln2seq(A);

    addrandinit(time (NULL));


    pos=aln2pos_simple ( A, A->nseq);

    /* 1 Apply neutral substitutions    */

    if ( neutral_substitution)
        {
	for (  c=0; c< neutral_substitution; c++)
	    {
	    for (  a=0; a< A->nseq; a++)
                {

		    for ( b=0; b< A->len_aln; b++)
		        {

			if (pos[a][b]<=0)continue;
			ps=MAX(0,pos[a][b]-(pos[a][b]-1)%3-1);


			n1=(A->S)->seq[a][pos[a][b]-1];
			r1=translate_dna_codon ( (A->S)->seq[a]+ps, 'o');

			n2=nuc_list[(int)addrand((unsigned long) 4)];
			(A->S)->seq[a][pos[a][b]-1]=n2;
			r2=translate_dna_codon ( (A->S)->seq[a]+ps, 'o');


			if ( r1==r2 && r1!='o')A->seq_al[a][b]=n2;

			else (A->S)->seq[a][pos[a][b]-1]=n1;
			}
		}
	    }
	}

    /* 2 Apply         substitutions    */
     if ( random_substitution)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue;
		    if (addrand ((unsigned long) 100)>random_substitution)continue;

		    n1=nuc_list[(int)addrand((unsigned long)4)];
		    (A->S)->seq[a][pos[a][b]-1]=n1;
		    A->seq_al[a][b]=n1;
		    }
	    }
	}

    /* 3 Apply amino acid substitutions */
      if ( amino_acid_substitution)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b+=3)
		    {
		    if (pos[a][b]<=0)continue;
		    if (addrand ((unsigned long) 100)>amino_acid_substitution)continue;
		    ps=MAX(0,pos[a][b]-(pos[a][b]-1)%3-1);

		    r1=translate_dna_codon ( (A->S)->seq[a]+ps, 'o');
		    new_codon=mutate_amino_acid(r1, "clustalw_col");

		    for ( c=ps; c<ps+3; c++)(A->S)->seq[a][c]=new_codon[c-ps];
		    }
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue;
		    else A->seq_al[a][b]=(A->S)->seq[a][pos[a][b]-1];
		    }
	    }
	}
    /* 3 Apply amino acid deletions     */
     if ( amino_acid_deletion)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b+=3)
		    {
		    if (pos[a][b]<=0)continue;
		    if (addrand ((unsigned long) 1000)>amino_acid_deletion)continue;
		    ps=MAX(0,pos[a][b]-(pos[a][b]-1)%3-1);
		    n=addrand ((unsigned long) 4)+1;

		    for ( c=ps; c<ps+(3*n) && c<A->len_aln; c++)(A->S)->seq[a][c]='-';
		    }
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue;
		    else A->seq_al[a][b]=(A->S)->seq[a][pos[a][b]-1];
		    }
	    }
	}
    /* 4 Apply amino acid insertions    */

/*FRAMESHIFT MUTATIONS*/
    /* 5 Apply nucleotide deletions*/
     if ( random_deletion)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue;
		    if (addrand ((unsigned long) 1000)>random_deletion)continue;

		    n1='-';
		    (A->S)->seq[a][pos[a][b]-1]=n1;
		    A->seq_al[a][b]=n1;
		    }
	    }
	}
    /* 6 Apply nucleotide deletions*/
     free_int (pos, -1);
   return A;

}

Alignment* clean_est  ( Alignment *A)
        {
	  /*Rules are as follow:
	    Internal Gap > 30% Requences ----> -
	    Best Residue < 50% Residues  ----> 'N'
	  */
	  int a, b,c;
	  int best;
	  int tot;

	  for ( a=0; a< A->len_aln; a++)
	    {

	      for (tot=0, b=0; b<4; b++)tot+=(A->P)->count[b][a];
	      best=best_int (5,1, &c, (A->P)->count[0][a],(A->P)->count[1][a],(A->P)->count[2][a],(A->P)->count[3][a],(A->P)->count[4][a]);

	      if ( tot==0)
		{
		  fprintf ( stderr, "\nWARNING: POSITION WITH NO INFORMATION [clean_est:%s]", PROGRAM);
		  A->seq_al[0][a]='-';
		}
	      else if (((A->P)->count[4][a]*100)/tot >30)A->seq_al[0][a]='-';
	      else if ( (best*100)/tot<50)A->seq_al[0][a]='n';

	    }
	return A;
	}



char **make_symbols ( char *name, int *n)
    {
    char **symbol;

    symbol=declare_char ( STRING, STRING);

    if ( strcmp (name, "3d_ali")==0)
        {
	sprintf ( symbol[0], "gih");
	sprintf ( symbol[1], "eb");
	sprintf ( symbol[2], "x");
	sprintf ( symbol[3], "#l");
	n[0]=4;
	}

    else if ( strcmp (name, "all")==0)
        {
	  int a, i;
	  for ( i=0,a=0; a<26; a++)
	    {
	      sprintf ( symbol[i++], "%c%c", 'a'+a, 'a'+a);
	      sprintf ( symbol[i++], "%c%c", 'A'+a, 'A'+a);
	    }
	  sprintf ( symbol[i++], "--");
	  n[0]=i;
	}

    else if ( strcmp (name, "set1")==0)
        {
	sprintf ( symbol[0], "ilvmfywhktcagH");
	sprintf ( symbol[1], "reqdnsP");
	sprintf ( symbol[2], "--");
	sprintf ( symbol[3], "#l");
	n[0]=4;
	}
    else if ( strcmp (name, "set2")==0)
        {
	n[0]=0;
	sprintf ( symbol[n[0]++], "gsacT");
	sprintf ( symbol[n[0]++], "ndtvpS");
	sprintf ( symbol[n[0]++], "ilkreqL");
	sprintf ( symbol[n[0]++], "--");
	sprintf ( symbol[n[0]++],"#l");
	}
    else if ( strcmp ( name, "any")==0)
        {
	sprintf ( symbol[0], "*x");
	n[0]=1;
       	}




    return symbol;
    }

char *testdna2gene (char *dna)
{
  int *w,a,l;
  l=strlen (dna);

  w=(int*)vcalloc(l+1, sizeof (int));
  for (a=0; a<l; a++)
    {
      w[a]=isupper (dna[a])?1:-1;
    }
  dna=dna2gene (dna,w);
  vfree (w);
  return dna;
}

Sequence *dnaseq2geneseq (Sequence *S, int **w)
{
  Sequence *PS;
  int a;
  char *p;

  PS=duplicate_sequence (S);
  for (a=0; a<S->nseq; a++)
    {
      p=dna2gene (S->seq[a], w[a]);
      if (strstr (p, "F"))
	{
	  HERE ("----FRAMESHIFT: %s", S->name[a]);
	}
      vfree (PS->seq[a]);
      PS->len[a]=strlen(p);
      PS->seq[a]=(char*)vcalloc (PS->len[a]+1, sizeof (char));
      sprintf ( PS->seq[a], "%s", p);
      vfree (p);
    }
  PS=reset_sequence_len (PS);
  return PS;
}

char *dna2gene (char *dna, int *w)
{
  int a, b, c, ns,l,od;
  int I1, I2, I3, START, NCE, NCS;
  int C1, S1_1, S2_1, S3_1, S4_1,NC1;
  int C2, S1_2, S2_2, S3_2, S4_2,NC2;
  int C3, S1_3, S2_3, S3_3, S4_3,NC3;
  int ST;
  int st;

  double p_C1, p_C2;
  double **C1_mat, **C2_mat;
  double *tb, **sc_mat, **tb_mat;
  double **em, **trans;
  double avg_w=0;

  char **string;
  double forbiden   =-100000;
  double frameshift1;
  double frameshift2;
  double exon_penalty;
  double exon_reward;
  double nostop_penalty;
  double shiftw;
  int frameshift_symbol='F';
  char *out_dna;
  int max=0;

  char *three_dna;

  three_dna=translate_dna_seq_on3frame (dna, 'x', NULL);
  lower_string(three_dna);
  l=strlen (dna);
  for (a=0; a<l; a++){max=MAX(max,w[a]);avg_w+=(double)w[a];}
  avg_w/=(double)l;
  shiftw=avg_w*-2;

  exon_penalty=-100*avg_w;
  exon_reward=avg_w;
  nostop_penalty=-100 *avg_w;
  frameshift1=forbiden;
  frameshift2=frameshift1;

  out_dna=(char*)vcalloc ( 2*strlen (dna)+1, sizeof (char));
  sprintf (out_dna, "%s", dna);
  ns=0;
  START=ns++;  I1=ns++;I2=ns++;I3=ns++;NCE=ns++;NCS=ns++;

  C1=ns++; S1_1=ns++;S2_1=ns++;NC1=ns++;S3_1=ns++;S4_1=ns++;
  C2=ns++; S1_2=ns++;S2_2=ns++;NC2=ns++;S3_2=ns++;S4_2=ns++;
  C3=ns++; S1_3=ns++;S2_3=ns++;NC3=ns++;S3_3=ns++;S4_3=ns++;
  ST=ns++;

  string=declare_char ( ns+1, 10);
  sprintf (string [S1_1], "S1_1");
  sprintf (string [S2_1], "S2_1");
  sprintf (string [S3_1], "S3_1");
  sprintf (string [S4_1], "S4_1");

  sprintf (string [S1_2], "S1_2");
  sprintf (string [S2_2], "S2_2");
  sprintf (string [S3_2], "S3_2");
  sprintf (string [S4_2], "S4_2");

  sprintf (string [S1_3], "S1_3");
  sprintf (string [S2_3], "S2_3");
  sprintf (string [S3_3], "S3_3");
  sprintf (string [S4_3], "S4_3");

  sprintf (string [START], "START");
  sprintf (string [NCE], "NCE");
  sprintf (string [NCS], "NCS");
  sprintf (string [NC1], "NC1");
  sprintf (string [NC2], "NC2");
  sprintf (string [NC3], "NC3");
  sprintf (string [I1], "I1");
  sprintf (string [I2], "I2");
  sprintf (string [I3], "I3");

  sprintf (string [C1], "C1");
  sprintf (string [C2], "C2");
  sprintf (string [C3], "C3");

  sprintf (string [ST], "STOP");

  trans=declare_double(ns,ns);
  em=declare_double   (ns,256);
  tb=(double*)vcalloc ( l+1, sizeof (double));
  sc_mat=declare_double (l+1, ns);
  tb_mat=declare_double (l+1, ns);
  C1_mat=declare_double (l+1, ns);
  C2_mat=declare_double (l+1, ns);

  for (a=0; a<ns; a++)
    {
      for (b=0; b<ns; b++)trans[a][b]=forbiden;
      for (b=0; b<256; b++)em[a][b]=forbiden;
    }

  trans[START][I1]=0;
  trans[START][NCS]=0;
  trans[NCS][NCS]=0;
  trans[NCS][NCE]=0;//allow sequence entirely non coding

  trans[NCS][I1]=0;
  trans[I1][I2]=0;
  trans[I2][I3]=0;
  trans[I3][C1]=0;

  trans[C1][C2]=0;
  trans[C1][C3]=frameshift1;
  trans[C1][C1]=frameshift2;

  trans[C2][C3]=0;
  trans[C2][C1]=frameshift1;
  trans[C2][C2]=frameshift2;


  trans[C3][C1]=exon_reward;
  trans[C3][C2]=frameshift1;
  trans[C3][C3]=frameshift2;

  trans[C3][NCE]=nostop_penalty;
  trans[C3][ST] =0; //normal termination
  trans[ST][NCE]=0;
  trans[NCE][NCE]=0;

  trans[C1][S1_1]=exon_penalty;
  trans[S1_1][S2_1]=0;
  trans[S2_1][NC1]=0;
  trans[NC1][NC1]=0;
  trans[NC1][S3_1]=0;
  trans[S3_1][S4_1]=0;
  trans[S4_1][C2]=0;
  trans[S4_1][C3]=frameshift1;
  trans[S4_1][C1]=frameshift2;


  trans[C2][S1_2]=exon_penalty;
  trans[S1_2][S2_2]=0;
  trans[S2_2][NC2]=0;
  trans[NC2][NC2]=0;
  trans[NC2][S3_2]=0;
  trans[S3_2][S4_2]=0;
  trans[S4_2][C3]=0;
  trans[S4_2][C1]=frameshift1;
  trans[S4_2][C2]=frameshift2;

  trans[C3][S1_3]=exon_penalty;
  trans[S1_3][S2_3]=0;
  trans[S2_3][NC3]=0;
  trans[NC3][NC3]=0;
  trans[NC3][S3_3]=0;
  trans[S3_3][S4_3]=0;
  trans[S4_3][C1]=0;
  trans[S4_3][C2]=frameshift1;
  trans[S4_3][C3]=frameshift2;

  em[I1]['A']=1;
  em[I2]['T']=1;
  em[I3]['G']=1;

  em[S1_1]['G']=1;
  em[S2_1]['T']=1;
  em[S3_1]['A']=1;
  em[S4_1]['G']=1;

  em[S1_2]['G']=1;
  em[S2_2]['T']=1;
  em[S3_2]['A']=1;
  em[S4_2]['G']=1;

  em[S1_3]['G']=1;
  em[S2_3]['T']=1;
  em[S3_3]['A']=1;
  em[S4_3]['G']=1;


  for (a=0; a<ns; a++)sc_mat[0][a]=tb_mat[0][a]=forbiden;
  sc_mat[0][START]=tb_mat[0][START]=0;

  for (a=1; a<=l ;a++)
    {
      int r;
      r=toupper (dna[a-1]);
      for (b=0; b<ns; b++)
	{
	  double best_sc,e,lw;
	  int best_t;

	  lw=(double)w[a-1]+shiftw;

	  if (b==ST && three_dna[a-1]=='x')e=0;
	  else if (b==C1 || b == C2 || b== C3)e=lw;
	  else if ( b==NC1|| b==NC2 || b==NC3 || b==NCE || b==NCS)e=-lw;
	  else e=em[b][r];

	  if (e==forbiden)
	  {
		sc_mat[a][b]=forbiden;
		tb_mat[a][b]=0;
	  }
	  else
	     {
	       for (best_sc=forbiden,best_t=0,c=0; c<ns; c++)
		 {
		   double tr, sc, p_sc;
		   tr  = trans[c][b];
		   p_sc=sc_mat[a-1][c];

		   //Frameshift handling

		   if ( tr== forbiden || p_sc==forbiden);
		   else if (tr!=forbiden)
		     {
		       if     (b==C2 && c!=C1 && c!=S4_1){p_C1='N'; p_C2=r;  }
		       else if(b==C3 && c!=C2 && c!=S4_2){p_C1='N'; p_C2='N';}
		       else
			 {
			   p_C1=C1_mat[a-1][c];
			   p_C2=C2_mat[a-1][c];
			 }

		       if (b==C3 && is_stop (p_C1, p_C2,r)){tr=forbiden;}
		       else
			 {
			   sc=tr+e+p_sc;
			   if (c==0 || sc>best_sc)
			     {
			       best_sc=sc;
			       best_t =c;
			     }
			 }
		     }
		 }

	       C1_mat[a][b]=(b==C1)?r:C1_mat[a-1][best_t];
	       C2_mat[a][b]=(b==C2)?r:C2_mat[a-1][best_t];
	       sc_mat[a][b]=best_sc;
	       tb_mat [a][b]=best_t;
	     }

	}
    }

  a=l;
  b=NCE;
  c=sc_mat[a][NCE];

  while (a>0)
    {
//      HERE ("**%d [%s] %c in %d", b,string[b], dna[a-1], a);
      tb[a]=b;
      b=tb_mat[a][b];
      a--;
    }

  od=0;
  st=0;
  for (a=0;a<l;a++)
    {
      int r,t,pt, coding;
      t =tb[a+1];
      pt=tb[a ];
      r=dna[a];
      if ( st || t==ST)st++;
      coding=(t==C1 || t==C2 || t==C3 || t==I1 || t==I2 || t==I3 ||(st && st<=3))?1:0;



      if      (t==C1 && (pt==C2 || pt==S4_2)){out_dna[od++]=frameshift_symbol;}
      else if (t==C1 && (pt==C1 || pt==S4_1)){out_dna[od++]=frameshift_symbol;out_dna[od++]=frameshift_symbol;}

      else if (t==C2 && (pt==C3 || pt==S4_3)){out_dna[od++]=frameshift_symbol;}
      else if (t==C2 && (pt==C2 || pt==S4_2)){out_dna[od++]=frameshift_symbol;out_dna[od++]=frameshift_symbol;}

      else if (t==C3 && (pt==C1 || pt==S4_1)){out_dna[od++]=frameshift_symbol;}
      else if (t==C3 && (pt==C3 || pt==S4_3)){out_dna[od++]=frameshift_symbol;out_dna[od++]=frameshift_symbol;}

      if (coding)out_dna[od++]=toupper (r);
      else out_dna[od++]=tolower(r);

    }

  free_double (tb_mat, -1);
  free_double (sc_mat, -1);
  free_double (trans, -1);
  free_double (em, -1);
  free_double (C1_mat, -1);
  free_double (C2_mat, -1);
  vfree (tb);
  vfree (three_dna);
  return out_dna;
}


int res_weights2avg(Sequence *R, int **w)
{
  int a, b;
  double avg=0;
  int n=0;

  for (a=0; a<R->nseq; a++)
    for (b=0; b<R->len[a]; b++){avg+=w[a][b];n++;}
  return avg/n;
}
int res_weights2min(Sequence *R, int **w)
{
  int a,b;
  int min=w[0][0];

  for (a=0; a<R->nseq; a++)
    for (b=0; b<R->len[a]; b++)min=MIN(min,(w[a][b]));
  return min;
}
int res_weights2max(Sequence *R, int **w)
{
  int a, b;
  int max=w[0][0];

  for (a=0; a<R->nseq; a++)
    for (b=0; b<R->len[a]; b++)max=MAX(max,(w[a][b]));
  return max;
}
int scan_res_weights4ac (Sequence *R, int **w, int start, int end, int step)
{
  int best_t, a;
  float best_ac=0;
  float *count;
  float *acc;
  int avg;

  avg=res_weights2avg(R,w);
  best_t=start;
  for (a=start; a<=end; a+=step)
    {

      count=res_weights2accuracy_counts (R,w,a,NULL);
      acc=counts2accuracy (count);

      if (acc[3]>best_ac)
	{
	  best_ac=acc[3];
	  best_t=a;
	}
      vfree (count);vfree (acc);
    }
  count=res_weights2accuracy_counts (R,w,best_t,NULL);
  acc=counts2accuracy (count);
  fprintf (stderr, "\nBest_T: %d ", best_t);
  display_accuracy (count,stderr);

  count=res_weights2accuracy_counts (R,w,2*avg,NULL);
  acc=counts2accuracy (count);
  fprintf (stderr, " Avg_T: %d ", 2*avg);
  display_accuracy (count,stderr);

  vfree (acc); vfree (count);
  return best_t;
}
int ** shift_res_weights ( Sequence *R, int **w, int shift)
{
  int a, b;
  for (a=0; a<R->nseq; a++)
    for (b=0; b<R->len[a]; b++)
      w[a][b]+=shift;
  return w;
}
float *res_weights2accuracy_counts ( Sequence *R, int **w,int T, float *result)
{
  int a, b, coding,pcoding;

  if (!result)result=(float*)vcalloc (4, sizeof (float));

  for (a=0; a<R->nseq; a++)
    {
      for (b=0; b<R->len[a]; b++)
	{
	  coding=(isupper(R->seq[a][b]))?1:0;
	  pcoding=(w[a][b]>T)?1:0;
	  if      (  coding &&  pcoding)result[0]++;//TP
	  else if ( !coding && !pcoding)result[1]++;//TN
	  else if ( !coding &&  pcoding)result[2]++;//FP
	  else if (  coding && !pcoding)result[3]++;//FN
	}
    }
  return result;
}

void genepred_seq2accuracy_counts4all (Sequence *R, Sequence *T)
{
  int a,b;
  float *result =(float*)vcalloc (4, sizeof (float));

  fprintf ( stderr, "\n");

  for (a=0; a<R->nseq; a++)
  {
    fprintf ( stderr, "gene: %s ", R->name[a]);
    for (b=0; b<T->nseq; b++)
    {
      if ( strm (R->name[a], T->name[b]) && hasupper(R->seq[a]))
      {
	vfree (display_accuracy (genepred2accuracy_counts (R->seq[a], T->seq[b], NULL),stderr));
	break;
      }
    }
  }
  vfree(result);
}

float* genepred_seq2accuracy_counts (Sequence *R, Sequence *T,float *result)
{
  int a,b;

  if (!result)result=(float*)vcalloc (4, sizeof (float));

  for (a=0; a<R->nseq; a++)
    for (b=0; b<T->nseq; b++)
      if ( strm (R->name[a], T->name[b]) && hasupper(R->seq[a]))
	genepred2accuracy_counts (R->seq[a], T->seq[b], result);
  return result;
}

float* genepred2accuracy_counts      (char *ref,  char *target , float *result)
{
  char *ref2, *target2;
  int l,a;
  if ( !result) result=(float*)vcalloc (4, sizeof (float));
  ref2=(char*)vcalloc ( strlen (ref)+1, sizeof (char));
  sprintf ( ref2, "%s", ref);

  target2=(char*)vcalloc ( strlen (target)+1, sizeof (char));
  sprintf ( target2, "%s", target);

  remove_charset (ref2, "Ff");
  remove_charset (target2, "Ff");

  if ( strlen (target2) != strlen (ref2))
    {fprintf (stderr, "ERROR: Gene and target have different length [FATAL]\n"); myexit (EXIT_FAILURE);}

  l=strlen (ref2);
  for (a=0; a<l; a++)
    {
      int coding, pcoding;
      coding =isupper (ref2[a]);
      pcoding=isupper (target2[a]);
      if      ( coding  &&  pcoding)result[0]++;//TP
      else if ( !coding && !pcoding)result[1]++;//TN
      else if ( !coding &&  pcoding)result[2]++;//FP
      else if (  coding && !pcoding)result[3]++;//FN
    }

  vfree (ref2);
  vfree (target2);
  return result;
 }

int is_stop( char r1, char r2, char r3)
{
  char codon[4];

  if (!r2 || !r3) return 0;
  else if (tolower (r1)=='n' || tolower(r2)=='n' || tolower(r3)=='n') return 0;
  else
    {
      sprintf (codon, "%c%c%c", tolower(r1), tolower(r2), tolower(r3));
      if (translate_dna_codon (codon, 'x')=='x')return 1;
      else return 0;
    }
}


char * translate_dna_seq_on3frame (  char *dna_seq, char stop, char *prot)
       {
	  int a, l;
	  char *buf;

	  l=strlen (dna_seq);
	  if ( prot==NULL)prot=(char*)vcalloc ( l+2, sizeof (char));

	   buf=(char*)vcalloc (l+4, sizeof (char));
	   sprintf (buf, "%s", dna_seq);
	   lower_string ( buf);
	   for ( a=0; a< l; a++)buf[a]=(buf[a]=='t')?'u':buf[a];

	   for (a=0; a< l; a++)
	       prot[a]=translate_dna_codon (buf+a, stop);
	   vfree (buf);
	   prot[a]='\0';

	   return prot;
       }
char * translate_dna_seq ( char *dna_seq, int frame, char stop, char *prot)
       {
	 //frame: 1->3
	   int a, b, l;
	   char *buf;

           l=strlen (dna_seq);
	   if ( prot==NULL)prot=(char*)vcalloc ( l, sizeof (char));
	   frame--;
	   buf=(char*)vcalloc (l+4, sizeof (char));
	   sprintf (buf, "%s", dna_seq);
	   lower_string ( buf);
	   for ( a=0; a< l; a++)buf[a]=(buf[a]=='t')?'u':buf[a];


	   for ( b=0,a=0+frame; a< l; a+=3,b++)
	     {

	       prot[b]=translate_dna_codon (buf+a, stop);
	     }
	   vfree (buf);
	   prot[b]='\0';
	   upper_string (buf);
	   return prot;
       }
char * back_translate_dna_codon ( char aa, int deterministic)
        {
	static char *r;
	int choice;

	vsrand(0);
	if ( r==NULL)r=(char*)vcalloc (4, sizeof (char));

	if (!is_gap(aa))aa=tolower(aa);

	if (is_gap(aa))sprintf (r, "---");
	else if ( aa>=0 && aa<=9)
	  {
	    sprintf (r, "%d%d%d", aa, aa,aa);
	  }
	else if ( aa>='0' && aa<='9')
	  {
	     sprintf (r, "%c%c%c", aa, aa,aa);
	  }
	else if ( aa=='a')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if      ( choice==0)sprintf (r, "gca");
	    else if ( choice==1)sprintf (r, "gcg");
	    else if ( choice==2)sprintf (r, "gcc");
	    else if ( choice==3)sprintf (r, "gct");
	  }
	else if ( aa=='c')
	  {
	   choice=(deterministic)?0:rand()%2;
	    if      ( choice==0)sprintf (r, "tgc");
	    else if ( choice==1)sprintf (r, "tgt");
	  }
	else if ( aa=='d')
	  {
	  choice=(deterministic)?0:rand()%2;
	  if ( choice==0)sprintf (r, "gac");
	  else if ( choice==1)sprintf (r, "gat");
	  }

	else if ( aa=='e')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if ( choice==0)sprintf (r, "gaa");
	    else sprintf (r, "gag");
	  }
	else if ( aa=='f')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if ( choice==0)sprintf (r, "ttc");
	    else sprintf (r, "ttt");
	  }
	else if ( aa=='g')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "gga");
	    else if ( choice==1) sprintf (r, "ggg");
	    else if ( choice==2) sprintf (r, "ggc");
	    else if ( choice==3) sprintf (r, "ggt");
	  }
	else if ( aa=='h')
	  {
	    choice =rand()%2;
	    if ( choice==0)sprintf (r, "cac");
	    else sprintf (r, "cat");
	  }
	else if ( aa=='i')
	  {
	    choice=(deterministic)?0:rand()%3;
	    if  ( choice==0)     sprintf (r, "ata");
	    else if ( choice==1) sprintf (r, "atc");
	    else if ( choice==2) sprintf (r, "att");
	  }
	else if ( aa=='k')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "aaa");
	    else if ( choice==1) sprintf (r, "aag");

	  }
	else if ( aa=='l')
	  {
	    choice=(deterministic)?0:rand()%6;
	    if  ( choice==0)     sprintf (r, "cta");
	    else if ( choice==1) sprintf (r, "ctg");
	    else if ( choice==2) sprintf (r, "ctc");
	    else if ( choice==3) sprintf (r, "ctt");
	    else if ( choice==4) sprintf (r, "tta");
	    else if ( choice==5) sprintf (r, "ttg");
	  }
	else if ( aa=='m')sprintf ( r, "atg");
	else if ( aa=='n')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "aac");
	    else if ( choice==1) sprintf (r, "aat");
	  }
	else if ( aa=='p')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "cca");
	    else if ( choice==1) sprintf (r, "ccg");
	    else if ( choice==2) sprintf (r, "ccc");
	    else if ( choice==3) sprintf (r, "cct");
	  }
	else if ( aa=='q')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "caa");
	    else if ( choice==1) sprintf (r, "cag");
	  }
        else if ( aa=='r')
	  {
	    choice=(deterministic)?0:rand()%6;
	    if  ( choice==0)     sprintf (r, "cga");
	    else if ( choice==1) sprintf (r, "cgg");
	    else if ( choice==2) sprintf (r, "cgc");
	    else if ( choice==3) sprintf (r, "cgt");
	    else if ( choice==4) sprintf (r, "aga");
	    else if ( choice==5) sprintf (r, "agg");

	  }
	else if ( aa=='s')
	  {
	    choice=(deterministic)?0:rand()%6;
	    if  ( choice==0)     sprintf (r, "tca");
	    else if ( choice==1) sprintf (r, "tcg");
	    else if ( choice==2) sprintf (r, "tcc");
	    else if ( choice==3) sprintf (r, "tct");
	    else if ( choice==4) sprintf (r, "agt");
	    else if ( choice==5) sprintf (r, "agc");

	  }
	else if ( aa=='t')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "aca");
	    else if ( choice==1) sprintf (r, "acg");
	    else if ( choice==2) sprintf (r, "acc");
	    else if ( choice==3) sprintf (r, "act");
	  }
	else if ( aa=='v')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "gta");
	    else if ( choice==1) sprintf (r, "gtg");
	    else if ( choice==2) sprintf (r, "gtc");
	    else if ( choice==3) sprintf (r, "gtt");
	  }
	else if ( aa=='w')
	  {
	    sprintf (r, "tgg");
	  }
	else if ( aa=='y')
	  {
	     choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "tac");
	    else if ( choice==1) sprintf (r, "tat");
	  }
	else
	  {
	    sprintf (r, "nnn");
	  }
	return r;

	}
int translate_dna_codon ( char *sequence, char stop)
        {
	char seq[4];
	int b;
	int upper;
	int ret;

	if ( strlen (sequence)<1) return 'n';
	upper=isupper (sequence[0])?1:0;

	if ( (b=strlen (sequence))<3)
	  ret='x';
	else
	  {
	    seq[0]=tolower(sequence[0]);
	    seq[1]=tolower(sequence[1]);
	    seq[2]=tolower(sequence[2]);
	    seq[3]='\0';

	    seq[0]=(seq[0]=='u')?'t':seq[0];
	    seq[1]=(seq[1]=='u')?'t':seq[1];
	    seq[2]=(seq[2]=='u')?'t':seq[2];
	    if (strm (seq, "---"))return '-';
	    else if ( strm5(seq, "gca", "gcg", "gcc", "gct","gcn"))ret='a';
	    else if ( strm2(seq, "tgc","tgt"))ret='c';
	    else if ( strm2(seq, "gac","gat"))ret='d';
	    else if ( strm2(seq, "gaa","gag"))ret='e';
	    else if ( strm2(seq, "ttc","ttt"))ret='f';
	    else if ( strm5(seq, "gga","ggg","ggc", "ggt", "ggn"))ret='g';
	    else if ( strm2(seq, "cac","cat"))ret='h';
	    else if ( strm3(seq, "ata","atc","att"))ret='i';
	    else if ( strm2(seq, "aaa","aag"))ret= 'k';
	    else if ( strm6(seq, "cta","ctg","ctc", "ctt", "tta", "ttg"))ret='l';
	    else if ( strm (seq, "ctn"))ret='l';
	    else if ( strm (seq, "atg"))ret='m';
	    else if ( strm2(seq, "aac","aat"))ret= 'n';
	    else if ( strm5(seq, "cca","ccg","ccc", "cct","ccn"))ret='p';
	    else if ( strm2(seq, "cag","caa"))ret='q';
	    else if ( strm6(seq, "cga","cgg","cgc", "cgt","aga","agg"))ret='r';
	    else if ( strm (seq, "cgn"))ret= 'r';
	    else if ( strm6(seq, "tca","tcg","tcc", "tct","agc","agt"))ret='s';
	    else if ( strm (seq, "ccn"))ret='s';
	    else if ( strm5(seq, "aca","acg","acc", "act", "acn"))ret='t';
	    else if ( strm5(seq, "gta","gtg","gtc", "gtt", "gtn"))ret='v';
	    else if ( strm (seq, "tgg"))ret='w';
	    else if ( strm2(seq, "tac","tat"))ret='y';
	    else if ( strm3(seq, "tag","taa","tga"))ret=stop;
	    else 
	      {
		ret='x';
	      }

	    ret= (upper)?toupper(ret):ret;
	  }
	return ret;
	}

int extend_seqaln (Sequence *S, Alignment *A)
{
  char **s;
  int n,a;
  if (S){s=S->seq;n=S->nseq;}
  else if (A){s=A->seq_al;n=A->nseq;}
  else return 0;

  for (a=0; a<n;a++){extend_seq(s[a]);}
  return 1;
}
int unextend_seqaln (Sequence *S, Alignment *A)
{
  char **s;
  int n, a;
  if (S){s=S->seq;n=S->nseq;}
  else if (A){s=A->seq_al;n=A->nseq;}
  else return 0;

  for (a=0; a<n;a++){unextend_seq(s[a]);}
  return 1;
}


char *extend_seq (char *seq)
{
  char *buf, *ebuf;
  int l, lb, a, b, upper,v;
  char r1, r2;

  l=strlen (seq);
  buf =(char*)vcalloc ( l+1, sizeof (char));
  ebuf=(char*)vcalloc ( l+1, sizeof (char));
  sprintf (  buf, "%s", seq);
  sprintf ( ebuf, "%s", seq);

  ungap ( buf);
  ungap (ebuf);
  lb=strlen (buf);

  for (a=0; a<lb-1; a++)
    {
      r1=buf[a];
      r2=buf[a+1];

      upper=(isupper(r1))?1:0;
      r1=tolower(r1);
      r2=tolower(r2);

      r1=(r1=='u')?'t':r1;
      r2=(r2=='u')?'t':r2;
      if (r1=='x' || r1=='n')v='x';
      else if (r2=='n' || r2=='x')v=r1;

      else if (r1=='a' && r2=='a')v='d';
      else if (r1=='a' && r2=='c')v='e';
      else if (r1=='a' && r2=='g')v='f';
      else if (r1=='a' && r2=='t')v='h';

      else if (r1=='c' && r2=='a')v='i';
      else if (r1=='c' && r2=='c')v='k';
      else if (r1=='c' && r2=='g')v='l';
      else if (r1=='c' && r2=='t')v='m';

      else if (r1=='g' && r2=='a')v='n';
      else if (r1=='g' && r2=='c')v='p';
      else if (r1=='g' && r2=='g')v='q';
      else if (r1=='g' && r2=='t')v='r';

      else if (r1=='t' && r2=='a')v='s';
      else if (r1=='t' && r2=='c')v='v';
      else if (r1=='t' && r2=='g')v='w';
      else if (r1=='t' && r2=='t')v='y';
      else
	{
	  v='j';
	}
      ebuf[a]=(upper)?toupper(v):v;
    }

  for (b=0,a=0; a<l; a++)
    {
      if ( !is_gap(seq[a]))seq[a]=ebuf[b++];
    }
  vfree (ebuf);
  vfree (buf);
  return seq;
}
char *unextend_seq (char *seq)
{
  char *buf, *ebuf;
  int l, lb, a, b, upper,v;
  char r1, r2;

  l=strlen (seq);
  buf =(char*)vcalloc ( l+1, sizeof (char));
  ebuf=(char*)vcalloc ( l+1, sizeof (char));
  sprintf (  buf, "%s", seq);
  sprintf ( ebuf, "%s", seq);

  ungap ( buf);
  ungap (ebuf);
  lb=strlen (buf);

  for (a=0; a<lb-1; a++)
    {
      r1=buf[a];
      upper=(isupper(r1))?1:0;
      r1=tolower(r1);
      r1=(r1=='u')?'t':r1;

      if (r1=='x')v='n';
      else if (r1=='a' || r1=='d' || r1=='e' || r1 == 'f' || r1 == 'h')v='a';
	  else if (r1=='c' || r1=='i' || r1=='k' || r1 == 'l' || r1 == 'm')v='c';
	  else if (r1=='g' || r1=='n' || r1=='p' || r1 == 'q' || r1 == 'r')v='g';
	  else if (r1=='t' || r1=='s' || r1=='v' || r1 == 'w' || r1 == 'y')v='t';
      else
	  {
		  v='j';
	  }

      ebuf[a]=(upper)?toupper(v):v;
    }

  for (b=0,a=0; a<l; a++)
    {
      if ( !is_gap(seq[a]))seq[a]=ebuf[b++];
    }
  vfree (ebuf);
  vfree (buf);
  return seq;
}



Alignment * mutate_aln ( Alignment *A, char *r)
{
  int a, b, c, mut,type, ratio;
  char alp[30];
  int alp_size;
  Sequence *S;
  Alignment*B;
  int n_mut, tot;

  vsrand(0);
  if ( r[0]=='\0')ratio=0.01*RAND_MAX;
  else ratio=atof(r)*RAND_MAX;

  S=aln2seq(A);
  S=get_sequence_type(S);



  if ( strm(S->type, "DNA") ||  strm(S->type, "RNA"))sprintf (alp, "AGCT");
  else if (  strm(S->type, "PROTEIN"))sprintf (alp, "ACDEFGHIKLMNPQRSTVWY");

  alp_size=strlen(alp);

  B=copy_aln (A,NULL);
  B=realloc_aln(B, B->len_aln*2+1);

  for ( a=0, b=0; a< A->len_aln; a++, b+=2)
    {
      for ( c=0; c< A->nseq; c++)
	{
	  B->seq_al[c][b]=tolower(A->seq_al[c][a]);
	  B->seq_al[c][b+1]='~';
	}
    }

  for ( c=0; c< A->nseq; c++)B->seq_al[c][b]='\0';
  B->len_aln=A->len_aln*2;



  tot=n_mut=0;
  for (a=0; a< B->len_aln; a+=2)
    for ( b=0; b<B->nseq; b++)
      {
	if ( is_gap(B->seq_al[b][a]))continue;
	mut=((rand()%RAND_MAX)>ratio)?0:1;
	tot++;
	n_mut+=mut;

	if (mut)
	  {
	    type=rand()%2;
	    if (type==0)/*deletion*/
	      {
		B->seq_al[b][a]='.';
	      }
	    else if ( type==1)
	      {
		B->seq_al[b][a+1]=alp[rand()%alp_size];
	      }
	    else if (type==2)
	      {
		B->seq_al[b][a]=alp[rand()%alp_size];
	      }

	  }
      }
  ungap_aln (B);


  free_sequence (S, S->nseq);
  free_aln (A);
  return B;

}

char* mutate_amino_acid ( char aa, char *mode)

     {
	 int a, b, c, d;
	 char nucleotide[]="agct";
	 char amino_acid[]="acdefghiklmnpqrstvwy";
	 static char **triplet;
	 static char **cw_col;
	 int ng_cw_col;
	 static int **amino_acid_list;
	 static int *lu;
	 char a1, a2;
	 char *mat;

	 aa=tolower(aa);
	 declare_name(mat);
	 if ( !mode)sprintf (mat, "clustalw_col");
	 else sprintf (mat, "%s", mode);
	 if (!triplet)
	    {
		triplet=declare_char ( 64, 4);
		for (d=0, a=0; a< 4;a++)
		    for ( b=0; b< 4; b++)
			for ( c=0; c< 4; c++, d++)
			    {
				triplet[d][0]=nucleotide[a];
				triplet[d][1]=nucleotide[b];
				triplet[d][2]=nucleotide[c];
			    }
	    }
	 if ( !cw_col)cw_col=make_group_aa ( &ng_cw_col,mat);
	 if ( !amino_acid_list)
	    {
		amino_acid_list=declare_int ( 20, 65);
		for ( a=0; a< 20; a++)
		    for ( b=0; b< 64; b++)
		        {
			    a1=translate_dna_codon ( triplet[b], 'x');
			    a2=amino_acid[a];
			    for ( d=0; d< ng_cw_col; d++)
				if ( is_in_set ( a1, cw_col[d]) && is_in_set ( a2, cw_col[d]))
				   {
				       amino_acid_list[a][++amino_acid_list[a][0]]=b;
				   }
			}
		lu=(int*)vcalloc ( 26, sizeof (int));
		for ( a=0; a<20; a++)
		    {
			lu[amino_acid[a]-'a']=a;
		    }
		/*
		for ( a=0; a< 20; a++)
		    {
			fprintf ( stderr, "\n%c", amino_acid[a]);
			for ( b=1; b<=amino_acid_list[a][0]; b++)
			    fprintf ( stderr, "\n\t%s %c", triplet[amino_acid_list[a][b]], translate_dna_codon (triplet[amino_acid_list[a][b]], 'x'));
		    }
		*/
	    }

	 return triplet [addrand((unsigned long)amino_acid_list[lu[aa-'a']][0])+1];
     }

/**************************************************************************************************/
/********************************                      ********************************************/
/********************************    PROCESSING        ********************************************/
/********************************                      ********************************************/

int ls_compare (const void * a, const void * b)
{
	return strcmp(*(char**)a,*(char**)b);
}

void modify_data  (Sequence_data_struc *D1in, Sequence_data_struc *D2in, Sequence_data_struc *DSTin, char **action_list,int n_actions, Action_data_struc *RAD)
     {
       Sequence  *COOR=NULL, *NS=NULL,*BUFS=NULL, *OUT_S=NULL;
       Constraint_list *CL;
       char *s;
       int value,upper_value, lower_value, start, end, a, b,c;
       int *count_table=NULL;
       char *action;
       Sequence_data_struc *D1;
       Sequence_data_struc *D2;
       Sequence_data_struc *DST;
       int s1, s2, r1, r2;
       
       Alignment *BUF;
       
       //Switches
       static int evaluate2tree;
       static int clean_flag;
       static int gtree=0;
       /*Switches*/
       
       action=action_list[0];

       if (action[0]=='2')
	 {

	   D1=D2in;
	   D2=D1in;
	   DST=DSTin;
	   action++;
	 }
       else if ( action[0]=='1')
	 {
	   D1=D1in;
	   D2=D2in;
	   DST=DSTin;
	   action++;
	 }
       else if ( action[0]=='3')
	 {
	   D1=DSTin;
	   D2=D1in;
	   DST=DSTin;
	   action++;
	 }
       else
	 {
	   D1=D1in;
	   D2=D2in;
	   DST=DSTin;
	 }
       if (!D1->A)D1->A=copy_aln (D1in->A, NULL);

       if (  strm(action, "seqnos"))
	 {
	  (D1->A)->output_res_num=1;
	 }
       else if      ( strm (action, "pdb2color"))
	  {
	    int nn=1;
	    int sl=strlen ((D1->S)->seq[0]);
	    char  *file=(D1->S)->file[0];
	    float *v=(float*)vcalloc (count_n_line_in_file(file), sizeof (float));
	    
	    nn=1;
	    while (ACTION(nn))
	      {
		int coor=atoi (action_list[nn++]);
		v[coor]=99.99;
			       
	      }
	    
	    bfactor2x_in_pdb (file, "stdout", v);
	    exit (0);
	  }
       else if ( strm (action,"aln2bootstrap"))
	 {
	   (D1->A)=aln2bootstrap (D1->A, ATOI_ACTION (1));
	   D1->S=aln2seq (D1->A);
	 }
      
       else if ( strm (action,"aln2sample"))
	 {
	   (D1->A)=aln2sample (D1->A, ATOI_ACTION (1));
	   D1->S=aln2seq (D1->A);
	 }
       else if ( strm (action,"aln2random_aln"))
	 {
	   (D1->A)=aln2random_aln (D1->A, ACTION (1));
	   D1->S=aln2seq (D1->A);
	 }
       else if ( strm (action, "or_scan"))
	 {
	   HERE ("OR SCAN");
	   D1->A=or_scan(D1->A, D2->A, ACTION(1));
	   D1->S=aln2seq (D1->A);
	 }
       else if ( strm (action, "or_sar"))
	 {
	   D1->A=or_sar(D1->A, D2->A, ACTION(1), PRINT);
	   D1->S=aln2seq (D1->A);
	 }
       else if ( strm ( action, "sar2subsar"))
	 {
	   /*in->sequences
	     in2->sar data
	   */
	   Alignment *subA, *subS;

	   if ( n_actions==1)
	     {
	       fprintf ( stderr, "\nin=aln, in2=sar sar2subsar [filter value compound1 compound2...] | [jack1] | [file]\n");
	       myexit (EXIT_FAILURE);
	     }

	   sarset2subsarset ( D1->A, D2->A, &subA, &subS, main_read_aln (action_list[2], NULL));
	   D1->A=subA;D2->A=subS;
	 }
       else if ( strm (action, "display_sar"))
	 {
	   D1->A=display_sar (D1->A, D2->A, action_list[1]);
	 }
       else if ( strm ( action, "sar2simpred"))
	 {
	   /*in->sequences
	     in2->sar data
	   */
	   sar2simpred ( D1->A, D2->A, action_list[1], action_list[2], atoi(action_list[3]), atoi (action_list[4]));
	 }
       else if ( strm ( action, "sar2simpred2"))
	 {
	   /*in->sequences
	     in2->sar data
	   */
	   if ( n_actions!=5)
	     {
	       fprintf ( stderr, "\nERROR: +sar2simpred2 seqnamesfile posfile compound limit");
	       myexit (EXIT_FAILURE);
	     }
	   sar2simpred2 ( D1->A, D2->A, action_list[1], action_list[2], action_list[3], atoi (action_list[4]));
	 }
        else if ( strm ( action, "sar_analyze"))
	 {
	   /*in->sequences
	     in2->sar data
	   */
	   sar_analyze ( D1->A, D2->A,action_list[1]);
	 }
       	else if ( strm ( action, "simple_sar_predict"))
	  {
	    //displays each column with ist score;
	    simple_sar_predict (D1->A, D2->A,ACTION(1));
	    myexit (EXIT_SUCCESS);
	  }
	else if ( strm ( action, "display_sar_analyze"))
	  {
	    //displays each column with ist score;
	    display_simple_sar_analyze_col (D1->A, D2->A,ACTION(1));
	    myexit (EXIT_SUCCESS);
	  }
       else if ( strm ( action, "display_sar_analyze_pc"))
	  {
	    //displays each column with ist score;
	    display_simple_sar_analyze_pair_col (D1->A, D2->A,ACTION(1));
	    myexit (EXIT_SUCCESS);
	  }
       else if ( strm ( action, "weight2sar"))
	 {
	   /*in->sequences
	     in2->sar data
	   */
	   if ( n_actions!=3)
	     {
	       fprintf ( stderr, "\nERROR: +weight2sar <weight_file> <limit>");
	       myexit (EXIT_FAILURE);
	     }
	   D1->A=weight2sar ( D1->A,D2->A, action_list[1], atoi(action_list[2]));

	 }
	else if ( strm ( action, "sar_weight"))
	 {
	   /*in->sequences
	     in2->sar data
	   */
	   if ( n_actions!=3)
	     {
	       fprintf ( stderr, "\nERROR: +sar_weight <sar_analyze> <compound>");
	       myexit (EXIT_FAILURE);
	     }
	   D1->A=aln2weighted_sar_score ( D1->A,D2->A, action_list[1], action_list[2]);
	   D1->S=aln2seq ( D1->A);
	 }

       else if ( strm (action, "name2unique_name"))
	 {
	   char *tmp1, *tmp2;
	   char command[1000];
	   tmp1=vtmpnam (NULL); tmp2=vtmpnam (NULL);

	   output_fasta_aln (tmp1,D1->A);
	   free_aln (D1->A);free_sequence (D1->S, -1);
	   sprintf ( command, "fasta_aln2fasta_aln_unique_name.pl %s >%s", tmp1, tmp2);
	   my_system ( command);
	   D1->S=get_fasta_sequence ( tmp2, NULL);
	   D1->A=seq2aln (D1->S,NULL, 1);
	 }
       else if ( strm (action, "rm_tag") || strm (action, "rm_template"))
	 {

	   char **temp_name=NULL,**temp_list=NULL, temp_nseq=0;
	   int z;

	   if ( D1 && D1->A){temp_name=(D1->A)->name;temp_nseq=(D1->A)->nseq;}
	   else if ( D1 && D1->S){temp_name=(D1->S)->name;temp_nseq=(D1->S)->nseq;}
           temp_list=rm_name_tag (temp_name,temp_nseq, NULL);
	   if ( n_actions>1 && strm (action_list[1], "template"))
	      {

               for ( z=0; z<temp_nseq; z++)
		{
		if (temp_list[z][0])
			{fprintf (stdout, "%s\n", temp_list[z]);}
	       	}
	      	myexit (EXIT_SUCCESS);
	      }
	 }
       else if (strm (action, "add_template") || strm (action, "swap_header") || strm (action, "seq2template"))
	 {
	   int n=1;
	   while (ACTION(n))
	     {
	       D1->S=seq2template_seq (D1->S, action_list[n++], NULL);
	     }
	   D1->A=seq2aln(D1->S, NULL, 1);
	 }
       else if ( strm ( action, "seq2year"))
	 {
	   D1->S=seq2year (D1->S, (n_actions>1)?atoi(action_list[1]):1);
	   D1->A=seq2aln(D1->S, NULL, 1);
	 }
       else if ( strm (action, "swap_lib_header"))
	 {
	   Sequence *S;
	   S=main_read_seq (action_list[1]);
	   (D1->CL)->S=S;

	 }
       else if ( strm (action, "weight_lib"))
	 {
	   int l;
	   int w;
	   w=atoi (action_list[1]);
	   if ( D1->CL)
	     {
	       int s1, s2,r1,r2;
	       Sequence *S=(D1->CL)->S;
	       int ***r=(D1->CL)->residue_index;

	       for (s1=0; s1<S->nseq; s1++)
		 for (r1=1; r1<=S->len[s1]; r1++)
		   for (b=1; b<r[s1][r1][0]; b+=3)
		     {
		       r[s1][r1][b+2]=w;
		     }
	     }
	 }
       else if ( strm (action, "struc2nb"))
	 {
	   int c;
	   for ( c=0; c< (D1->S)->nseq; c++)
	     {
	       struclist2nb ((D1->S)->name[c],(D1->S)->seq[c], (D1->S)->seq_comment[c], atof(action_list[1]),ACTION(2),ACTION(3) );
	     }
	   myexit (EXIT_SUCCESS);
	 }
       
       else if ( strm (action, "pdb2contacts"))
	 {
	   float D=0;
	   char mode[100];
	   char type[100];
	   Sequence *T=D2?(D2->S):NULL;
	   
	   
	   
	   if (!ACTION (1))
	     {
	       fprintf ( stderr, "+pdb2contacts <contact mode> <contact type> <D> -output contact_lib\n");
	       fprintf ( stderr, "+pdb2contacts def => <all contacts 1.2> => reports all contacts less than 1.2 A appart\n");
	       fprintf ( stderr, "    mode: intra       => Residues less than D appart with chain\n");
	       fprintf ( stderr, "    mode: inter       => Residues less than D appart between chains\n");
	       fprintf ( stderr, "    mode: all         => intra+inter\n\n");
	       
	       fprintf ( stderr, "    type: distances   => report CA distances between residue less than D appart\n");
	       fprintf ( stderr, "    type: contacts    => reports all AA with ATOMS less than D appart\n");
	       fprintf ( stderr, "    type: closest     => report only closest contact less than D appart\n");
	       fprintf ( stderr, "    type: count       => report the contact count\n");
	       fprintf ( stderr, "    type: best        => report the best contact for each AA pair\n");
	       fprintf ( stderr, "    type: find_pair   => RNA PDB\n");
	       fprintf ( stderr, "    type: find_pair-p => RNA PDB\n");
	       fprintf ( stderr, "    type: x3dna-ssr   => RNA PDB\n");
	       fprintf ( stderr, "    type: RNAplfold   => RNA sequence\\nn");
	       fprintf ( stderr, "    D   : distance in Angstrom (optional)\n");
	       
	       myexit (EXIT_SUCCESS);
	     }
	   else if (ACTION (1) && strm (ACTION(1), "def"))
	     {
	       sprintf (mode, "all");
	       sprintf (type, "contacts");
	       D=1.2;
	     }
	   else
	     {
	       sprintf (mode, "%s", ACTION(1));
	       sprintf (type, "%s", ACTION (2));
	       if (ACTION (3))D=atof (ACTION(3));
	     }
	   if (!T && !strm (type, "RNAplfold"))
	     {
	       printf_exit (EXIT_FAILURE,stderr, "+pdb2contacts %s %s requires a template:  -in2 <template_file> [FATAL]", mode, type);
	     }
	   
	   
	   ungap_seq(D1->S);
	   D1->CL=pdb2contacts (D1->S,T,D1->CL,mode,type,D);
	 }
       else if ( strm (action, "seq2contacts"))
	 {
	   //param1: mode
	   ungap_seq(D1->S); 
	   D1->CL=seq2contacts (D1->S, D2?(D2->S):NULL, D1->CL, ACTION (1));
	 }
       else if ( strm(action, "redundate"))
	 {
	   char *seq;
	   char *tree;

	   seq=(char*)vcalloc (100, sizeof (char));
	   tree=(char*)vcalloc (100, sizeof (char));
	   sprintf ( seq, "%s.redundated", (D1->S)->file[0]);
	   sprintf (tree, "%s.redundated", (D2->T)->file);

	   HERE ("%s %s", tree, seq);
	   redundate (D1->S, D2->T,seq, tree);
	 }
       else if ( strm(action, "treelist_prune")|| strm(action, "prune_treelist"))
	 {
	   Sequence *TS;
	   if (D2 && D2->S)TS=D2->S;
	   else TS=treelist2sub_seq((D1->S),ATOI_ACTION(1));
	   treelist2prune_treelist ( D1->S,TS, NULL);
	   D1->A=seq2aln (D1->S, NULL, NO_PAD);
	 }
       else if ( strm (action, "tree2unresolved_nodes"))
	 {
	   int ns;
	   int *l;
	   ns=tree2nseq (D1->T);
	   l=(int*)vcalloc (ns, sizeof (int));
	   tree2nnode_unresolved (D1->T, l);
	   for ( a=0; a<ns; a++)if (l[a])fprintf ( stdout, "SIZE: %d COUNT: %d\n", a, l[a]);
	   vfree (l);
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm(action, "tree_prune") || strm(action, "prune_tree"))
	 {
	   D1->T=prune_tree ( D1->T, D2->S);
	 }
       else if ( strm ( action, "tree2seq"))
	 {
	   D1->S=tree2seq(D1->T, NULL);
	   D1->A=seq2aln (D1->S, D1->A, 1);
	   (D1->A)->len_aln=1;
	   for ( a=0; a< (D1->A)->nseq; a++)sprintf ( (D1->A)->seq_al[a], "sequence");
	 }
       else if ( strm (action, "seq2dpatree"))
	 {
	   D1->T= seq2dpa_tree(D1->S,"ktup");
	 }
      
       else if ( strm (action, "tree2dpatree"))
	 {
	   D1->T= tree2dpa_tree(D1->T,(D2 && D2->A)?D2->A:D1->A, const_cast<char*>( (n_actions==1)?"idmat":action_list[1]) );
	 }
       else if ( strm (action, "tree2group"))
	 {
	   vfclose (tree2group (D1->T, (tree2seq(D1->T,NULL)), atoi(action_list[1]), atoi(action_list[2]),(n_actions==4)?action_list[3]:NULL, stdout));
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm(action, "unroot"))
	 {
	   D1->T=unroot_tree(D1->T);
	 }

       else if ( strm(action, "treelist2group")|| strm(action, "treelist2groups") )
	 {
	   Sequence *TS;

	   if (D2 && D2->S)TS=D2->S;
	   else TS=treelist2seq((D1->S));
	   treelist2groups (D1->S, TS, ACTION(1), stdout);
	   myexit (EXIT_SUCCESS);

	   //	   treelist2groups (D1->S,(D2)?D2->S:NULL, ACTION(1), stdout );
	   //exit (EXIT_SUCCESS);
	 }
       else if ( strm(action, "splits2tree"))
	  {

	   D1->T=split2tree ((D2)?D2->T:NULL,D1->S, ACTION(1));

	 }
       else if ( strm(action, "count_splits"))
	 {

	   count_splits ((D2)?D2->T:NULL,D1->S, ACTION(1));
	   myexit (EXIT_SUCCESS);
	 }
        else if ( strm(action, "count_groups"))
	 {
	   count_tree_groups (D1->S, ACTION(1));
	 }
       else if ( strm (action, "tree2dist"))
	 {
	   int ta, tb, ***td;
	   Sequence *TS;

	   TS=(D2)?D2->S:NULL;
	   td=tree2dist (D1->T,TS, NULL);
	   if (!TS)TS=tree2seq(D1->T, NULL);
	   for (ta=0; ta<TS->nseq; ta++)
	     {
	       fprintf ( stdout, "%-15s ",TS->name[ta]);
	       for ( tb=0; tb<TS->nseq; tb++)
		 {
		   int n=0;
		   if ( ACTION(1) && strm (ACTION(1), "length"))n=1;

		   fprintf (stdout, " %4d", td [n][ta][tb]);
		 }
	       fprintf ( stdout, "\n");
	     }
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "treelist2lti"))
	 {
	   Sequence *TS;
	   if (D2 && D2->S)TS=D2->S;
	   else TS=treelist2sub_seq((D1->S),ATOI_ACTION(2));
	   treelist2lti (D1->S,TS, (int)ATOI_ACTION(1), stdout );
	   myexit (0);
	 }
       else if ( strm (action,"treelist2frame"))
	 {
	   Sequence *TS;
	   if (D2 && D2->S)TS=D2->S;
	   else TS=treelist2sub_seq((D1->S),ATOI_ACTION(1));
	   treelist2frame (D1->S, TS);
	   myexit (EXIT_SUCCESS);
	 }

       else if ( strm (action, "treelist2seq"))
	 {
	   D1->S=treelist2sub_seq (D1->S,ATOI_ACTION(1));
	   D1->A=seq2aln(D1->S, NULL, 1);
	 }
       else if ( strm (action, "treelist2leafgroup"))
	 {
	   treelist2leafgroup (D1->S, (D2)?D2->S:NULL, ACTION(1));
	   myexit (0);
	 }
       else if ( strm(action, "treelist2splits"))
	 {
	   if (D1->T)D1->S=add_file2file_list ((D1->T)->file, NULL);
	   treelist2splits (D1->S, (D2)?D2->S:NULL);
	 }

       else if ( strm(action, "treelist2dmat"))
	 {
	   treelist2dmat (D1->S);
	 }
       else if ( strm(action, "tree2collapse") )
	 {
	   char *string;
	   int x,ng,l;
	   
	   if (!D1->T && (D1->A)->Tree)D1->T=newick_string2tree (((D1->A)->Tree)->seq_al[0]);
	   
	   l=0;
	   if ( strm (ACTION(1), "groups"))
	     {
	       ng=atoi(ACTION (2))+1;
	       l=1;
	       string=(char*)vcalloc ( 1000, sizeof (char));
	     }
	   else 
	     ng=n_actions;
	   
	   for (x=1; x<ng; x++)
	     {
	       if (!l)
		 {
		   string=ACTION(x);
		   string=substitute_char (string, '\\', 0);
		 }
	       else
		 {
		   sprintf (string, "-%d", x);
		 }
	       D1->T=collapse_tree (D1->T, NULL,string);
	     }
	   D1->S=tree2seq(D1->T, NULL);
	   D1->A=seq2aln (D1->S, NULL, RM_GAP);
	 }
      
       else if ( strm(action, "tree2node") )
	 {
	   print_node_list ( D1->T,(DST)?DST->S:NULL);
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm(action, "newick2randomize4dpa"))
	 {
	   // +newick2randomize4dpa dpa_nseq n_replicates
	   
	   int shuf=atoi (ACTION(1));
	   int n=atoi (ACTION (2));
	   newick2random4dpa(D2->S,D1->T, n, shuf); 
	   myexit(EXIT_SUCCESS);
	 }
       else if ( strm(action, "newick_randomize"))
	 {
	   int shuf=atoi (ACTION(1));
	   int zz;
	   for (zz=0; zz<shuf; zz++)
	     {
	       no_rec_print_tree_randomize (D1->T, stdout);
	       fprintf (stdout, ";\n");
	     }
	   myexit(EXIT_SUCCESS);
	 }
       else if ( strm(action, "newick_shuffle"))
	 {
	   int shuf=atoi (ACTION(1));
	   int zz;
	   for (zz=0; zz<shuf; zz++)
	     {
	       no_rec_print_tree_shuffle (D1->T, stdout);
	       fprintf (stdout, ";\n");
	     }
	   myexit(EXIT_SUCCESS);
	 }
	   
       else if ( strm(action, "tree_cmp_list") )
	 {
	   D1->T=main_compare_trees_list ( D1->T, D2->S, stdout);
	 }
       else if ( strm(action, "tree_cmp") || strm (action, "tree_compare"))
	 {
	   D1->T=main_compare_trees ( D1->T, D2->T, stdout);
	 }
       else if ( strm (action, "tree_scan"))
	 {
	   D1->T=tree_scan (D1->A, D2->T, ACTION(1), ACTION(2));
	 }
       else if ( strm (action, "split_cmp"))
	 {
	   main_compare_splits (D1->T, D2->T, ACTION(1), stdout);
	 }

       else if ( strm(action, "node_sort"))
	 {
	   node_sort ( action_list[1], D1->T);
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "mafftnewick2newick"))
	 {
	   if (!D2->S)printf_exit (EXIT_FAILURE,stderr,"ERROR action +mafftnewick2newick requires -in2=<seqfile> [FATAL]");
	   seqindex2seqname4tree (D1->T, D2->S);
	 }
       else if ( strm (action, "newick2mafftnewick"))
	 {
	   
	   if (!D2->S)printf_exit (EXIT_FAILURE,stderr,"ERROR action +newick2mafftnewick requires -in2=<seqfile> [FATAL]");
	   seqname2seqindexname4tree (D1->T, D2->S);
	   
	 }
       else if ( strm ( action, "treelist2bs") ||strm ( action, "tree2bs") )
	 {
	   if (ACTION(1) && strm (ACTION(1), "best"))treelist2node_support_best (D1->A);
	   else if (ACTION(1) && strm (ACTION(1), "cons"))treelist2cons (D1->A);
	   else treelist2node_support (D1->A);

	   if (!int_variable_isset ("print_replicates"))(D1->A)->nseq=1;
	 }
       else if (strm (action, "print_replicates"))
	 {
	   set_int_variable ("print_replicates",1);
	 }
       else if ( strm (action, "treelist2bs_compare"))
	 {
	   
	   treelist2bs_compare (D1->A, (D2)?D2->A:NULL);
	   exit (EXIT_SUCCESS);
	 }
       else if ( strm ( action, "tree2nni"))
	 {
	   tree2nni (D1->T, NULL);
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm ( action, "tree2ns"))
	 {
	   treelist2ns (D1->T, D2->S, ACTION(1));
	 }
	     
       else if ( strm ( action, "print"))
	 {
	   int x;
	   
	   for (x=1; x<n_actions;x++)
	     {
	       if ( strm (ACTION(x), "bs"))
		 {
		   float bs;
		   if (D1->T)bs=tree2avg_bs(D1->T);
		   else if (D1->A && (D1->A)->Tree)bs=newick2avg_bs (((D1->A)->Tree)->seq_al[0]);
		   else bs=newick2avg_bs ((D1->A)->seq_al[0]);
		   fprintf ( stdout, "AVERAGE_BS: %.2f\n", bs);
		 }
	       else if ( strm (ACTION(x), "nseq"))
		 {
		   fprintf (stdout, "NSEQ: %d\n", (D1->A)->nseq);
		 }
	     }
	 }
       else if ( strm (action, "genepred2acc"))
	 {
	   //D2->S=reference
	   //D1->S=prediction
	   vfree (display_accuracy (genepred_seq2accuracy_counts (D2->S, D1->S, NULL),stderr));
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "tree_cog_cmp"))
	 {
	   main_compare_cog_tree (D1->T,action_list[1]);
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "tree_aln_cmp"))
	 {
	   main_compare_aln_tree (D1->T, D2->A, stdout);
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm(action, "change_bootstrap"))
	 {
	   D1->T=reset_boot_tree ( D1->T, (n_actions>=2)?atoi(action_list[1]):0);
	 }
       else if ( strm(action, "change_distances"))
	 {
	   D1->T=reset_dist_tree ( D1->T, (n_actions>=2)?atof(action_list[1]):0.00);
	 }
       else if ( strm(action, "seq2dnd"))
	 {
	   D1->T=seq2dnd (D1->S, ACTION(1));
	 }
       else if ( strm(action, "aln2tree"))
	 {
	   D1->T=tree_compute (D1->A, n_actions-1, action_list+1);
	 }
       else if  ( strm(action, "aln2km_tree"))
	 {
	   if (!ACTION(1))myexit(fprintf_error (stderr,"-aln2km_tree <mode:diaa|triaa|aln> <nboot>"));
	   D1->T= aln2km_tree(D1->A, (ACTION(1)), ATOI_ACTION(2));
	 }
       else if ( strm(action, "similarities2tree"))
	 {
	   D1->T=similarities_file2tree (ACTION(1));
	 }

       else if (  strm(action, "original_seqnos"))
	 {
	  (D1->A)->output_res_num=2;
	 }
      	   
       else if ( strm (action, "aln2pred"))
	 {
	   aln2pred (D1->A, D2->A, ACTION (1));
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm(action, "color"))
	 {
	   cputenv ("TREE_MODE_4_TCOFFEE=%s",ACTION(1));
	 }
       else if ( strm(action, "gtree"))
	 {
	   gtree=1;
	   
	 }
       else if ( strm(action, "tree"))
	 {
	   if      (!ACTION(1))cputenv ("REPLICATES_4_TCOFFEE=1");
	   else if (is_number(ACTION(1)))cputenv ("REPLICATES_4_TCOFFEE=%s",ACTION(1));
	   else
	     {
	       int j;
	       for (j=1; j<n_actions; j+=2)
		 {
		   if      (strm (action_list[j], "mode")) cputenv ("TREE_MODE_4_TCOFFEE=%s",action_list[j+1]);
		   else if (strm (action_list[j], "gap" )) cputenv ("TREE_GAP_4_TCOFFEE=%s" ,action_list[j+1]);
		   else if (strm (action_list[j], "replicates"))cputenv ("REPLICATES_4_TCOFFEE=%s" ,action_list[j+1]);
		   else if (strm (action_list[j], "group"))cputenv ("SGROUP_4_TCOFFEE=%s" ,action_list[j+1]);
		   
		   else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: %s is not a known +tree parameter (replicates <int>|mode <string>|gap <float>|goup <seqfile>)[FATAL]",action_list[j]);
		 }
	     }
	   
	 }
       else if (strm(action, "columns4tree"))
	 {
	   char *f=ACTION(1);
	   set_string_variable ("columns4treeF",f);
	 }
       
       else if ( strm(action, "evaluateGroup"))
	 {
	   
	   if(ACTION(1))D1->A=evaluate_tree_group ((D1->A), main_read_seq(ACTION(1)));
	   else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: -evaluateGroup requires a sequence list in FASTA formal [FATAL]") ;
	 }
       else if ( strm(action, "evaluateTree"))
	 {
	   if (ACTION(1))
	     {
	       Sequence *G=main_read_seq (ACTION(1));
	       DST->A=treealn_evaluate4tcoffee (D1->A,G);
	     }
	   else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: -evaluateTree requires a sequence list in FASTA formal [FATAL]") ;
	 }
       
       else if ( strm(action, "evaluate3DM"))
	 {
	   int enb;
	   float max=0;
	   char *strikem=NULL;
	   Constraint_list *CL;
	   Alignment *A;
	   int na=1;
	   char *ev3d;

	   if (ACTION(na) && strm (ACTION(na), "group"))
	     {
	        cputenv ("SGROUP_4_TCOFFEE=%s" ,ACTION(na+1));
		cputenv ("REPLICATES_4_TCOFFEE=columns");
		na+=2;
	     }
	   if (strm (ACTION(na), "strike"))
	     {
	       //no param or: <max> <enb> <strikeMatrixFile>, values can be left to "def"
	       ev3d="strike";
	       enb=3;
	       max=1.2;
	       strikem=(char*)vcalloc (100, sizeof (char));
	       sprintf (strikem, "strike");
	       if (ACTION (na+1))
		 {
		   if (strm (ACTION(na+1), "def"));
		   else max=atof(ACTION(na+1));
		   
		   if (strm (ACTION(na+2), "def") );
		   else enb=atoi(ACTION(na+2));
		   		   
		   if (strm (ACTION(na+3), "def" ));
		   else sprintf (strikem, "%s", ACTION(na+3));
		 }
	     }
	   else if (strm (ACTION(na), "distances"))
	     {
	       ev3d="distances";
	       enb=3;
	       max=15;//Angstrom
	       
	       if (ACTION(na+1)&& !strm (ACTION(na+1), "def"))max=atof(ACTION(na+1));
	       if (ACTION(na+2)&& !strm (ACTION(na+2), "def"))enb=atoi(ACTION(na+2));
	     }
	   else if (strm (ACTION(na), "contacts"))
	     {
	       ev3d="contacts";
	       enb=3;
	       max=1.2;
		   
	       if (ACTION(na+1) && !strm (ACTION(na+1), "def"))max=atof(ACTION(na+1));
	       if (ACTION(na+2) && !strm (ACTION(na+2), "def"))enb=atoi(ACTION(na+2));
	     }
	   DST->A=msa_list2struc_evaluate4tcoffee (D1->S,D2->S,DST->S,ev3d,max,enb, strikem);
	   D1->A=(DST->A)->A;
	   (DST->A)->A=NULL;
	 }
       else if ( strm(action, "msa2distances"))
	 {
	   float radius;
	   float threshold;
	   int min;
	   
	   Constraint_list *CL;
	   Alignment *A;
	   int na=1;
	   DST->A=copy_aln (D1->A, NULL);
	   DST->S=aln2seq(DST->A);
	   	  
	   
	   	   
	   if (!D2)CL=D1->CL;
	   else if (!D2->CL)CL=D1->CL;
	   else CL=D2->CL;
	   if (!CL)//do a default contact based evaluation
	     {
	       
	       ungap_seq(D1->S);
	       D1->CL=pdb2contacts (D1->S, D2?(D2->S):NULL,D1->CL, "intra","distances",-1);
	       CL=D1->CL;
	     }
	   

	   if (ACTION(1))
	     {
	       radius=atof(ACTION(1));
	     }
	   else
	     {
	       radius=20;
	     }
	   if (ACTION(2))
	     {
	       threshold=atof(ACTION(2));
	     }
	   else 
	     threshold=0.5;
	   
	   if (ACTION(3))
	     {
	       min=atoi(ACTION(3));
	     }
	   else 
	     min=4;
	   
	   DST->A=msa2distances (D1->A,CL,radius,threshold, min);
	 }
       else if ( strm(action, "evaluate3D"))
	 {
	   int enb;
	   float max=0;
	   char *strikem=NULL;
	   Constraint_list *CL;
	   Alignment *A;
	   int na=1;
	   DST->A=copy_aln (D1->A, NULL);
	   DST->S=aln2seq(DST->A);
	   char *ev3d;

	   
	   if (ACTION(na) && strm (ACTION(na), "group"))
	     {
	        cputenv ("SGROUP_4_TCOFFEE=%s" ,ACTION(na+1));
		cputenv ("REPLICATES_4_TCOFFEE=columns");
		na+=2;
	     }
	   if (strm (ACTION(na), "strike"))
	     {
	       //no param or: <max> <enb> <strikeMatrixFile>, values can be feft to "def"
	       ev3d="strike";
	       enb=3;
	       max=1.2;
	       strikem=(char*)vcalloc (100, sizeof (char));
	       sprintf (strikem, "strike");
	       if (ACTION (na+1))
		 {
		   if (strm (ACTION(na+1), "def"));
		   else max=atof(ACTION(na+1));
		   
		   if (strm (ACTION(na+2), "def") );
		   else enb=atoi(ACTION(na+2));
		   		   
		   if (strm (ACTION(na+3), "def" ));
		   else sprintf (strikem, "%s", ACTION(na+3));
		 }
	     }
	   
	   else if (strm (ACTION(na), "distances"))
	     {
	       ev3d="distances";
	       enb=3;
	       max=150000;//Angstrom
	       
	       if (ACTION(na+1)&& !strm (ACTION(na+1), "def"))max=atof(ACTION(na+1));
	       if (ACTION(na+2)&& !strm (ACTION(na+2), "def"))enb=atoi(ACTION(na+2));
	       
	     }
	   else if (strm (ACTION(na), "contacts"))
	     {
	       ev3d="contacts";
	       enb=3;
	       max=1.2;
		   
	       if (ACTION(na+1) && !strm (ACTION(na+1), "def"))max=atof(ACTION(na+1));
	       if (ACTION(na+2) && !strm (ACTION(na+2), "def"))enb=atoi(ACTION(na+2));
	     }
	   
	   if (!D2)CL=D1->CL;
	   else if (!D2->CL)CL=D1->CL;
	   else CL=D2->CL;
	   
	   if (!CL)//do a default contact based evaluation
	     {
	       	       
	       ungap_seq(D1->S);
	       if (strm (ev3d, "distances"))
		 D1->CL=pdb2contacts (D1->S, D2?(D2->S):NULL,D1->CL, "intra","distances",max);
	       else if (strm (ev3d, "contacts"))
		 D1->CL=pdb2contacts (D1->S, D2?(D2->S):NULL,D1->CL, "intra","contacts",0);
	       else
		 {
		   
		   D1->CL=pdb2contacts (D1->S, D2?(D2->S):NULL,D1->CL,"intra", "contacts",0);
		   D1->CL=seq2contacts (D1->S, D2?(D2->S):NULL,D1->CL, NULL);
		 }
	       CL=D1->CL;
	     }


	   if (gtree)struc_evaluate4tcoffee4gt (D1->A,CL,ev3d,max,enb, strikem);
	   DST->A=struc_evaluate4tcoffee (D1->A,CL,ev3d,max,enb, strikem);
	 }
       else if ( strm(action, "hot"))
	 {
	   hot (D1->S, (D2)?D2->T:NULL, ACTION(1),0, ACTION(2),0);
	   exit(0);
	 }
       else if ( strm(action, "shot"))
	 {
	   hot (D1->S, (D2)?D2->T:NULL, ACTION(1),1, ACTION(2),0);
	   exit (0);
	 }
        else if ( strm(action, "hotshot"))
	 {
	   
	   float tot=0;
	   float nf=0;
	   hotshot (D1->S, (D2)?D2->T:NULL, ACTION(1),&tot, &nf);
	   fprintf ( stdout, "##HOTSHOT Global: %.3f\n",tot/nf); 
	   exit (0);
	 }
       else if ( strm(action, "sh"))
	 {
	   //method clustalw, clustalo, mafft
	   
	   float results=shuff (D1->S, ACTION (1),ATOI_ACTION(2));
	   exit (0);
	 }
	else if ( strm(action, "hotshot2"))
	 {
	   float *results=(float*)vcalloc (100, sizeof (float));
	   hotshot2 (D1->S, (D2)?D2->T:NULL, ACTION(1),results);
	   if (results[0]>0)
	     {
	       float t=results[0];
	       fprintf (stdout,"##GLOBAL HOT: %7.3f SHOT: %7.3f HOT-SHOT: %7.3f HOT>SHOT %7.3f SHOT>HOT %7.3f\n", results[2]/t,results[4]/t, results[6]/t, results[8]/t, results[10]/t);
	     }
	   
	   exit (0);
	 }
       else if ( strm(action, "evaluate"))
	 {
	   Alignment *A;
	   DST->A=copy_aln (D1->A, NULL);
	   DST->S=aln2seq(DST->A);
	   
	   
	   if    (strm (action_list[1], "id2"))
	     {
	       fprintf ( stdout, "ID2: %d\n", aln2sim2(D1->A));
	       exit (EXIT_SUCCESS);
	     }
	   else if (strm (action_list[1], "res2del") || strm (action_list[1], "del2ins"))
	     {
	       int s, c;
	       float threshold=0;
	       double*nr=(double*)vcalloc ((D1->A)->len_aln, sizeof (double));
	       double*ri=(double*)vcalloc ((D1->A)->nseq, sizeof (double));
	       double*gi=(double*)vcalloc ((D1->A)->nseq, sizeof (double));
	       Alignment *A=D1->A;
	       if (action_list[2])threshold=atof (action_list[2]);
	       double sum_gi, sum_gi2, avg_gi, sd_gi;
	       double sum_ri,sum_ri2, avg_ri, sd_ri;

	       sum_gi=sum_gi2=avg_gi=sd_gi=0;
	       sum_ri=sum_ri2=avg_ri=sd_ri=0;
	       
	       
	       for ( c=0; c< A->len_aln; c++)
		 {
		   for (s=0; s<A->nseq;s++)
		     {
		       if (A->seq_al[s][c]!='-')nr[c]++;
		     }
		 }
	       for ( c=0; c<A->len_aln; c++)
		 {
		   double r,g;
		   r=g=0;
		   if (nr[a]>0.1)    r=(A->nseq-nr[c])/nr[c];
		   if (nr[a]<A->nseq)g=(nr[c])/(A->nseq-nr[c]);
		   for (s=0; s<A->nseq; s++)
		     {
		       if (A->seq_al[s][c]!='-')ri[s]+=r;
		       else gi[s]+=g;
		     }
		 }
	       for (s=0; s<A->nseq; s++)
		 {
		   sum_gi +=gi[s];
		   sum_gi2+=gi[s]*gi[s];

		   sum_ri +=ri[s];
		   sum_ri2+=ri[s]*ri[s];
		 }
	       avg_gi=sum_gi/(double)A->nseq;
	       sd_gi=sqrt((sum_gi2-((sum_gi*sum_gi)/(double)A->nseq))/(double)A->nseq);

	       avg_ri=sum_ri/(double)A->nseq;
	       sd_ri=sqrt((sum_ri2-((sum_ri*sum_ri)/(double)A->nseq))/(double)A->nseq);
	       
	       
	      		   
	       for (s=0; s<A->nseq; s++)
		 {
		   double z;
		   if (strm (action_list[1], "res2del"))z=FABS((avg_ri-ri[s]))/sd_ri;
		   else z=FABS((avg_gi-gi[s]))/sd_gi;

		   if (threshold <0.00001 || z<(double)threshold)
		     {
		       fprintf ( stdout, ">%s del2ins: %.2f res2del: %.2f Z%s %.2f\n",A->name[s],gi[s],ri[s], action_list[1], (float)z);
		       fprintf ( stdout, "%s\n", A->seq_al[s]);
		     }
		 }
	       exit (EXIT_SUCCESS);
	     }
			   
	   else if    (strm (action_list[1], "id") || is_matrix(ACTION(1)))
	     {
	       Constraint_list *CL;
	       CL=(Constraint_list*)vcalloc (1, sizeof (Constraint_list));
	       if (is_matrix(ACTION(2)))CL->M=read_matrice (ACTION(2));
	       else if (is_matrix (ACTION(1)))CL->M=read_matrice (ACTION(1));
	       else CL->M=read_matrice ("idmat");


	       DST->A=matrix_evaluate_output (D1->A, CL);
	       
		    
	       
	       (D1->A)->score=(D1->A)->score_aln=(DST->A)->score=(DST->A)->score_aln;
	       
	       //fprintf ( stdout, "ID: %d\n", aln2sim ((D1->A), "idmat"));
	       //exit (EXIT_SUCCESS);
	       free_int (CL->M, -1);
	       vfree (CL);
	     }
	   else if (n_actions>1 && strm (  action_list[1], "categories"))
	     {
	       CL=declare_constraint_list ( DST->S,NULL, NULL, 0,NULL, read_matrice("pam250mt"));
	       DST->A=  main_coffee_evaluate_output(DST->A, CL, "categories");
	     }
	   else if (n_actions>1 && strm (  action_list[1], "sar"))
	     {
	       CL=declare_constraint_list ( DST->S,NULL, NULL, 0,NULL, read_matrice("pam250mt"));
	       DST->A=  main_coffee_evaluate_output(DST->A, CL, "sar");
	       (D1->A)->score=(D1->A)->score_aln=(DST->A)->score=(DST->A)->score_aln;
	       (DST->A)->score_seq[(D1->A)->nseq]*=10;
	     }
	   else if (n_actions>1 && strstr (  action_list[1], "boxshade"))
	     {
	       char color_mode[1000];
	       sprintf (color_mode,"boxshade_%d", atoi(ACTION2(2,"30")));
	       CL=declare_constraint_list ( DST->S,NULL, NULL, 0,NULL, read_matrice("pam250mt"));
	       DST->A=  main_coffee_evaluate_output(DST->A, CL, color_mode);
	     }
	  
	   else
	     {
	       printf_exit ( EXIT_FAILURE,stderr, "\nERROR: +evaluate mode [%s] is unknown[FATAL]", action_list[1]);
	     }
	   
	   DST->S=aln2seq ( DST->A);

	   A=D1->A;

	   //sprintf ( A->name[A->nseq], "cons");
	   //sprintf ( A->seq_al[A->nseq], "%s", aln2cons_seq_mat (A, "idmat"));

	 }
       else if ( strm (action, "sp_evaluate"))
	 {
	   fprintf ( stdout, "SP Score: %.2f", sum_pair ((DST && DST->A)?DST->A:D1->A,ACTION(1),atoi(ACTION2(2,"0")),atoi(ACTION2(3,"0"))));
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "lat_evaluate"))
	 {
	   float score;
	   score=lat_sum_pair ( D1->A, action_list[1]);
	   fprintf ( stdout, "\nLAT_SCORE: %.2f", score);
	   myexit (EXIT_SUCCESS);

	 }
       else if ( strm (action, "add_scale"))
	 {
	   D1->A=aln2scale (D1->A, ACTION(1));
	 }
       else if ( strm (action, "RNAfold_cmp"))
	 {
	   D1->A=compare_RNA_fold (D1->A, D2->A);
	 }
       else if ( strm (action, "aln2alifold"))
	 {
	   D1->A=aln2alifold (D1->A);
	   D1->S=aln2seq ( D1->A);
	 }


       else if ( strm (action, "add_alifold"))
	 {
	   D1->A=add_alifold2aln (D1->A, (D2)?D2->A:NULL);

	 }
       else if ( strm (action, "alifold2analyze"))
	 {
	   D1->A=alifold2analyze (D1->A, (D2)?D2->A:NULL, ACTION(1));
	   D1->S=aln2seq(D1->A);
	 }
       else if ( strm (action, "aln2conservation"))
	 {
	   D1->A=aln2conservation ( D1->A, ATOI_ACTION (1), ACTION (2));
	   myexit (EXIT_FAILURE);
	 }
       else if ( strm (action, "aln2cons"))
	 {
	   char *cons_seq;
	   char *cons_name;
	   cons_name=(char*)vcalloc (100, sizeof (char));
	   sprintf(cons_name, "%s", (n_actions<=2)?"Cons":action_list[2]);
	   cons_seq=aln2cons_seq_mat (D1->A, const_cast<char*>( (n_actions==1)?"blosum62mt":action_list[1]) );
	   free_aln (D1->A);free_sequence(D1->S, -1);
	   D1->S=fill_sequence_struc (1, &cons_seq, &cons_name, NULL);
	   /*keep the gaps*/
	   (D1->S)->len[0]=strlen (cons_seq); sprintf ( (D1->S)->seq[0], "%s", cons_seq);
	   D1->A=seq2aln (D1->S, NULL, KEEP_GAP);
	   vfree (cons_name);vfree (cons_seq);
	 }
       else if ( strm (action, "seq2filter"))
	 {
	   D1->S=seq2filter ( D1->S, atoi(action_list[1]), atoi(action_list[2]));

	 }
       else if ( strm (action, "aln2resindex"))
	 {
	   //-in: aln, file: ref_seq ref_res target_seq
	   //-in2 target sequences
	   aln2resindex (D1->A, (D2)?D2->A:NULL, stdout);
	   myexit (EXIT_SUCCESS);
	 }
       else if (strm(action, "keep_name"))
	 {
	   RAD->keep_name=1-RAD->keep_name;
	 }
        else if (strm(action, "use_consensus") ||strm(action, "use_cons") )
	 {
	   RAD->use_consensus=1-RAD->use_consensus;
	 }
       else if ( strm(action, "ungap"))
	 {
	   seq2aln (D1->S, D1->A, 1);
	 }
       else if ( strm2(action, "rmlower", "rm_lower"))
	 {
	   
	   RmLowerInAln(D1->A, ACTION(1));
	   D1->S=aln2seq ( D1->A);
	   (D1->A)->S=D1->S;
	 }
       else if ( strm2(action, "rmgap", "rm_gap"))
	 {
	   if (big())
	     {
	       sprintf(D1->file, "%s",ungap_fastaF_big(D1->file,NULL,(n_actions==1)?100:atoi(action_list[1])));
	       
	     }
	   ungap_aln_n (D1->A, (n_actions==1)?100:atoi(action_list[1]));
	   //free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	   (D1->A)->S=D1->S;
	 }
       else if ( strm(action, "rmgap_col"))
	 {
	   D1->A=remove_gap_column ( D1->A,action_list[1]);
	 }
       else if ( strm(action,"random"))
	 {

	   D1->A= make_random_aln(NULL,(n_actions==1)?1:atoi(action_list[1]),(n_actions==2)?100:atoi(action_list[2]),"acdefghiklmnpqrstvwy");

	   D1->S=aln2seq ( D1->A);
	 }

       else if ( strm(action, "landscape"))
	  {

	    set_landscape_msa ((n_actions==1)?0:atoi(action_list[1]));
	  }
       else if ( strm(action, "clean_maln"))
	  {
	    if ( !DST)
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL:%s]\n", PROGRAM);
		   myexit(EXIT_FAILURE);
		   }
	    (DST->A)=aln2number (DST->A);
	    D1->A=clean_maln(D1->A, DST->A,(n_actions==1)?1:atoi(action_list[1]),(n_actions==1)?1:atoi(action_list[2]));
	  }
       else if ( strm (action, "extract"))
	 {

	   COOR=get_pir_sequence  (RAD->coor_file, NULL);
	   D1->S=extract_sub_seq ( COOR, D1->S);
	   free_aln (D1->A);
	   D1->A=declare_Alignment(D1->S);
	   seq2aln (D1->S, D1->A, RAD->rm_gap);
	   free_sequence (COOR, COOR->nseq);
	 }
       else if ( strm (action, "reorder_column"))
	 {



	       Alignment *RO1, *RO2;
	       Sequence *OUT_S;
	       int s;

	       RO1=rotate_aln (D1->A,NULL);
	       if (ACTION(1) && strm (ACTION(1), "tree"))
		 {
		   D1->T=tree_compute (RO1,n_actions-2, action_list+2);
		    OUT_S=tree2seq(D1->T, NULL);
		    RO1=reorder_aln(RO1, OUT_S->name, OUT_S->nseq);
		  }
	       else if ( ACTION(1) && strm (ACTION(1), "random"))
		 {
		   RO1=reorder_aln ( RO1, NULL, RO1->nseq);
		 }

	       RO2=rotate_aln (RO1, NULL);
	       for (s=0; s< RO2->nseq; s++)
		 sprintf ( RO2->name[s], "%s", (D1->A)->name[s]);
	       free_aln (RO1);
	       free_aln (D1->A);
	       D1->A=RO2;
	       D1->S=aln2seq(D1->A);
	 }
       else if ( strm (action, "shuffle"))
	 {
	   D1->A=shuffle_aln (D1->A, ATOI_ACTION(1), ACTION(2), ACTION (3));
	 }
       else if ( strm (action, "reorder"))
	 {
	   
	   if ( n_actions==2 && strm (action_list[1], "random"))
	     {
	       D1->A=reorder_aln ( D1->A, NULL, (D1->A)->nseq);
	     }
	   else if (n_actions==2 && strm (action_list[1], "invert"))
	     {
	       char **nname;
	       int z, y;

	       nname=declare_char ((D1->A)->nseq, 100);
	       for ( z=0,y=(D1->A)->nseq-1; z<(D1->A)->nseq; z++, y--)
		 {
		   sprintf (nname[z], "%s",(D1->A)->name[y]);
		 }

	       D1->A=reorder_aln ( D1->A, nname, (D1->A)->nseq);
	       free_char (nname, -1);
	     }
	   else if (n_actions==2 && strm (action_list[1], "scramble"))
	     {
	       D1->A=aln2scramble_seq(D1->A);
	     }

	   else if ( n_actions==2 && strm (action_list[1], "tree"))
	     {

	       OUT_S=tree2seq (D2->T, NULL);
	       D1->A=reorder_aln(D1->A, OUT_S->name, OUT_S->nseq);
	       free_sequence (D1->S,(D1->S)->nseq);
	       D1->S=aln2seq (D1->A);
	     }
	   else
	     {
	       (D2->A)->S=aln2seq (D2->A);
	       (D1->A)->S=aln2seq (D1->A);
	       OUT_S=trim_aln_seq_name(D2->A, D1->A);
	       D1->A=reorder_aln(D1->A, OUT_S->name, OUT_S->nseq);
	       free_sequence (D1->S,(D1->S)->nseq);
	       D1->S=aln2seq (D1->A);
	     }
	 }
       else if ( strm (action, "aln2replicate"))
	 {
	   aln2N_replicate (D1->A, ACTION(1), ACTION(2));
	 }
       else if ( strm (action, "paralogous_cat"))
	 {
	   D1->A=orthologous_concatenate_aln (D1->A,D2->S, ACTION (1));
	 }

       else if ( strm (action, "cat_aln"))
	 {
	   /*D1->A=aln_cat ( D1->A, D2 ->A);*/

	   if (D2 && D2->A && !ACTION(1))
	     D1->A=concatenate_aln (D1->A, D2->A, ACTION(1));
	   else if (ACTION(1) && is_aln(ACTION(1)))
	     {
	         Alignment *B;
		 int n=1;

		 while (ACTION(n))
		   {

		     B=main_read_aln (ACTION(n), NULL);
		     D1->A=concatenate_aln (D1->A, B, NULL);
		     n++;
		   }
		 D1->S=aln2seq(D1->A);
	     }

	   else
	     {
	       Alignment *A, *B;

	       A=main_read_aln ((D1->A)->name[0], NULL);

	       for ( a=1; a<(D1->A)->nseq; a++)
		 {
		   B=main_read_aln ((D1->A)->name[a], NULL);
		   A=concatenate_aln (A, B, ACTION(1));

		 }
	       D1->A=A;
	       D1->S=aln2seq(D1->A);
	     }
	 }

       else if ( strm ( action, "msalist2cat_pwaln"))
	 {
	   int a, b, c;
	   int sim, min, max;

	   if (n_actions!=3)
	     {
	       min=0;
	       max=100;
	     }
	   else
	     {
	       min=atoi(action_list[1]);
	       max=atoi(action_list[2]);
	     }

	   fprintf ( stdout, ">A\n");
	   for (a=0;a<(D1->S)->nseq; a++)
	     {
	       Alignment *A;
	       HERE ("process %s",  (D1->S)->name[a]);
	       A=main_read_aln((D1->S)->name[a],NULL);
	       for (b=0; b<A->nseq-1; b++)
		 {
		   for ( c=b+1; c<A->nseq; c++)
		     {
		       sim=get_seq_sim (A->seq_al[b], A->seq_al[c], "-", "");
		       if (sim>=min && sim<=max)fprintf (stdout, "xxx%s", A->seq_al[b]);
		     }
		 }
	       free_aln (A);
	     }
	   fprintf ( stdout, "\n>B\n");
	   for (a=0;a<(D1->S)->nseq; a++)
	     {
	       Alignment *A;
	       HERE ("process %s",  (D1->S)->name[a]);
	       A=main_read_aln((D1->S)->name[a],NULL);
	       for (b=0; b<A->nseq-1; b++)
		 {
		   for ( c=b+1; c<A->nseq; c++)
		     {
		       sim=get_seq_sim (A->seq_al[b], A->seq_al[c], "-", "");
		       if (sim>=min && sim<=max)fprintf (stdout, "xxx%s", A->seq_al[c]);
		     }
		 }
	       free_aln (A);
	     }

	   fprintf ( stdout, "\n");
	   myexit (EXIT_SUCCESS);
	 }

       else if ( strm (action, "collapse_tree"))
	 {
	   D1->T=tree2collapsed_tree (D1->T, n_actions-1, action_list+1);
	 }
       else if ( strm (action, "collapse_aln"))
	 {
	   D1->A=aln2collapsed_aln (D1->A, n_actions-1, action_list+1);
	 }
       else if ( strm (action, "extract_aln"))
	 {
	   D1->A=aln2sub_aln_file (D1->A, n_actions-1, action_list+1);
	   myexit (EXIT_SUCCESS);
	 }



       else if ( strm (action, "remove_aa"))
	 {
	   int pos,len, n;
	   pos=atoi(action_list[1]);
	   len=atoi(action_list[2]);
	   n=atoi (action_list[3]);
	   if ( atoi (action_list[4])==1)len=-len;
	   if (pos && n>1)
	     {
	       fprintf ( stderr, "\nWARNING: rm_aa, position (pos) and iteration number (n) simulatneously defined. Iteration number reset to 1 [%s]\n", PROGRAM);
	       n=1;
	     }
	   for ( a=0; a< n; a++)
	     D1->A=probabilistic_rm_aa (D1->A, pos, len);
	 }
       else if ( strm (action, "remove_nuc"))
	 {
	   int pos;
	   pos=atoi(action_list[1]);

	   if ( pos>3 || pos<1)
	     printf_exit (EXIT_FAILURE, stderr, "Remove_nuc: indicate a number between 1 and 3\n");

	   pos--;
	   for ( c=0,a=0; a<(D1->A)->len_aln; a++, c++)
	     {
	       if (c==3)c=0;
	       for (b=0; b<(D1->A)->nseq; b++)
		 {
		 if (c==pos)
		   {
		     (D1->A)->seq_al[b][a]='-';
		   }
		 }
	     }

	   D1->S=aln2seq (D1->A);
	 }

       else if (strm ( action, "conserved_positions"))
	 {
	   Alignment *A;
	   int  a, b, c;
	   int *cache=NULL;


	   A=D1->A;
	   for ( a=0; a< A->nseq && !cache; a++)
	     {
	       if ( strm (action_list[1], A->name[a]))
		 {
		   cache=(int*)vcalloc ( A->len_aln+1, sizeof (int));
		   for ( c=0,b=0; b<A->len_aln; b++)
		     {
		       if ( is_gap (A->seq_al[a][b]))cache[b]=-1;
		       else cache[b]=++c;
		     }
		 }
	     }

	   for ( a=0; a< A->len_aln; a++)
	     {
	       r1=A->seq_al[0][a];
	       if ( is_gap(r1))continue;
	       for ( c=0,b=0; b<A->nseq; b++)
		 {
		   r2=A->seq_al[b][a];
		   c+=(r1==r2)?1:0;
		 }
	       if ( (c*100)/A->nseq>=atoi(action_list[2]))
		 fprintf ( stdout, "COL: %d Res: %c %s %d\n", a+1, r1, action_list[1], cache[a]+atoi(action_list[3]));
	     }
	   myexit (EXIT_FAILURE);
	 }
       else if (strm ( action, "extract_block") )
	 {

	   BUF=copy_aln (D1->A, NULL);
	   
	   if ( check_file_exists(action_list[1]))
	     BUF=extract_aln3(BUF,action_list[1]);
	   else
	     {
	       
	     BUF=extract_aln2(BUF,atoi(action_list[2]),atoi(action_list[3]),action_list[1]);
	     }
	   D1->A=copy_aln (BUF,D1->A);

	 }
       else if ( strm ( action, "extract_pos_list"))
	 {
	   D1->A=alnpos_list2block (D1->A, n_actions-1, action_list+1);
	 }
       else if ( strm ( action, "seq2msa"))
	 {
	   D1->A=simple_progressive_aln ( D1->S, NULL, NULL, action_list[1]);
	 }
       else if ( strm ( action, "realign_block") )
	 {
	   D1->A=realign_block ( D1->A, atoi (action_list[1]), atoi (action_list[2]), (n_actions==4)?action_list[3]:NULL);
	 }
       else if ( strm (action, "extract_seq"))
	 {
	   if (big())
	     {
	       char *tmp1;
	       D1->file,"%s", trim_fastaF_big (D1->file,action_list[1],tmp1=vtmpnam (NULL), NULL, NULL,NULL);
	       sprintf (D1->file, "%s",tmp1);
	     }
	   else
	     {
	       int is_file;
	       if ( check_file_exists (action_list[1])&& format_is_fasta (action_list[1]))
		 {
		   is_file=1;
		   BUFS=main_read_seq (action_list[1]);
		   action_list=BUFS->name;
		   n_actions=BUFS->nseq;
		 }
	       else
		 {
		   is_file=0;
		   action_list++;
		   n_actions--;
		 }
	       
	       for ( a=0; a< n_actions;)
		 {
			 s=action_list[a];

			 if ( n_actions==1 || is_file==1)
			 {
				 start=1;
				 end=0;
				 a+=1;
			 }
			 else
			 {

				 start=(strm2 (s,"#","*"))?1:(atoi(action_list[a+1]));
				 end=  (strm2 (action_list[a+2],"#","*"))?0:(atoi(action_list[a+2]));
				 a+=3;
			 }

			 if ( strm2 (s, "#", "*"))
			 {
				 OUT_S=extract_one_seq((D1->A)->name[0],start, end, D1->A, RAD->keep_name);
				 for (b=1; b< (D1->A)->nseq; b++)
				 {
					 NS=extract_one_seq((D1->A)->name[b],start, end, D1->A, RAD->keep_name);
					 if (count_n_res_in_array(NS->seq[0], -1))
						 OUT_S=add_sequence ( NS,OUT_S, 0);
				 }
			 }
			 else
			 {
				 if ( a==1)OUT_S=extract_one_seq(s,start, end, D1->A, RAD->keep_name);
				 else
				 {
					 NS=extract_one_seq(s,start, end, D1->A, RAD->keep_name);
					 OUT_S=add_sequence ( NS,OUT_S, 0);
				 }
			 }
		 }
		 D1->S=OUT_S;
		 free_aln (D1->A);
		 D1->A=declare_Alignment(D1->S);
		 seq2aln (D1->S, D1->A, RAD->rm_gap);
	     }
	 }
	 else if ( strm (action, "ls_extract_seq"))
	 {
		 // if given a file, read it
		 if ( check_file_exists (action_list[1]) && format_is_fasta (action_list[1]))
		 {

			 BUFS=main_read_seq (action_list[1]);
			 action_list=BUFS->name;
			 n_actions=BUFS->nseq;
		 }
		 else
		 {
			 action_list++;
			 n_actions--;
		 }

		//sort names for binary search
		char **seq_found;
		qsort(action_list, n_actions, sizeof(char*), ls_compare);
		char *names;
		size_t i,pos = 0;
		size_t n_seqs = D1->S->nseq;
		int max_len=0;
		int min_len=INT_MAX;
		for (i=0; i<n_seqs; ++i)
		{
			seq_found = (char**)bsearch(&D1->S->name[i], action_list, n_actions, sizeof(char*), ls_compare);
			if (seq_found != NULL)
			{	//copy values of a sequence to unused position
				D1->S->name[pos]=D1->S->name[i];
				D1->S->seq[pos]=D1->S->seq[i];
				D1->S->len[pos]=D1->S->len[i];
				if (max_len < D1->S->len[pos])
					max_len=D1->S->len[pos];
				if (min_len > D1->S->len[pos])
					min_len=D1->S->len[pos];
				D1->S->seq_comment[pos]=D1->S->seq_comment[i];
				D1->S->file[pos]=D1->S->file[i];
				D1->S->T[pos]=D1->S->T[i];
				if (D1->S->genome_co != NULL)
					D1->S->genome_co[pos]=D1->S->genome_co[i];
				++pos;
			}
			else
			{	//free memory of deleted sequences
				vfree(D1->S->name[i]);
				vfree(D1->S->seq[i]);
				vfree(D1->S->seq_comment[i]);
				vfree(D1->S->T[i]);
				vfree(D1->S->file[i]);
			}
		}

		//update values


		D1->S->max_nseq=pos;
		D1->S->nseq=pos;
		D1->S->max_len=max_len;
		D1->S->min_len=min_len;

		//free memory
		D1->S->name=(char**)vrealloc(D1->S->name, pos*sizeof(char*));
		D1->S->seq=(char**)vrealloc(D1->S->seq, pos*sizeof(char*));
		D1->S->seq_comment=(char**)vrealloc(D1->S->seq_comment, pos*sizeof(char*));
		D1->S->file=(char**)vrealloc(D1->S->file, pos*sizeof(char*));
		D1->S->T=(Template**)vrealloc(D1->S->T, pos*sizeof(Template*));
		D1->S->len=(int*)vrealloc(D1->S->len, pos*sizeof(int));
		if (D1->S->genome_co != NULL)
			vrealloc(D1->S->genome_co, pos*sizeof(Genomic_info));

		D1->A=declare_Alignment(D1->S);
		seq2aln (D1->S, D1->A, RAD->rm_gap);
	 }
	 else if ( strm (action, "extract_lib_list"))
	   {
	     D1->CL=constraint_list2sub_constraint_list (D1->CL,main_read_seq (action_list[1]));
	   }
       else if ( strm (action, "extract_seq_list"))
	 {
	   
	   if (D1->CL)
	      D1->CL=constraint_list2sub_constraint_list (D1->CL,main_read_seq (action_list[1]));
	   else
	     {
	       int nadded=0;
	       if ( check_file_exists (action_list[1]) && format_is_fasta (action_list[1]))
		 {
		   
		   BUFS=main_read_seq (action_list[1]);
		   action_list=BUFS->name;
		   n_actions=BUFS->nseq;
		 }
	       else
		 {
		   action_list++;
		   n_actions--;
		 }
	       
	       for ( a=0; a< n_actions;a++)
		 {
		   if ( (name_is_in_list (action_list[a], (D1->S)->name, (D1->S)->nseq, 100))!=-1)
		     {
		       nadded++;
		       HERE ("%s", action_list[a]);
		       NS=extract_one_seq(action_list[a],1,0, D1->A, KEEP_NAME);
		       OUT_S=add_sequence ( NS,OUT_S, 0);
		     }
		   else 
		     {
		       fprintf (stderr, "WARNING: %s could not be extracted\n", action_list[a]);
		     }
		 }
	       if (nadded==0)exit (0);
	       D1->S=OUT_S;
	       free_aln (D1->A);
	       D1->A=declare_Alignment(D1->S);
	       seq2aln (D1->S, D1->A, RAD->rm_gap);
	     }
	 }
       else if ( strm (action, "remove_seq") || strm (action, "rm_seq"))
	 {
	   char *buf=NULL;
	   char **list;
	   int n;
	   int l;

	   list=declare_char ((D1->S)->nseq, 200);
	   for ( n=0,a=0; a< (D1->A)->nseq; a++)
	     {

	       buf=csprintf(buf, "%s",(D1->S)->seq[a]);
	       ungap (buf);
	       l=strlen(buf);

	       for (c=1, b=1; b< n_actions; b++)
		 {
		   if ( strm (action_list[b], (D1->S)->name[a])){(D1->S)->seq[a]=NULL;break;}
		   else if ( strm (action_list[b], "empty") && l==0)
		     {
		       fprintf ( stderr, "WARNING: Sequence %s does not contain any residue: automatically removed from the set [WARNING:%s]\n",(D1->S)->name[a], PROGRAM);
		       (D1->S)->seq[a]=NULL;break;
		     }
		   else if ( strm (action_list[b], "unique"))
		     {
		       if ( name_is_in_list ((D1->S)->name[a], list,n, 100)!=-1)
			 {
			   (D1->S)->seq[a]=NULL;break;
			 }
		       else
			 {
			   sprintf ( list[n++], "%s", (D1->S)->name[a]);
			 }
		     }
		 }
	     }
	   D1->S=duplicate_sequence (D1->S);
	   free_aln (D1->A);
	   free_char ( list, -1);
	   D1->A=declare_Alignment(D1->S);
	   seq2aln (D1->S, D1->A, RAD->rm_gap);
	 }

       else if (  strm (action, "aln2overaln")|| strm (action,"overaln_param"))
	 {
	   //mode (lower|number|uanlign) Penalty (0-100) Thresold (0-9)
	   int  p1,p2,p3,f, t;
	   char *s;
	   int eb=0;
	   char clean_mode[100];
	   OveralnP *F;

	   F=(OveralnP*)vcalloc (1, sizeof (OveralnP));
	   if ( D2 && D2->A)
	     {
	       D1->A=mark_exon_boundaries (D1->A, D2->A);
	       eb=1;
	     }
	   else if ( (s=get_string_variable ("exon_boundaries")))
	     {
	      Sequence *S;
	      Alignment *EB;
	      EB=seq2aln(S=main_read_seq(s),NULL, 0);
	      D1->A=mark_exon_boundaries (D1->A, EB);
	      free_sequence (S, S->nseq); free_aln (EB);
	      eb=1;
	     }


	   if (ACTION(1)==NULL)sprintf (F->mode, "lower");
	   else if (strstr (ACTION(1), "h"))
	     {
	       fprintf ( stdout, "aln2unalign lower|number|unalign|uanlign2 F P1 P2 P3 T\n");
	       myexit (EXIT_SUCCESS);
	     }
	   else sprintf (F->mode, "%s", ACTION(1));

	   F->t=ATOI_ACTION(2);
	   F->f=ATOI_ACTION(3);
	   F->p1=ATOI_ACTION(4);
	   F->p2=ATOI_ACTION(5);
	   F->p3=ATOI_ACTION(6);
	   F->p3=ATOI_ACTION(7);

	   if (int_variable_isset ("overaln_target"))f=get_int_variable ("overaln_target");
	   if (int_variable_isset ("overaln_threshold"))t=get_int_variable ("overaln_threshold");
	   if (eb)sprintf (F->model, "fsa2");
	   else   sprintf (F->model, "fsa1");

	   D1->A=aln2clean_pw_aln (D1->A, F);

	 }
       else if (  strm (action, "unalign_groups"))
	 {
	   //unalign everything in lower case
	   unalign_aln_2 (D1->A, NULL, 0);
	 }
       else if (  strm (action,"aln2unalign"))
	 {
	   Alignment *SA;
	   Sequence *SS;
	   SA=copy_aln (D1->A, NULL);
	   SS=aln2seq(SA);

	   thread_seq_struc2aln (SA, SS);
	   D1->A=unalign_aln (D1->A,SA, ATOI_ACTION(1));
	   D1->S=aln2seq ( D1->A);
	 }
       else if (  strm (action, "clean_cdna"))
	 {
	   Alignment *A;
	   A=D1->A;
	   for (a=0; a< A->nseq; a++)
	     {
	       char *d, *buf, f;

	       d=A->seq_al[a];
	       f=get_longest_frame (d, 3);
	       buf=(char*)vcalloc ( strlen (d)+1, sizeof (char));
	       sprintf (buf, "%s", d+f);
	       sprintf (d, "%s", buf);
	       vfree (buf);
	     }
	 }
       else if ( strm (action, "clean_cdna2"))
	 {
	   D1->A=clean_cdna_aln ( D1->A);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if ( strm  (action, "aln2short_aln"))
	   {
	     D1->A=aln2short_aln (D1->A, action_list[1], action_list[2], atoi(action_list[3]));
	     free_sequence ( D1->S, (D1->S)->nseq);
	     D1->S=aln2seq ( D1->A);
	   }
       else if ( strm ( action, "complement"))
	 {
	   D1->A=complement_aln (D1->A);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if ( strm ( action, "extend"))
	 {
	   extend_seqaln( NULL,D1->A);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if ( strm ( action, "unextend"))
	 {
	   unextend_seqaln( NULL,D1->A);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if ( strm ( action, "translate"))
	 {
	   D1->A=translate_dna_aln( D1->A,(n_actions==1)?0:atoi(action_list[1]));
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if (strm2 ( action, "back_translate","backtranslate"))
	 {
	  D1->A=back_translate_dna_aln( D1->A);
	  free_sequence ( D1->S, (D1->S)->nseq);
	  D1->S=aln2seq ( D1->A);
	 }
       else if (strm ( action, "rotate"))
	 {
	   D1->A=rotate_aln( D1->A, action_list[1]);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if (strm ( action, "invert"))
	 {
	  D1->A=invert_aln( D1->A);
	  free_sequence ( D1->S, (D1->S)->nseq);
	  D1->S=aln2seq ( D1->A);
	 }
       else if (strm ( action, "test_dna2gene"))
	 {
	   testdna2gene ((D1->S)->seq[0]);
	 }
       else if (strm ( action, "code_dna_aln"))
	 {
	  D1->A=code_dna_aln( D1->A);
	  free_sequence ( D1->S, (D1->S)->nseq);
	  D1->S=aln2seq ( D1->A);
	 }

       else if ( strm ( action, "mutate"))
	 {
	   D1->A=mutate_aln( D1->A, const_cast<char*>( (n_actions==1)?"0":action_list[1]) );
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       else if ( strm ( action, "thread_profile_on_msa"))
	 {
	   (D1->A)->S=NULL;
	   D1->A=thread_profile_files2aln (D1->A, action_list[1], NULL);
	   D1->S=aln2seq(D1->A);
	 }
       else if ( strm ( action, "thread_dna_on_prot_aln"))
	 {
	   if (D1->S)fast_get_sequence_type(D1->S);
	   if (D2->S)fast_get_sequence_type(D2->S);
	   
	   if (D1->S && strm ((D1->S)->type, "DNA"))
	     D1->A=thread_dnaseq_on_prot_aln (D1->S, D2->A);
	   else if (D2->S && strm ((D2->S)->type, "DNA"))
	     D1->A=thread_dnaseq_on_prot_aln (D2->S, D1->A);
	   else
	     printf_exit (EXIT_FAILURE, stderr, "Error: +thread_dna_on_prot_aln requires -in=<prot_aln> -in2=<dna sequence>");
	   
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       
       //else if ( strm ( action, "thread_dna_on_prot_aln"))
       //{
       //   D1->A=thread_dnaseq_on_prot_aln (D1->S, D2->A);
       //   free_sequence (D1->S,(D1->S)->nseq);
       //   D1->S=aln2seq (D1->A);
       // }
       else if ( strm ( action, "thread_struc_on_aln"))
	 {
	   thread_seq_struc2aln ( D2->A, D1->S);
	   D1->A=copy_aln(D2->A, NULL);
	   D1->S=aln2seq (D1->A);
	 }
      
       else if ( strm (action, "sim_filter"))
	 {
	   D1->A=sim_filter (D1->A, action_list[1], ACTION (2));
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       else if ( strm (action, "protein_db")||  strm (action, "db"))
	 {
	   cputenv ("protein_db_4_TCOFFEE=%s", action_list[1]);
	 }
       else if ( strm (action, "compress"))
	 {
	   cputenv ("compress_4_TCOFFEE=1");
	   

	 }
       else if ( strm (action, "thread"))
	 {
	   int nproc=(atoi(action_list[1]));
	   if (!nproc)nproc=get_nproc();
	   cputenv ("thread_4_TCOFFEE=%d", nproc);
	 }
       else if ( strm (action, "psiJ") || strm (action, "num_irtrations"))
	 {
	   cputenv ("num_iterations_4_TCOFFEE=%s", action_list[1]);
	 }
       
       else if ( strm (action, "outfmt"))
	 {
	   cputenv ("outfmt_4_TCOFFEE=%s", action_list[1]);
	 }
       else if ( strm (action, "outdir"))
	 {
	   cputenv ("cache_4_TCOFFEE=%s", action_list[1]);
	 }
       else if ( strm (action, "cache"))
	 {
	   cputenv ("cache_4_TCOFFEE=%s", action_list[1]);
	 }
       else if ( strm (action, "seq2blast"))
	 {
	   int a;
	   if (!getenv("thread_4_TCOFFEE"))cputenv ("thread_4_TCOFFEE=1");
	   D1->S=seq2blast (D1->S);
	   if (D1->A)
	     {
	       for ( a=0; a<(D1->S)->nseq; a++)
		 {
		   (D1->A)->seq_comment[a] =csprintf ((D1->A)->seq_comment[a], "%s",(D1->S)->seq_comment[a]);
		   (D1->A)->aln_comment[a] =csprintf ((D1->A)->aln_comment[a], "%s",(D1->S)->aln_comment[a]);
		 }
	     }
	 }
       else if ( strm (action, "prot_min_cov"))
	 {
	   cputenv ("prot_min_cov_4_TCOFFEE=%s", action_list[1]);
	 }
        else if ( strm (action, "prot_min_sim"))
	 {
	   cputenv ("prot_min_sim_4_TCOFFEE=%s",action_list[1]);
	 }
	else if ( strm (action, "prot_max_sim"))
	 {
	   cputenv ("prot_max_sim_4_TCOFFEE=%s",action_list[1]);
	 }
        else if ( strm (action, "psitrim"))
	 {
	   cputenv ("psitrim_4_TCOFFEE=%s",action_list[1]);
	 }
       else if ( strm (action, "psitrim_tree"))
	 {
	   cputenv ("psitrim_tree_4_TCOFFEE=%s",action_list[1]);
	 }
       else if ( strm (action, "psitrim_mode"))
	 {
	   cputenv ("psitrim_mode_4_TCOFFEE=%s",action_list[1]);
	 }
       
       else if ( strm (action, "seq2prf"))
	 {
	   D1->S=seq2prf (D1->S);
	   if (D1->A)
	     {
	       for ( a=0; a<(D1->S)->nseq; a++)
		 {
		   (D1->A)->seq_comment[a] =csprintf ((D1->A)->seq_comment[a], "%s",(D1->S)->seq_comment[a]);
		   (D1->A)->aln_comment[a] =csprintf ((D1->A)->aln_comment[a], "%s",(D1->S)->aln_comment[a]);
		 }
	     }
	 }
       else if ( strm (action, "kmeans"))
	 {
	   //k, mode: diaa,triaa, name
	   if (!ATOI_ACTION(1))myexit(fprintf_error (stderr,"-kmeans <nclusters> <mode:diaa|triaa> <outputname>"));
	   km_seq (D1->A,ATOI_ACTION (1),ACTION(2),ACTION(3));
	 }
       else if ( strm (action, "gap_trim"))
	 {
	   D1->A=gap_trim (D1->A,ATOI_ACTION(1));

	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       else if ( strm (action, "trimRNA"))
	 {
	   D1->A=trim_RNA(D1->A, D2->S, ATOI_ACTION(1));
	 }
       else if ( strm (action, "phylotrim"))
	 {
	   
	   D1->A=phylotrim (D1->A,(D2)?D2->T:NULL, ACTION(1), ACTION (2), ACTION (3));

	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       
       else if ( strm (action, "regtrim"))
	 {
	   int ns, p;
	   char *tt=get_string_variable ("treemode");

	  
	   if (strstr (action_list[1], "%"))
	     {
	       sscanf (action_list[1], "%d%%", &p);
	       ns=(D1->A)->nseq*((float)p/(float)100);
	     }
	   else
	     {
	       ns=atoi (action_list[1]);
	     }
	   if (!D2 && strstr(tt, "msa"))
	     {
	       (D1->S)->seq=(D1->A)->seq_al;
	     }
	   
	   D1->S=regtrim (D1->S, (D2)?D2->T:NULL,ns);
	   D1->A=seq2aln(D1->S, NULL, 0);
	 }
       else if ( strm (action, "trim"))
	 {
	   D1->A=simple_trimseq (D1->A,(D2)?D2->A:NULL, action_list[1], ACTION (2), NULL);
	   
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	   
	 }

       else if (strm ( action, "trimTC"))
	 {
	   value=(n_actions==1)?10:atoi(action_list[1]);

	   D1->A=tc_trimseq(D1->A,D1->S,action_list[1]);
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       else if (strm ( action, "trimTC2"))
	 {
	   char *group_file;
	   Alignment *B=NULL;
	   char trim_mode[100];
	   if ( n_actions==1 || !(strm (action_list[1], "NSEQ") ||strm (action_list[1], "MINID")) )
	     {
	       fprintf ( stderr, "\nTrimTC2 <NSEQ | MINID>  <number sequences| minimum identity> (<matrix>)\n");
	       myexit (EXIT_FAILURE);
	     }
	   sprintf (trim_mode, "%s", action_list[1]);action_list+=2; n_actions-=2;

	   if ( strm ( trim_mode, "NSEQ"))
	     {
	       group_file=tree2Ngroup( (D1)?D1->A:NULL, (D2)?D2->T:NULL, atoi (action_list[0]), vtmpnam(NULL), const_cast<char*>( (n_actions==1)?"idmat":action_list[1]) );
	     }
	   else
	     {
	       group_file=tree2Ngroup( (D1)?D1->A:NULL, (D2)?D2->T:NULL, -1*atoi (action_list[0]), vtmpnam(NULL), const_cast<char*>( (n_actions==1)?"idmat":action_list[1]) );
	     }

	   B=copy_aln (D1->A, B);
	   B=aln2sub_aln_file (B,1,&group_file);
	   B=aln2sub_seq (B, 1, &group_file);
	   D1->A=extract_sub_aln2 (D1->A, B->nseq, B->name);
	 }
       else if ( strm (action, "chain"))
	 {
	   D1->A=seq2seq_chain (D1->A,D2->A, ACTION(2));
	 }


       else if (strm ( action, "master_trim"))
	 {
	   value=(n_actions==1)?10:atoi(action_list[1]);

	   D1->A=master_trimseq(D1->A,D1->S,action_list[1]);
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
        else if ( strm (action, "force_aln"))
	 {
	   char ***rlist=NULL;
	   int count=0;

	   if ( n_actions==2)
	     {
	       if (!is_lib_02(action_list[1]))
		 {
		   fprintf ( stderr, "\nERROR: force_aln requires files in TC_LIB_FORMAT_02 [FATAL:%s]", PROGRAM);
		   myexit (EXIT_FAILURE);
		 }
	       else
		   rlist=file2list (action_list[1], " ");
	     }
	   else
	     {
	       rlist=(char***)declare_arrayN(3, sizeof (char),3,7, 10);

	       strcat (rlist[1][1],action_list[1]);strcat (rlist[1][3],action_list[2]);
	       strcat (rlist[1][4],action_list[3]);strcat (rlist[1][6],action_list[4]);
	       sprintf ( rlist[2][0], "-1");
	     }
	   count=1;
	   while (rlist[count] && atoi(rlist[count][0])!=-1)
	     {
	       char st1[100], st2[100], st3[100], st4[100];

	       sprintf ( st1, "%s", rlist[count][1]);sprintf ( st2, "%s", rlist[count][3]);
	       sprintf ( st3, "%s", rlist[count][4]);sprintf ( st4, "%s", rlist[count][6]);
	       fprintf ( stderr, "\nFORCE: %s %s %s %s", st1, st2, st3, st4);

	       if (is_number (st1))s1=atoi (st1)-1;
	       else s1=name_is_in_list (st1,(D1->A)->name, (D1->A)->nseq, 100);
	       if ( s1<0 || s1>= (D1->A)->nseq)crash ("wrong sequence index");
	       r1=atoi (st2)-1;

	       if (is_number (st3))s2=atoi (st3)-1;
	       else s2=name_is_in_list (st3,(D1->A)->name, (D1->A)->nseq, 100);
	       if ( s2<0 || s2>= (D1->A)->nseq)crash ("wrong sequence index");
	       r2=atoi (st4)-1;

	       (D1->A)=add_constraint2aln ((D1->A), s1, r1, s2, r2);
	       count++;
	     }
	   fprintf ( stderr, "\n");
	   free_arrayN((void*)rlist,3);
	 }

        else if (strm ( action, "grep"))
	  {
	    D1->A=grep_seq (D1->A, ACTION(1),ACTION(2), ACTION(3));
	    if (D1->A==NULL) myexit (EXIT_SUCCESS);
	    else D1->S=aln2seq (D1->A);
	  }

	else if (strm (action, "find"))
	  {
	    int r, l;
	    char *search_string;

	    search_string=(char*)vcalloc ( 30, sizeof (char));
	    if ( strm (action_list[1], "lower"))sprintf ( search_string, "abcdefghijklmnopqrstuvwxyz");
	    else if ( strm ( action_list[1], "upper"))sprintf ( search_string, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	    else
	      {
		vfree (search_string);search_string=(char*)vcalloc ( strlen (action_list[1])+1, sizeof (char));
		sprintf (search_string, "%s", action_list[1]);
	      }

	    for (a=0; a<(D1->A)->nseq; a++)
	      for ( l=0,b=0; b< (D1->A)->len_aln; b++)
		{
		  r=(D1->A)->seq_al[a][b];
		  l+=!is_gap(r);
		  if ( r!='\0' && strrchr (search_string, r))
		    {
		      /*fprintf ( stdout, "%-15s res %c alnpos %4d seqpos %4d\n", (D1->A)->name[a], r, b+1, l);*/
		      fprintf ( stdout, "%s %d %d\n", (D1->A)->name[a], l, l+1);
		    }
		}
	    myexit (EXIT_SUCCESS);
	  }
        else if ( strm (action, "merge_annotation"))
	  {
	    D1->A=merge_annotation (D1->A, DST?DST->A:NULL, ACTION(1));
	    D1->S=aln2seq (D1->A);
	  }
	else if ( strm  (action, "color_residue"))
	  {
	    int i;
	    Alignment *A;
	    A=D1->A;

	    DST->A=copy_aln (D1->A, NULL);
	    DST->S=aln2seq (DST->A);
	    for (a=0; a< (DST->S)->nseq; a++)ungap ((DST->S)->seq[a]);

	    if (n_actions>2)
	      {
		for (a=1; a<n_actions; a+=3)
		  {
		    i=name_is_in_list(action_list[a], (D1->A)->name, (D1->A)->nseq, 100);
		    if (i!=-1)
		      {
			(DST->S)->seq[i][atoi(action_list[a+1])-1]='0'+atoi(action_list[a+2])-1;
		      }
		    else fprintf (stderr, "\nWARNING: Could not find Sequence %s", action_list[a]);
		  }
	      }
	    else
	      {
		char name[1000];
		int pos, val;
		FILE *fp;

		fp=vfopen (action_list[1], "r");
		while (fscanf (fp, "%s %d %d\n", name, &pos, &val)==3)
		  {

		     i=name_is_in_list(name, (D1->A)->name, (D1->A)->nseq, 100);
		     if (i!=-1)(DST->S)->seq[i][pos-1]='0'+val;
		     else fprintf (stderr, "\nWARNING: Could not find Sequence %s", action_list[a]);
		  }
		vfclose (fp);
	      }
	    DST->A=seq2aln (DST->S, NULL, 1);
	  }
       else if ( strm  (action, "edit_residue"))
	  {
	    FILE *fp;
	    int i, pos;
	    int **p;
	    char mod[100], name[100];
	    Alignment *A;

	    A=D1->A;

	    p=aln2inv_pos (A);
	    if (n_actions>2)
	      {
		for (a=1; a<n_actions; a+=3)
		  {

		    i=name_is_in_list(action_list[a], (D1->A)->name, (D1->A)->nseq, 100);
		    if (i!=-1)
		      {
			pos=atoi(action_list[a+1]);

			pos=p[i][pos]-1;
			sprintf (mod, "%s", action_list[a+2]);
			if ( strm (mod, "upper"))(D1->A)->seq_al[i][pos]=toupper((D1->A)->seq_al[i][pos]);
			else if ( strm (mod, "lower"))(D1->A)->seq_al[i][pos]=tolower((D1->A)->seq_al[i][pos]);
			else (D1->A)->seq_al[i][pos]=mod[0];
		      }
		    else fprintf (stderr, "\nWARNING: Could not find Sequence %s", action_list[a]);

		  }
	      }
	    else
	      {
		fp=vfopen (action_list[1], "r");
		while (fscanf (fp, "%s %d %s\n", name, &pos, mod)==3)
		  {

		     i=name_is_in_list(name, (D1->A)->name, (D1->A)->nseq, 100);
		     if (i!=-1)
		       {
			 pos=p[i][pos]-1;
			 if ( strm (mod, "upper"))(D1->A)->seq_al[i][pos]=toupper(A->seq_al[i][pos]);
			 else if ( strm (mod, "lower"))A->seq_al[i][pos]=tolower(A->seq_al[i][pos]);
			 else A->seq_al[i][pos]=mod[0];
		       }
		      else fprintf(stderr, "\nWARNING: Could not find Sequence %s", action_list[1]);
		  }
		vfclose (fp);
	      }
	    D1->S=aln2seq (D1->A);
	  }
       else if ( strm (action, "clean_flag"))
	 {
	   clean_flag=1-clean_flag;
	 }
       else if ( strm  (action, "aln2case"))
	 {
	   D1->A=aln2case_aln (D1->A, ACTION(1), ACTION(2));
	   D1->S=aln2seq(D1->A);
	 }

       else if ( strm5 (action, "convert","upper","lower", "keep", "switchcase"))
	 {
	   b=1;

	   if ( n_actions>1 && is_number (action_list[b]))
	     {
	       lower_value=upper_value=atoi(action_list[b++]);
	     }
	   else if ( n_actions>1 && strm (action_list[b], "gap"))
	     {
	       DST=(Sequence_data_struc*)vcalloc (1,sizeof(Sequence_data_struc));
	       DST->A=aln2gap_cache (D1->A,0);
	       lower_value=0;
	       upper_value=0;
	       b++;
	     }
	   else if (n_actions>1 && action_list[b] && action_list[b][0]=='[')

	     {
	       lower_value=atoi(strtok (action_list[b]+1, "-[]"));
	       upper_value=atoi(strtok (NULL, "-[]"));

	       b++;
	     }
	   else
	     {
	       lower_value=upper_value=-1;
	     }

	   if ( n_actions >b ||strm (action, "keep") )
	     {
	       if ( !RAD->symbol_list)RAD->symbol_list=declare_char (STRING, STRING);
	       RAD->n_symbol=0;
	       if ( strm (action, "keep") )sprintf ( RAD->symbol_list[RAD->n_symbol++], "#-");
	       else
		 {
		   for (a=b; a< n_actions; a++)
		     {
		       sprintf ( RAD->symbol_list[RAD->n_symbol], "%s", action_list[a]);
		       RAD->n_symbol++;
		     }
		 }
	     }

	   for ( value=0; value<=9; value++)
	     {
	       if ( lower_value==-1)value=-1;

	       if ( (value>=lower_value && value<=upper_value)|| value==-1)
		 {
		   if (strm(action,"convert")) D1->A=filter_aln_convert (D1->A, DST?DST->A:NULL,RAD->use_consensus,value,RAD->n_symbol, RAD->symbol_list);
		   else if (strm(action,"upper"))D1->A=filter_aln_lower_upper (D1->A, DST?DST->A:NULL,RAD->use_consensus,value);
		   else if (strm(action,"lower"))D1->A=filter_aln_upper_lower (D1->A, DST?DST->A:NULL,RAD->use_consensus,value);
		   else if (strm(action,"switchcase"))D1->A=filter_aln_switchcase (D1->A, DST?DST->A:NULL,RAD->use_consensus,value);
		 }
	       else
		 {
		   if (strm(action,"keep")) D1->A=filter_aln_convert (D1->A, DST?DST->A:NULL,RAD->use_consensus,value,RAD->n_symbol, RAD->symbol_list);
		 }
	       if (value==-1)break;

	     }

	   /*free_sequence (D1->S,(D1->S)->nseq);*/
	   if (!D1->S)D1->S=aln2seq (D1->A);
	 }
	else if ( strm ( action, "count_pairs"))
	  {
	    int a, b,c,v, **matrix;
	    Alignment *A;
	    matrix=declare_int (300,300);
	    A=D1->A;
	    for ( a=0; a< A->nseq-1; a++)
	      for (b=0; b< A->nseq; b++)
		for (c=0; c<A->len_aln; c++)
		  matrix[(int)A->seq_al[a][c]][(int)A->seq_al[b][c]]++;
	    for ( a=0; a<255; a++)
	      for ( b=a; b<256; b++)
		{
		  v=matrix[a][b]+matrix[b][a];
		  if (v)fprintf ( stdout, "\n%c %c %d", a, b, v);
		}
	    myexit (EXIT_SUCCESS);
	  }
	else if ( strm (action, "count_misc"))
	  {
	    count_misc (D1->A, (!D2)?NULL:D2->A);
	  }
       else if ( strm (action, "count"))
	 {
	   b=1;
	   if ( n_actions>1 && is_number (action_list[b]))
	     {
	       lower_value=upper_value=atoi(action_list[b++]);
	     }
	   else if (n_actions>1 && action_list[b] && action_list[b] && action_list[b][0]=='[')

	     {
	       lower_value=atoi(strtok (action_list[b]+1, "-[]"));
	       upper_value=atoi(strtok (NULL, "-[]"));

	       b++;
	     }
	   else
	     {
	       lower_value=upper_value=-1;
	     }
	   if ( n_actions >b)
	     {
	       if ( !RAD->symbol_list)RAD->symbol_list=declare_char (STRING, STRING);
	       RAD->n_symbol=0;
	       for (a=b; a< n_actions; a++)
		 {
		   sprintf ( RAD->symbol_list[RAD->n_symbol], "%s", action_list[a]);
		   RAD->n_symbol++;
		 }
	     }
	   for ( value=lower_value; value<=upper_value; value++)
	     {
	       count_table=count_in_aln (D1->A, DST?DST->A:NULL,value,RAD->n_symbol, RAD->symbol_list, count_table);
	     }
	   for ( a=0; a<RAD->n_symbol; a++)
	     {
	       fprintf ( stdout, "%s %d\n", RAD->symbol_list[a], count_table[a]);
	     }
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	   vfree(count_table);
	   myexit(EXIT_SUCCESS);
	 }
       else if ( strm (action, "species_weight"))
	 {
	   seq_weight2species_weight (D1->A, D2->S);
	   exit (0);
	 }
       else if ( strm (action, "aln2voronoi"))
	 {
	   aln2voronoi_weights (D1->A);

	 }
       else if ( strm (action, "msa_weight"))
	 {
	   int random_value;
	   char command [LONG_STRING];
	   char aln_name[FILENAMELEN];
	   char tree_name[FILENAMELEN];
	   char dist_matrix_name[FILENAMELEN];
	   char weight_name[FILENAMELEN];
	   char method_4_msa_weights[1000];

	   if ( n_actions==1)
	     {
	       fprintf ( stderr, "\nError: msa_weight requires a weight_method");
	     }

	   sprintf ( method_4_msa_weights, "%s", (get_env_variable ("METHOD_4_MSA_WEIGHTS",NO_REPORT))?get_env_variable ("METHOD_4_MSA_WEIGHTS",NO_REPORT):METHOD_4_MSA_WEIGHTS);

	   /*1 Computation of the tree and the distance matrix*/
	   random_value=addrand ((unsigned long) 100000)+1;
	   sprintf (aln_name, "%d.aln", random_value);
	   sprintf (tree_name, "%d.ph", random_value);
	   sprintf (dist_matrix_name, "%d.dst", random_value);
	   sprintf (weight_name, "%d.weight", random_value);
	   output_fasta_aln (aln_name, D1->A);

	   sprintf ( command, "clustalw -infile=%s -tree -outputtree=dist %s", aln_name, TO_NULL_DEVICE);
	   my_system ( command);
	   sprintf ( command, "%s -method %s -aln %s -tree %s -dmatrix %s -weightfile %s %s",method_4_msa_weights, action_list[1],aln_name, tree_name, dist_matrix_name,weight_name, TO_NULL_DEVICE);
	   my_system ( command);

	   (D1->A)->S=aln2seq (D1->A);
	   ((D1->A)->S)->W=read_seq_weight ( (D1->A)->name, (D1->A)->nseq,weight_name);
	   vremove (weight_name);
	   vremove (aln_name);
	   vremove (tree_name);
	   vremove (dist_matrix_name);
	 }
       else if ( strm (action, "pavie_seq2random_seq"))
	 {
	   D1->S=pavie_seq2random_seq (D1->S, action_list[1]);
	   D1->A=seq2aln (D1->S,NULL,1);
	 }
       else if ( strm ( action, "pavie_seq2noisy_seq"))
	 {
	   /*<amount of noise: 0-100> (<alp>)*/

	   D1->S=pavie_seq2noisy_seq (D1->S, atoi(action_list[1]),ACTION(2));
	   D1->A=seq2aln (D1->S,NULL,1);
	 }
       else if ( strm (action, "pavie_seq2pavie_mat"))
	 {

	   pavie_seq2trained_pavie_mat ( D1->S, (n_actions==2)?action_list[1]:NULL);
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "pavie_seq2pavie_aln"))
	 {

	   pavie_seq2pavie_aln ( D1->S, action_list[1], ACTION(2));
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "pavie_seq2pavie_dm"))
	 {
	    if (strstr (ACTION2(2,""), "_MSA_"))
	      D1->S=aln2seq_main(D1->A, KEEP_GAP);


	   pavie_seq2pavie_aln ( D1->S, action_list[1],  const_cast<char*>( (n_actions==3)?action_list[2]:"_MATDIST_") );
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action, "pavie_seq2pavie_msa"))
	 {
	   D1->A=pavie_seq2pavie_msa ( D1->S, action_list[1], (n_actions==3)?action_list[2]:NULL);
	 }
       else if ( strm (action, "pavie_seq2pavie_tree"))
	 {
	   D1->T=pavie_seq2pavie_tree ( D1->S, action_list[1], (n_actions==3)?action_list[2]:NULL);
	 }
       else if ( strm (action, "pavie_seq2pavie_sort"))
	 {
	   D1->A=pavie_seq2pavie_sort ( D1->S, action_list[1], (n_actions==3)?action_list[2]:NULL);
	 }

       else if ( strm (action, "aln2mat_diaa"))
	 {
	   aln2mat_diaa (D1->S);
	 }
       else if ( strm (action, "aln2proba_mat"))
	 {
	   aln2proba_mat(D1->S);
	 }

       else if ( strm (action, "aln2mat"))
	 {
	   aln2mat(D1->S);
	 }

       else if ( strm (action, "seq2latmat"))
	 {
	   seq2latmat ( D1->S, "stdout");
	   myexit (EXIT_SUCCESS);
	 }
       else if ( strm (action , "rm_target_pdb"))
	 {
	   int i, j;
	   char *buf;

	   for (i=0; i< (D1->A)->nseq; i++)
	     {
	       j=1;buf=(D1->A)->name[i];
	       while (buf[j]!='_' && buf[j-1]!='_' && buf[j]!='\0')j++;
	       buf[j]='\0';
	     }
	 }
       else if ( strm ( action, "mat2cmp"))
	 {
	   double *r;
	   r=mat2cmp (D1->M, D2->M);
	   fprintf ( stdout, "\nMATRIX COMPARISON: R=%.3f R2=%.3f On %d pairs of values\n", (float)r[0], (float)r[1], (int)r[2]);
	   myexit (EXIT_SUCCESS);
	 }
//Special modes
       else if ( strm ( action, "overaln_list"))
	 {
	   float *re, tre=0,sn, tsn=0, sp, tsp=0;
	   int p1,p2,p3, t, f;
	   FILE *fp;
	   char fname [100];
	   Alignment **LA;
	   Alignment **LB;

	   HERE ("F P1 P2 P3 T");

	   t=ATOI_ACTION(1);
	   f=ATOI_ACTION(2);
	   p1=ATOI_ACTION(3);
	   p2=ATOI_ACTION(4);
	   p3=ATOI_ACTION(5);



	   LA=(Alignment**)vcalloc ((D1->A)->nseq, sizeof (Alignment*));
	   LB=(Alignment**)vcalloc ((D2->A)->nseq, sizeof (Alignment*));
	   for (a=0; a<(D1->A)->nseq; a++)
	     {
	        LA[a]=main_read_aln ((D1->A)->name[a], NULL);
		LB[a]=main_read_aln ((D2->A)->name[a], NULL);
	     }

	   for ( a=0; a<(D1->A)->nseq; a++)
	     {
	       Alignment *A, *B;
	       A=LA[a];
	       B=LB[a];
	       re=analyze_overaln (A, B, "_case_l_",t,f,p1,p2,p3);
	       fprintf (stdout, "\n%d: sn: %.2f sp: %.2f re: %.2f F: %d P: %d P2: %d T: %d",a, re[0],re[1],re[2],f, p1,p2,t);
	       tsn+=re[0];
	       tsp+=re[1];
	       tre+=re[2];
	       vfree(re);
	     }
	   fprintf (stdout, "\nTOT: sn: %.2f sp: %.2f re: %.2f F: %d P: %d P2: %d T: %d", tsn/(D1->A)->nseq,tsp/(D1->A)->nseq, tre/(D1->A)->nseq,f,p1,p2,t);

	   myexit (0);
	 }
       else if ( strm ( action, "overaln_list_scan"))
	 {
	   float *re, tre=0, tsn=0, tsp;
	   int p1,p2, p3, t, f;
	   FILE *fp;
	   char fname [100];
	   Alignment **LA;
	   Alignment **LB;

	   if ( ACTION(1))sprintf ( fname, "%s", ACTION(1));
	   else sprintf ( fname, "scan_results.txt");

	   fprintf ( stdout, "SCAN Results will be ouput in %s\n", fname);


	   LA=(Alignment**)vcalloc ((D1->A)->nseq, sizeof (Alignment*));
	   LB=(Alignment**)vcalloc ((D2->A)->nseq, sizeof (Alignment*));
	   for (a=0; a<(D1->A)->nseq; a++)
	     {
	        LA[a]=main_read_aln ((D1->A)->name[a], NULL);
		LB[a]=main_read_aln ((D2->A)->name[a], NULL);
	     }
	   for (f=32; f<=40; f++)
	     {
	       for (p1=90; p1<=100; p1+=5)
		 {
		   for ( t=1; t<=3; t++)
		     {
		       for (p2=0; p2<=40; p2+=5)
			 {
			   for (p3=0;p3<=0;p3+=5)
			     {
			       tre=tsn=tsp=0;
			       for ( a=0; a<(D1->A)->nseq; a++)
				 {
				   Alignment *A, *B;
				   A=LA[a];
				   B=LB[a];
				   re=analyze_overaln (A, B, "_case_l_",t,f,p1,p2,p3);

				   tsn+=re[0];
				   tsp+=re[1];
				   tre+=re[2];
				   vfree (re);
				 }
			       fp=vfopen (fname, "a");
			       fprintf (fp, "\nTOT: sn: %.2f sp: %.2f re: %.2f P: %d P2: %d P3: %d T: %d F: %d", tsn/(D1->A)->nseq,tsp/(D1->A)->nseq, tre/(D1->A)->nseq, p1,p2, p3,t,f);
			       fprintf (stderr, "\nTOT: sn: %.2f sp: %.2f re: %.2f P: %d P2: %d P3: %d T: %d F: %d", tsn/(D1->A)->nseq,tsp/(D1->A)->nseq, tre/(D1->A)->nseq, p1,p2, p3,t,f);
			       vfclose (fp);
			     }
			 }
		     }
		 }
	     }
	   myexit (0);
	 }
       else if ( strm ( action, "overaln"))//Evaluate the capacity to predict over-aligned regions
	 {
	   OveralnP *F;
	   F=(OveralnP*)vcalloc (1, sizeof (OveralnP));
	   //al1: ref
	   //al2: alignment
	   //ATOI(1): P (0-100)
	   //ATOI(2): T (0-9)

	   float *r;
	   DST=(Sequence_data_struc*)vcalloc (1,sizeof(Sequence_data_struc));
	   DST->A=aln2gap_cache (D1->A,0);
	   lower_value=0;
	   upper_value=0;
	   D1->A=filter_aln_upper_lower (D1->A, DST->A, 0, 0);

	   sprintf (F->mode, "%s", ((s=get_string_variable ("overaln_mode")))?s:"lower");
	   if (!strm (F->mode, "lower") && !strstr (F->mode, "unalign"))printf_exit (EXIT_FAILURE,stderr,"\nERROR: unknown overal_mode in overal output [%s] [FATAL:%s]", F->mode, PROGRAM);

	   if (int_variable_isset ("overaln_threshold"))F->t=get_int_variable ("overaln_threshold");
	   if (int_variable_isset ("overaln_target"))F->f=get_int_variable ("overaln_target");
	   if (int_variable_isset ("overaln_P1"))F->f=get_int_variable ("overaln_P1");
	   if (int_variable_isset ("overaln_P1"))F->f=get_int_variable ("overaln_P2");
	   if (int_variable_isset ("overaln_P1"))F->f=get_int_variable ("overaln_P3");
	   if (int_variable_isset ("overaln_P1"))F->f=get_int_variable ("overaln_P4");//F P1 P2 P3 T;

	   D2->A=aln2clean_pw_aln (D2->A, F);
	   r=aln2pred (D1->A, D2->A,"case_l_");
	   fprintf ( stdout, "sn %.2f sp %.2f re %.2f\n", r[0], r[1], r[2]);
	   myexit (0);
	 }
       else if ( strm (action, "seq2ngs"))
	 {
	 
	   int cov=10;
	   int rl=50;
	   
	   int len,ni,nj,nk,nl;
	   nl=0;
	   for (ni=0; ni<(D1->S)->nseq; ni++)
	       {
		 len=strlen ((D1->S)->seq[ni]);
		 for (nj=0; nj<len-50; nj+=rl/cov)
		   {
		     fprintf (stdout, ">%d\n", ++nl);
		     for (nk=nj; nk<nj+rl && nk<len; nk++) 
		       {
			 fprintf (stdout, "%c", (D1->S)->seq[ni][nk]);
		       }
		     fprintf (stdout, "\n");
		   }
	       }
	   exit (0);
	 }
       		     
//JM_START
       else if ( strm ( action, "aln2hitMat"))
	 {
 		aln2hitMat(D1->A, ACTION(1));
 		myexit (EXIT_SUCCESS);
	 }
//JM_END

       else
	 {
	   fprintf ( stderr, "\nWARNING: ACTION %s UNKNOWN and IGNORED\n", action);
	 }

     }


void aln2mat_diaa (Sequence *S)
{
  int a, aa1, aa2, aa3, aa4;
  int s1, s2, p;
  Alignment *A;

  int ****m;
  int **c;
  int naa=0;
  int count=0;
  double Delta=0.00001;
  int *alp;
  int tot,u;
  double observed, expected, f_diaa1, f_diaa2, v;


  alp=(int*)vcalloc (256, sizeof (int));
  for (a=0; a<26; a++)alp[a+'a']=1;
  alp['b']=0;
  alp['j']=0;
  alp['o']=0;
  alp['u']=0;
  alp['x']=0;
  alp['z']=0;

  m=(int****)declare_arrayN (4,sizeof (int),26,26,26,26);
  c=(int**)declare_arrayN  (2,sizeof (int),26,26);

  for ( a=0; a< S->nseq; a++)
    {
      fprintf ( stderr, "%s\n", S->name[a]);
      A=main_read_aln (S->name[a],NULL);
      for (s1=0; s1<A->nseq; s1++)lower_string (A->seq_al[s1]);

      for ( s1=0; s1<A->nseq-1; s1++)
	for (s2=s1+1; s2<A->nseq; s2++)
	  {
	    for (p=0; p<A->len_aln-1; p++)
	      {

		u =alp[aa1=A->seq_al[s1][p]];
		u+=alp[aa2=A->seq_al[s1][p+1]];
		u+=alp[aa3=A->seq_al[s2][p]];
		u+=alp[aa4=A->seq_al[s2][p+1]];

		if ( u==4)
		  {
		    aa1-='a';aa2-='a';aa3-='a'; aa4-='a';

		    c[aa1][aa2]++;
		    c[aa3][aa4]++;
		    m[aa1][aa2][aa3][aa4]++;
		    count+=2;
		  }
	      }
	  }
      free_aln (A);
    }
  fprintf ( stdout, "# DIAA_MATRIX_FORMAT_01\n");
  naa=26;
  for (aa1=0; aa1<naa; aa1++)
    for (aa2=0; aa2<naa; aa2++)
      for (aa3=0; aa3<naa; aa3++)
	for (aa4=0; aa4<naa;aa4++)
	  {
	    u =alp[aa1+'a'];
	    u+=alp[aa2+'a'];
	    u+=alp[aa3+'a'];
	    u+=alp[aa4+'a'];

	    if ( u==4)
	      {
		tot=m[aa1][aa2][aa3][aa4]+m[aa3][aa4][aa1][aa2];
		observed=((double)tot)/(double)((double)count/(double)2);
		f_diaa1=(double)c[aa1][aa2]/(double)count;
		f_diaa2=(double)c[aa3][aa4]/(double)count;

		expected=f_diaa1*f_diaa2;
		if (expected<Delta)v=0;
		else if (observed<Delta)v=-100;
		else
		  {
		    v=log(observed/expected)*10;
		  }
	    // if (tot>0)fprintf ( stdout, "TEST C=%d expected=%.4f observed=%.4f v=%.4f [%d %d %d][%d] tot=%d\n", count, (float)expected, (float)observed, (float) v, c[aa1][aa2], c[aa3][aa4], count, m[aa1][aa2][aa3][aa4], tot);
		fprintf ( stdout, "%c%c %c%c %d %d\n", aa1+'a', aa2+'a', aa3+'a', aa4+'a', (int)v, m[aa1][aa2][aa3][aa4]+ m[aa3][aa4][aa1][aa2]);
	      }
	  }
  myexit (EXIT_SUCCESS);
}
void aln2proba_mat (Sequence *S)
{
  int a, aa1, aa3;
  int s1, s2, p;
  Alignment *A;
  int **mat;
  int **m;
  int *c;
  int naa=0;
  int count=0;
  double Delta=0.00001;
  int *alp;
  int tot,u;
  double observed, expected, f_diaa1, f_diaa2, v;
  char *balp;

  balp=(char*)vcalloc ( 256, sizeof (char));
  for (a=0; a<strlen (BLAST_AA_ALPHABET); a++)balp[BLAST_AA_ALPHABET[a]]=a;

  mat=declare_int (256, 256);
  alp=(int*)vcalloc (256, sizeof (int));
  for (a=0; a<26; a++)alp[a+'a']=1;
  alp['b']=0;
  alp['j']=0;
  alp['o']=0;
  alp['u']=0;
  alp['x']=0;
  alp['z']=0;

  m=(int**)declare_arrayN (2,sizeof (int),26,26);
  c=(int*)declare_arrayN  (1,sizeof (int),26);
  
  for ( a=0; a< S->nseq; a++)
    {
      fprintf ( stderr, "%s\n", S->name[a]);
      A=main_read_aln (S->name[a],NULL);
      for (s1=0; s1<A->nseq; s1++) lower_string (A->seq_al[s1]);

      for ( s1=0; s1<A->nseq-1; s1++)
	for (s2=s1+1; s2<A->nseq; s2++)
	  {
	    for (p=0; p<A->len_aln-1; p++)
	      {

		u =alp[aa1=A->seq_al[s1][p]];
		u+=alp[aa3=A->seq_al[s2][p]];

		if ( u==2)
		  {
		    aa1-='a';aa3-='a';

		    c[aa1]++;
		    c[aa3]++;
		    m[aa1][aa3]++;
		    count+=2;
		  }
	      }
	  }
      free_aln (A);
    }
  fprintf ( stdout, "# PROBABILITY_MATRIX_01\n");
  naa=26;
  
  for (aa1=0; aa1<naa; aa1++)
      for (aa3=0; aa3<naa; aa3++)
	  {
	    u =alp[aa1+'a'];
	    u+=alp[aa3+'a'];
	    fprintf (stdout,"%c %c %.5f\n",aa1+'a', aa3+'a', (float)(m[aa1][aa3]+m[aa3][aa1])/(float)count); 
	  }
 
  myexit (EXIT_SUCCESS);
}

void aln2mat (Sequence *S)
{
  int a, aa1, aa3;
  int s1, s2, p;
  Alignment *A;
  int **mat;
  int **m;
  int *c;
  int naa=0;
  int count=0;
  double Delta=0.00001;
  int *alp;
  int tot,u;
  double observed, expected, f_diaa1, f_diaa2, v;
  char *balp;

  balp=(char*)vcalloc ( 256, sizeof (char));
  for (a=0; a<strlen (BLAST_AA_ALPHABET); a++)balp[BLAST_AA_ALPHABET[a]]=a;

  mat=declare_int (256, 256);
  alp=(int*)vcalloc (256, sizeof (int));
  for (a=0; a<26; a++)alp[a+'a']=1;
  alp['b']=0;
  alp['j']=0;
  alp['o']=0;
  alp['u']=0;
  alp['x']=0;
  alp['z']=0;

  m=(int**)declare_arrayN (2,sizeof (int),26,26);
  c=(int*)declare_arrayN  (1,sizeof (int),26);

  for ( a=0; a< S->nseq; a++)
    {
      fprintf ( stderr, "%s\n", S->name[a]);
      A=main_read_aln (S->name[a],NULL);
      for (s1=0; s1<A->nseq; s1++)lower_string (A->seq_al[s1]);

      for ( s1=0; s1<A->nseq-1; s1++)
	for (s2=s1+1; s2<A->nseq; s2++)
	  {
	    for (p=0; p<A->len_aln-1; p++)
	      {

		u =alp[aa1=A->seq_al[s1][p]];
		u+=alp[aa3=A->seq_al[s2][p]];

		if ( u==2)
		  {
		    aa1-='a';aa3-='a';

		    c[aa1]++;
		    c[aa3]++;
		    m[aa1][aa3]++;
		    count+=2;
		  }
	      }
	  }
      free_aln (A);
    }
  fprintf ( stdout, "# MONOAA_MATRIX_FORMAT_01\n");
  naa=26;
  for (aa1=0; aa1<naa; aa1++)
      for (aa3=0; aa3<naa; aa3++)
	  {
	    u =alp[aa1+'a'];
	    u+=alp[aa3+'a'];

	    if ( u==2)
	      {
		tot=m[aa1][aa3]+m[aa3][aa1];
		observed=((double)tot)/(double)((double)count/(double)2);
		f_diaa1=(double)c[aa1]/(double)count;
		f_diaa2=(double)c[aa3]/(double)count;

		expected=f_diaa1*f_diaa2;
		if (expected<Delta)v=0;
		else if (observed<Delta)v=-100;
		else
		  {
		    v=log(observed/expected)/(log(2)/2);
		  }
	    // if (tot>0)fprintf ( stdout, "TEST C=%d expected=%.4f observed=%.4f v=%.4f [%d %d %d][%d] tot=%d\n", count, (float)expected, (float)observed, (float) v, c[aa1][aa2], c[aa3][aa4], count, m[aa1][aa2][aa3][aa4], tot);
		//fprintf ( stdout, "%c %c %d %d\n", aa1+'A', aa3+'A', (int)v, m[aa1][aa3]+ m[aa3][aa1]);
		mat[aa1][aa3]=(int)v;
	      }
	  }
  output_blast_mat (mat, "stdout");
  myexit (EXIT_SUCCESS);
}


int **seq2latmat ( Sequence *S, char *fname)
{
  int a, b, r0, r1;
  char *aa;
  int naa;
  int *count, tot;
  int **mat;
  double observed, expected;
  FILE *fp;

  fp=vfopen (fname, "w");

  count=(int*)vcalloc ( 256, sizeof (int));
  mat=declare_int (256, 256);

  naa=strlen ( BLAST_AA_ALPHABET);
  aa=(char*)vcalloc ( naa+2, sizeof (char));
  sprintf ( aa, "%s", BLAST_AA_ALPHABET);
  lower_string (aa);

  for ( tot=0,a=0; a< S->nseq; a++)
    {
      ungap (S->seq[a]);
      for ( b=1; b<S->len[a]; b++)
	{
	  r0=tolower(S->seq[a][b-1]);
	  r1=tolower(S->seq[a][b]);

	  mat[r0][r1]++;
	  //count[r0]++;
	  count[r1]++;
	  tot++;
	}
    }
  for ( a=0; a< naa; a++)
    for (b=0; b< naa; b++)
      {
	if ( aa[a]=='*' || aa[b]=='*');
	else
	  {
	    expected=((double)count[(int)aa[a]]/(double)tot)* ((double)count[(int)aa[b]]/(double)tot)*(double)tot;
	    observed=((double)mat[(int)aa[a]][(int)aa[b]]);

	    /*
	      fprintf ( stderr, "\n%c=%d %c=%d Tot=%d Obs=%d Exp=%d\n", aa[a],count[aa[a]], aa[b],count[aa[b]],tot, mat[aa[a]][aa[b]],(int)expected);
	      fprintf ( stderr, "\n%d", mat[aa[a]][aa[b]]);
	      fprintf ( stderr, "\n%d", mat[aa[a]][aa[b]]);
	    */
	    mat[(int)aa[a]][(int)aa[b]]=(expected==0 || observed==0)?0:((int)10*log((observed/expected)));
	  }
      }

  fprintf (fp,"# BLAST_MATRIX FORMAT\n#ALPHABET=%s\n#TRANSITION MATRIX TRAINED ON %d Sequence\n#", BLAST_AA_ALPHABET, S->nseq);
  for (a=0; a< naa; a++)fprintf ( fp, "%3c ", toupper(aa[a]));
  fprintf (fp,"\n");
  for (a=0; a< naa; a++)
    {

      fprintf (fp, "%c", toupper(aa[a]));
      for ( b=0; b< naa; b++)
	{
	  fprintf (fp, "%3d ", mat[(int)aa[a]][(int)aa[b]]);
	}
      fprintf ( fp, "\n");
    }
  vfclose (fp);
  vfree (count);
  vfree (aa);

  return mat;
}

double* mat2cmp ( int **mat1, int **mat2)
{
  int a, b, n, x, y;
  double **list, *r;
  if ( !mat1 || !mat2)
    {
      fprintf ( stderr, "\nERROR: mat2cmp needs two matrices [FATAL:%s]", PROGRAM);
      myexit (EXIT_FAILURE);
    }

  for (n=0, a=0; a< 256; a++)
    for ( b=0; b<256; b++)
      {
	x=mat1[a][b];
	y=mat2[a][b];
	if ( x|| y)n++;
      }
  if ( n==0) return 0;
  list=declare_double (n, 2);

  for (n=0, a=0; a<256; a++)
    for ( b=0; b<256; b++)
      {
	x=mat1[a][b];
	y=mat2[a][b];
	if ( x || y)
	  {
	    list[n][0]=x;
	    list[n][1]=y;
	    n++;
	  }
      }
  r=return_r (list, n);
  free_double(list, -1);
  return r;
}

int ** read_blast_matrix ( char *mat_name)
        {
	FILE *fp;
	int n_aa,aa1, aa2;
	int a, b, c;
	int **matrix;
	int value;
	char alp[257];
	char *buf=NULL;

	buf=(char*)vcalloc (300, sizeof (char));
	fp=vfopen (mat_name, "r");
	while ((c=fgetc(fp))=='#')buf=vfgets(buf, fp);
	buf=vfgets(buf,fp);
	
	a=n_aa=0;
	while ((c=buf[a++])!='\0')if (!isspace(c))alp[n_aa++]=c;
	alp[n_aa]='\0';
	
	matrix=declare_int (256,256);
	vfree ( matrix[30]);
	matrix[30]=(int*)vcalloc(10000, sizeof (int));
	
	for ( a=0; a< n_aa; a++)
	    {
	    fscanf ( fp, "%s ", buf);

	    aa1=buf[0];

	    if ( aa1!=alp[a])
		{
		fprintf ( stderr, "\nParsing_error when reading blast_matrix %s:\n%c %c",mat_name, aa1,alp[a]);
		fprintf ( stderr, "\n%c ", fgetc(fp));
		myexit (EXIT_FAILURE);
		}
	    for ( b=0; b<n_aa; b++)
	        {
		aa2=alp[b];
		fscanf ( fp, "%d ", &value);
		if (is_gap(aa1) || is_gap(aa2))
		  {
		    int c1, c2;
		    continue;
		    c1=(is_gap(aa1))?GAP_CODE:aa1;
		    c2=(is_gap(aa2))?GAP_CODE:aa2;
		    if ( c1==GAP_CODE && c2==GAP_CODE)
		      matrix[c1][c2]=value;
		    else if ( c1==GAP_CODE)
		      {
			matrix[c1][tolower(c2)]=value;
			matrix[c1][toupper(c2)]=value;
		      }
		    else
		      {
			matrix[tolower(c1)][c2]=value;
			matrix[toupper(c1)][c2]=value;
		      }
		  }
		else if ( aa1!='*' && aa2!='*')
		  {
		    matrix[tolower(aa1)-'a'][tolower(aa2)-'a']=value;
		  }
		}
	    fscanf(fp, "\n");
	    }
	fclose (fp);

	return matrix;
	}
int ** read_blast_matrix_old ( char *mat_name)
        {
	FILE *fp;
	int n_aa,aa1, aa2;
	int a, b, c;
	int **matrix;
	int value;
	char sbuf[VERY_LONG_STRING];
	char buf[2];
	char alp[257];

	for(a=0; a<256; a++)alp[0]='\0';
	matrix=declare_int (256,256);
	vfree ( matrix[30]);
	matrix[30]=(int*)vcalloc(10000, sizeof (int));
	fp=vfopen ( mat_name, "r");
	while ( (c=fgetc(fp))=='#' || isspace(c) )
	  {
	    char *p;
	    fgets ( sbuf, VERY_LONG_STRING, fp);
	    if ( (p=strstr (sbuf, "ALPHABET")))
	      sscanf (p, "ALPHABET=%s", alp);
	  }
	
	ungetc(c, fp);
	lower_string (alp);
	n_aa=strlen (alp);

	





	for ( a=0; a< n_aa; a++)
	    {
	    fscanf ( fp, "%s ", buf);

	    aa1=tolower(buf[0]);

	    if ( aa1!=alp[a])
		{
		fprintf ( stderr, "\nParsing_error when reading blast_matrix %s:\n%c %c",mat_name, aa1,alp[a]);
		fprintf ( stderr, "\n%c ", fgetc(fp));
		myexit (EXIT_FAILURE);
		}
	    for ( b=0; b<n_aa; b++)
	        {
		aa2=tolower ((char) alp[b]);
		fscanf ( fp, "%d ", &value);
		if (is_gap(aa1) || is_gap(aa2))
		  {
		    int c1, c2;
		    c1=(is_gap(aa1))?GAP_CODE:aa1;
		    c2=(is_gap(aa2))?GAP_CODE:aa2;
		    if ( c1==GAP_CODE && c2==GAP_CODE)
		      matrix[c1][c2]=value;
		    else if ( c1==GAP_CODE)
		      {
			matrix[c1][tolower(c2)]=value;
			matrix[c1][toupper(c2)]=value;
		      }
		    else
		      {
			matrix[tolower(c1)][c2]=value;
			matrix[toupper(c1)][c2]=value;
		      }
		  }
		else if ( aa1!='*' && aa2!='*')
		  {
		    matrix[tolower(aa1)-'A'][tolower(aa2)-'A']=value;
		    matrix[toupper(aa1)-'A'][toupper(aa2)-'A']=value;
		    matrix[tolower(aa1)-'A'][toupper(aa2)-'A']=value;
		    matrix[toupper(aa1)-'A'][tolower(aa2)-'A']=value;
		  }
		}
	    fscanf(fp, "\n");
	    }
	fclose (fp);

	return matrix;
	}

int output_blast_mat (int **mat, char *fname)
{
  return output_mat(mat, fname, BLAST_AA_ALPHABET, 'a');

}
int output_header_mat (int **mat, char *fname)
{
  int a, b,l;
  FILE *fp;
  int naa;

  char raa[]="ABCDEFGHIKLMNPQRSTVWXYZ";
  char *aa;

  naa=strlen (raa);
  aa=(char*)vcalloc ( naa+2, sizeof (char));
  sprintf ( aa, "%s",raa);
  lower_string (aa);

  fp=vfopen (fname, "w");
  fprintf ( fp, "int new_mat[]={\n");
  l=strlen (aa);
  for (a=0; a<naa; a++)
    {
      for (b=0; b<=a; b++)
	{
	  fprintf (fp, "%3d, ", mat[aa[a]-'a'][aa[b]-'a']);
	}
      fprintf (fp, "\n");
    }
  fprintf ( fp, "}");
  vfclose (fp);
}
int output_mat (int **mat, char *fname, char *alp, int offset)
{
  char *aa;
  int a,b, naa;
  FILE *fp;



  naa=strlen (alp);
  aa=(char*)vcalloc ( naa+2, sizeof (char));
  sprintf ( aa, "%s",alp);
  lower_string (aa);
  if (!(fp=vfopen (fname, "w")))return 0;
  fprintf (fp,"# BLAST_MATRIX FORMAT\n#ALPHABET=%s\n  ",alp);
  for (a=0; a< naa; a++)fprintf ( fp, "%5c ", toupper(aa[a]));
  fprintf (fp,"\n");
  for (a=0; a< naa; a++)
    {

      fprintf (fp, "%c", toupper(aa[a]));
      for ( b=0; b< naa; b++)
	{
	  if (aa[a]!='*' && aa[b]!='*')
	    fprintf (fp, " %5d", mat[aa[a]-offset][aa[b]-offset]);
	  else
	    fprintf (fp, " %5d", 0);
	}
      fprintf ( fp, "\n");
    }
  vfree (aa);
  vfclose (fp);
  return 1;
}

void output_pavie_mat (int **mat, char *fname, double gep, char *alp)
{
  int n, a, b;
  FILE *fp;

  n=strlen (alp);
  fp=vfopen (fname, "w");
  fprintf (fp,"# PAVIE_MATRIX FORMAT\n#ALPHABET=%s\n",alp);

  for(a=0; a< n; a++)
     {
       for ( b=a; b<n; b++)
	{
	  fprintf (fp, "%c %c %.3f\n", toupper(alp[a]), toupper(alp[b]), (float)mat[alp[a]-'A'][alp[b]-'A']/PAVIE_MAT_FACTOR);
	}
     }
   if ( gep!=UNDEFINED)fprintf ( fp, "- - %.3f\n", gep/PAVIE_MAT_FACTOR);
   vfclose(fp);
 }

int ** read_pavie_matrix ( char *mat_name)
        {
	FILE *fp;
	int c, n_aa;
	char aa1, aa2;
	float v;
	int **matrix;
	char sbuf[VERY_LONG_STRING];
	char alp[257];
	int gep=UNDEFINED;

	matrix=declare_int (256,256);


	fp=vfopen ( mat_name, "r");
	while ( (c=fgetc(fp))=='#' || isspace(c) )
	  {
	    fgets ( sbuf, VERY_LONG_STRING, fp);
	    if ( sscanf (sbuf, "ALPHABET=%s", alp)==1);
	  }
	ungetc(c, fp);

	n_aa=strlen (alp);
	while ( fgets ( sbuf, VERY_LONG_STRING, fp)!=NULL)
	    {
	      aa1=aa2='Z';
	      if (sscanf (sbuf, "%c %c %f",&aa1, &aa2, &v)==3)
		{
		  v*=PAVIE_MAT_FACTOR;
		  if (aa1=='-' && aa2=='-')gep=v;
		  else
		    {
		      matrix[tolower(aa1)-'A'][tolower(aa2)-'A']=v;
		      matrix[toupper(aa1)-'A'][toupper(aa2)-'A']=v;
		      matrix[tolower(aa1)-'A'][toupper(aa2)-'A']=v;
		      matrix[toupper(aa1)-'A'][tolower(aa2)-'A']=v;

		      matrix[tolower(aa2)-'A'][tolower(aa1)-'A']=v;
		      matrix[toupper(aa2)-'A'][toupper(aa1)-'A']=v;
		      matrix[tolower(aa2)-'A'][toupper(aa1)-'A']=v;
		      matrix[toupper(aa2)-'A'][tolower(aa1)-'A']=v;
		    }
	      }
	    }
	if ( gep!=UNDEFINED)
	  {
	    int a;
	    for (a=0; a< n_aa; a++)
	      {
		if (!matrix[tolower(alp[a])-'A'][GAP_CODE])
		  {
		    matrix[tolower(alp[a])-'A'][GAP_CODE]=gep;
		    matrix[toupper(alp[a])-'A'][GAP_CODE]=gep;
		  }
	      }
	  }
	vfclose (fp);
	return matrix;
	}

Sequence *seq2year ( Sequence *S, int modulo)
{
  int a, b, y;
  int first;
  char *s;
  char new_channel[100];

  sprintf( new_channel, "_agechannel%d",modulo);

  for ( a=0; a<S->nseq; a++)
    {
      if (S->seq_comment[a] && (s=strstr(S->seq_comment[a], "_FIRSTYEAR")))
	{
	  sscanf (s, "_FIRSTYEAR%d_", &first);
	}
      else first=1;

      for ( y=first,b=0; b<S->len[a]; b++)
	{
	  if ( !is_gap(S->seq[a][b]))
	    {
	      S->seq[a][b]='a'+((y/modulo))%10;
	      y++;
	    }
	}
      if ( (s=strstr ( S->name[a], "_agechannel")))
	   {
	     sprintf ( s, "%s", new_channel);
	   }
      else strcat (S->name[a], new_channel);
    }
  return S;
}

Sequence* output_n_pavie_age_channel (Sequence *S, char *name, int n)
{
  int x, a;
  if (!n)n=2;


  for ( x=1,a=0; a< n; a++, x*=10)
    {
      S=output_pavie_age_channel(S, name,x);
    }
return S;
}




Sequence* output_pavie_age_channel (Sequence *S, char *name, int modulo)
  {
    Alignment *A;
    FILE *fp;
    static int display;
    char mat_list_name[100];
    char seq_list[1000];
    char mat_name[1000];
    char *tmp;

    sprintf ( mat_list_name, "%s_pavie_age_matrix.mat_list", name);
    sprintf (seq_list, "%s_age_channel.fasta",name);

    if ( display==0 )
      {
	if (check_file_exists(seq_list))vremove (seq_list);
	if (check_file_exists(mat_list_name))vremove (mat_list_name);
      }
    sprintf (mat_name, "%s_age_mat_mod%d.mat",name, modulo);
    output_age_matrix ( mat_name, modulo);

    fp=vfopen  ( mat_list_name,"a");
    fprintf ( fp, "%s\n", mat_name);
    vfclose ( fp);

    S=seq2year (S,modulo);
    A=seq2aln (S, NULL, KEEP_GAP);
    output_fasta_seq (tmp=vtmpnam (NULL),A);
    file_cat ( tmp, seq_list);

    if ( display==0)
      {
	display_output_filename ( stdout, "AGE_MAT_LIST", "MAT_LIST", mat_list_name, CHECK);
	display_output_filename ( stdout, "AGE_SEQ", "FASTA", seq_list, CHECK);
	display=1;
      }
    fprintf ( stderr, "\nModulo:%d years", modulo);
    fprintf ( stderr, "\n");
    free_aln (A);
    return S;
  }
//
// Name MAnipulation
//

Alignment *clean_aln (Alignment *A)
{
  if ( A)
    {
      A->seq_comment=clean_string (A->nseq, A->seq_comment);
      A->aln_comment=clean_string (A->nseq, A->aln_comment);
      A->name=translate_names(A->nseq, A->name);
      (A->S)=clean_sequence ((A->S));
    }
  return A;
}
Sequence *clean_sequence ( Sequence *S)
{
  if ( !S) return S;

  S->seq_comment=clean_string (S->nseq, S->seq_comment);
  S->name=translate_names(S->nseq, S->name);
  return S;
}
char ** translate_names (int n, char **name)
{
  int a;
  for ( a=0; a<n; a++)
    name[a]=translate_name(name[a]);
  return name;
}
char * translate_name ( char *name)
	{

	int len;
	int a;
	char buf[1000];

	len=strlen (name);

	//if ( name[0]=='\'')return name;

	for ( a=0; a<len; a++)
		{
		if ( isspace(name[a]))name[a]='\0';
		else if (strchr (";(),:#><'�", name[a]))name[a]='_';

		}
	sprintf (buf,"%s",decode_name (name, DECODE));
	if ( strlen (buf)>read_array_size_new ((char *)name))
	  {
	    name=(char*)vrealloc (name, sizeof (char)*(strlen (buf)+1));
	  }
	sprintf (name, "%s", buf);

	return name;
	}
char *decode_name (char *name, int mode)
{
  static char ***name_list;
  static int n;
  static char tag[100];
  int a;

  if (mode==CLEAN)
    {
      for (a=0; a<n; a++)
	{
	  vfree (name_list[a][0]);
	  vfree (name_list[a][1]);
	  vfree (name_list[a]);
	}
      vfree (name_list);
      tag[0]='\0';
    }

  //spacial modes
  if ( mode == CODELIST)
    {
      char *file;
      file=vtmpnam (NULL);
      for (a=0; a< n; a++)
	printf_file(file, "a", "#CODE: %s <=> %s\n", name_list[a][0], name_list[a][1]);
      return file;
    }
  if (mode ==DECODE && name_list==NULL)return name;
  if ( name==NULL) return name;



  if (!tag[0])
    {
      vsrand (0);
      sprintf ( tag, "TCTAG_%d",rand ()%100000);
    }

  if ( mode==CODE)
    {
      for (a=0; a< n; a++)
	if ( strm (name, name_list[a][0]))return name_list[a][1];

      name_list=(char***)realloc (name_list, sizeof (char**)*(n+1));
      name_list[n]=(char**)vcalloc (2, sizeof (char*));
      name_list[n][0]=(char*)vcalloc (strlen (name)+1, sizeof (char));
      name_list[n][1]=(char*)vcalloc (100, sizeof (char));
      sprintf ( name_list[n][0], "%s", name);
      sprintf ( name_list[n][1], "%s_%d", tag,n+1);
      return name_list[n++][1];
    }
  else if ( mode ==DECODE)
    {
      char *p;
      int i;
      if ( !(p=after_strstr (name, tag)))return name;
      else
	{
	  sscanf (p, "_%d", &i);
	  return name_list[i-1][0];
	}
    }
  else
    {
      printf_exit (EXIT_FAILURE, stderr,"Unknown Mode for Decode_name [FATAL:%s]", PROGRAM);
    }
  return NULL;
}


FILE * display_sequences_names (Sequence *S, FILE *fp, int check_pdb_status, int print_templates)
        {
	    int a;
	    int max_len;
	    char *r;

	    if ( !S)
	       {
		   fprintf (fp,"\nERROR: NO SEQUENCE READ [FATAL:%s]\n", PROGRAM); myexit (EXIT_FAILURE);
	       }
	    for ( a=0, max_len=0; a< S->nseq; a++)max_len=MAX(max_len, strlen (S->name[a]));
	    fprintf ( fp, "\nINPUT SEQUENCES: %d SEQUENCES  [%s]", S->nseq,(S->type)?S->type:"Unknown type");
	    if (S->nseq<MAX_NSEQ_4_DISPLAY)
	      {
		for ( a=0; a< S->nseq; a++)
		  {
		    fprintf (fp, "\n  Input File %-*s Seq %-*s Length %4d type %s",max_len,S->file[a], max_len,S->name[a],(int)strlen ( S->seq[a]), S->type);
		    if (check_pdb_status)
		      {
			if ((r=seq_is_pdb_struc (S, a)))fprintf (fp, " Struct Yes PDBID %s", get_pdb_id(r));
			else fprintf (fp, " Struct No");
		      }
		    else fprintf (fp, " Struct Unchecked");
		    if ( print_templates)
		      fp=display_sequence_templates (S, a, fp);
		  }
		fprintf ( fp, "\n");
	      }

	    return fp;

	}
Sequence *add_file2file_list (char *name, Sequence *S)
{

  if (!S) S=declare_sequence (1,1,10);
  else S=realloc_sequence   (S,S->nseq+1,0);S->nseq=0;

  sprintf ( S->name[S->nseq++], "%s", name);
  return S;

}




