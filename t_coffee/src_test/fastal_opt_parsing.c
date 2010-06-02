
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "fastal_lib_header.h"


void
tree_parse(struct Fastal_arguments *arguments, char* param)
{
	char delims[] = "";
	arguments->tree_method = strtok(param,delims);
	if (arguments->tree_method == "parttree")
	{
		arguments->tree_param1 = 6;
		arguments->tree_param2 = 150;
	}
	char *tmp = strtok(NULL,delims);
	if (tmp != NULL)
	{
		arguments->tree_param1 = atoi(tmp);
		arguments->tree_param2 = atoi(strtok(NULL,delims));
	}
// 	printf("A: %s %i %i", arguments->tree_method, arguments->tree_param1, arguments->tree_param2);
}



void
arg_parse (int argc, char **argv, struct Fastal_arguments *arguments)
{
// 	default values

	arguments->diag_method = "blast";
	arguments->output_file = "out->aln";
	arguments->tree_file = NULL;
	arguments->gep = -1;
	arguments->gop = -10;
	arguments->method = "fast";
	arguments->tree_method = "oligotree";
	arguments->mat="dna_idmat";
	arguments->tree_only = 0;
	arguments->evaluate = 0;
	arguments->score = 0;
// 	arguments->retree = 0;
	arguments->agreement_score = 0;
	arguments->num_ref_aln = 0;
	arguments->is_dna = 1;
	arguments->aln_ref = NULL;
	arguments->aln2test = NULL;
	arguments->tree_out = 0;
	arguments->gap_iterate = 0;


	int i = 1;
	char *param;

	for (i = 0; i < argc; ++i)
	{
		param = argv[i];
		if ( ( !strcmp(param, "-d" )) || ( !strcmp(param, "--is_dna")))
		{
			arguments->is_dna = 1;
			arguments->mat = "dna_idmat";
		}
		else if ( ( !strcmp(param, "-a" )) || ( !strcmp(param, "--is_aa")))
		{
			arguments->is_dna = 0;
			arguments->mat = "blosum62mt";

		}
// 		printf("%s\n", arguments->mat);
	}
	if (arguments->is_dna)
	{
		arguments->tree_param1 = 2;
		arguments->tree_param2 = 5;
	}
	else
	{
		arguments->tree_param1 = 1;
		arguments->tree_param2 = 10;
	}


	i = 1;
	while (i < argc)
	{
		param = argv[i];
		if ( ( !strcmp(param, "-i" )) || ( !strcmp(param, "--in")))
		{
			arguments->sequence_file = argv[++i];
		}
		else if ( ( !strcmp(param, "-t" )) || ( !strcmp(param, "--tree_file")))
		{
			arguments->tree_file = argv[++i];
		}
		else if ( !strcmp(param, "--mat"))
		{
			arguments->mat = argv[++i];
		}
		else if ( !strcmp(param, "--tree_method"))
		{
			tree_parse(arguments, argv[++i]);
// 			arguments->tree_file = argv[++i];
		}
		else if ( ( !strcmp(param, "-o" )) || ( !strcmp(param, "--outfile")))
		{
			arguments->output_file = argv[++i];
		}
		else if ( ( !strcmp(param, "-m" )) || ( !strcmp(param, "--method")))
		{
			++i;
			if ( (!strcmp(argv[i], "fast")) || (!strcmp(argv[i], "nw")) || (!strcmp(argv[i], "gotoh")) ||  (!strcmp(argv[i], "udisc")) )
				arguments->method = argv[i];
			else
			{
				printf("Method %s unknown\n", argv[i]);
				exit(1);
			}
		}
		else if ( ( !strcmp(param, "-b" )) || ( !strcmp(param, "--diag_method")))
		{
			++i;
			if ( (!strcmp(argv[i], "blast")) || (!strcmp(argv[i], "blastz")) || (!strcmp(argv[i], "blat")) || (!strcmp(argv[i], "ktup")))
				arguments->diag_method = argv[i];
			else
			{
				printf("DIAG Method %s unknown\n", argv[i]);
				exit(1);
			}
		}
		else if ( ( !strcmp(param, "-d" )) || ( !strcmp(param, "--is_dna")))
		{
			arguments->is_dna = 1;
		}
		else if ( ( !strcmp(param, "-a" )) || ( !strcmp(param, "--is_aa")))
		{
			arguments->is_dna = 0;
		}
		else if ( ( !strcmp(param, "-g" )) || ( !strcmp(param, "--gop")))
		{
			arguments->gop = atof(argv[++i]);
		}
// 		else if ( ( !strcmp(param, "-r" )) || ( !strcmp(param, "--retree")))
// 		{
// 			arguments->retree = atoi(argv[++i]);
// 		}
		else if ( ( !strcmp(param, "-e" )) || ( !strcmp(param, "--gep")))
		{
			arguments->gep = atof(argv[++i]);
		}
		else if ( ( !strcmp(param, "-k" )) || ( !strcmp(param, "--p1")))
		{
			arguments->tree_param1 = atoi(argv[++i]);
		}
		else if ( ( !strcmp(param, "-c" )) || ( !strcmp(param, "--p2")))
		{
			arguments->tree_param2 = atoi(argv[++i]);
		}
		else if ( !strcmp(param, "--eval_aln"))
		{
			arguments->evaluate = 1;
		}
		else if ( ( !strcmp(param, "-s" )) || ( !strcmp(param, "--score")))
		{
			arguments->score = 1;
		}
		else if ( !strcmp(param, "--tree_out"))
		{
			arguments->tree_out = 1;
		}
		else if ( !strcmp(param, "--tree_only"))
		{
			arguments->tree_only = 1;
		}
		else if ( !strcmp(param, "--agreement_score"))
		{
			arguments->agreement_score = 1;
		}
		else if ( !strcmp(param, "--make_ref_aln"))
		{
			arguments->num_ref_aln = atoi(argv[++i]);
			arguments->num_seq_in_ref = atoi(argv[++i]);
		}
		else if ( !strcmp(param, "--aln_ref"))
		{
			arguments->aln_ref = argv[++i];
		}
		else if ( !strcmp(param, "--gap_iterate"))
		{
			arguments->gap_iterate = atoi(argv[++i]);
		}
		else if ( !strcmp(param, "--aln2test"))
		{
			arguments->aln2test = argv[++i];
		}
		else if ( ( !strcmp(param, "-h" )) || ( !strcmp(param, "--help")) || ( !strcmp(param, "-?" )))
		{
			printf("Fastal - a fast alignment tool\n");
			printf("-i    --in              The sequence_file\n");
			printf("-t    --tree_file       The treefile\n");
			printf("      --tree_method     The method to produce the tree:\n");
			printf("                        oligo - very fast method\n");
			printf("                        parttree - method like in mafft\n");
			printf("-o    --outfile         The output file\n");
			printf("-m    --method          The method to use:\n");
			printf("                        fast    - fast approximate algorithm [default]\n");
			printf("                        nw      - needleman-wunsch\n");
			printf("                        gotoh   - gotoh\n");
			printf("-d    --is_dna          Sequence are DNA\n");
			printf("-a    --is_aa           Sequences are amino acids\n");
			printf("-g    --gop             Gap opening costs\n");
			printf("-e    --gep             Gap extension costs\n");
			printf("-r    --retree          Number of reconstructions [default=0]\n");
			printf("      --tree_only       only tree is produced\n");
			printf("      --eval_aln        calculates only the sum-of-pairs-score of a given alignment\n");
			printf("-s    --score           calculates the sum of paris score after construction of an alignment. Can take a long time. [default: off]\n");
// 			printf("      --make_ref_aln    calculates the sum of paris score after construction of an alignment. Can take a long time. [default: off]\n");
			exit(0);
		}
		else
		{
			printf("Argument %s unknown\n", param);
			exit(1);
		}
		++i;
	}
}
