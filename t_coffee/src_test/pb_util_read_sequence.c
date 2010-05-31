/********* Sequence input routines for CLUSTAL W *******************/
/* DES was here.  FEB. 1994 */
/* Now reads PILEUP/MSF and CLUSTAL alignment files */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

/*
*	Prototypes
*/

extern Boolean linetype(char *,char *);
extern Boolean blankline(char *);
extern void warning(char *,...);
extern void error(char *,...);
extern char *	rtrim(char *);
extern char *	blank_to_(char *);
extern void 	getstr(char *,char *);

void fill_chartab(void);


static void         get_seq(char *,char *,int *,char *);
static void get_clustal_seq(char *,char *,int *,char *,int);
static void     get_msf_seq(char *,char *,int *,char *,int);
static void check_infile(int *);




static int count_clustal_seqs(void);
static int count_msf_seqs(void);

/*
 *	Global variables
 */

static FILE *fin;


char *amino_acid_codes   =    "ABCDEFGHIKLMNPQRSTUVWXYZ-";  /* DES */
char *nucleic_acid_order = 	  "ACGTUN";
static int seqFormat;
static char chartab[128];

void fill_chartab(void)	/* Create translation and check table */
{
	register int i;
	register int c;


	for(i=0;i<128;chartab[i++]=0);
	for(i=0,c=0;c<=amino_acid_codes[i];i++)
		chartab[c]=chartab[tolower(c)]=c;

}

static void get_msf_seq(char *sname,char *seq,int *len,char *tit,int seqno)
/* read the seqno_th. sequence from a PILEUP multiple alignment file */
{
	static char *line;
	int i,j,k;
	unsigned char c;
	if ( !line)line=vcalloc ( (MAXLINE+1), sizeof (char));

	fseek(fin,0,0); 		/* start at the beginning */

	*len=0;				/* initialise length to zero */
        for(i=0;;i++) {
		if(fgets(line,MAXLINE+1,fin)==NULL) return; /* read the title*/
		if(linetype(line,"/") ) break;		    /* lines...ignore*/
	}

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) {

			for(i=1;i<seqno;i++) fgets(line,MAXLINE+1,fin);
                        for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;
			for(k=j;k<=strlen(line);k++) if(line[k] == ' ') break;
			strncpy(sname,line+j,MIN(MAXNAMES,k-j)); 
			sname[MIN(MAXNAMES,k-j)]=EOS;
			rtrim(sname);
                       	blank_to_(sname);

			for(i=k;*len < SEQ_MAX_LEN;i++) {
				c=line[i];
				if(c == '.') c = '-';
				if(c == '*') c = 'X';
				if(c == '\n' || c == EOS) break; /* EOL */
				if( (c=chartab[c])) seq[++(*len)]=c;
			}
			if(*len == SEQ_MAX_LEN) return;

			for(i=0;;i++) {
				if(fgets(line,MAXLINE+1,fin)==NULL) return;
				if(blankline(line)) break;
			}
		}
	}
}


static void get_clustal_seq(char *sname,char *seq,int *len,char *tit,int seqno)
/* read the seqno_th. sequence from a clustal multiple alignment file */
{
	static char *line;
	int i,j;
	unsigned char c;
	if ( !line)line=vcalloc ( (MAXLINE+1), sizeof (char));
	fseek(fin,0,0); 		/* start at the beginning */

	*len=0;				/* initialise length to zero */
	fgets(line,MAXLINE+1,fin);	/* read the title line...ignore it */

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) {

			for(i=1;i<seqno;i++) fgets(line,MAXLINE+1,fin);
                        for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;
			strncpy(sname,line+j,MAXNAMES-j); /* remember entryname */
			sname[MAXNAMES]=EOS;
			rtrim(sname);
                       	blank_to_(sname);

			for(i=14;*len < SEQ_MAX_LEN;i++) {
				c=line[i];
				if(c == '\n' || c == EOS) break; /* EOL */
				if( (c=chartab[c])) seq[++(*len)]=c;
			}
			if(*len == SEQ_MAX_LEN) return;

			for(i=0;;i++) {
				if(fgets(line,MAXLINE+1,fin)==NULL) return;
				if(blankline(line)) break;
			}
		}
	}
}


static void get_seq(char *sname,char *seq,int *len,char *tit)
{
	static char *line;
	int i, offset, c=0;

	if ( !line)line=vcalloc ( (MAXLINE+1), sizeof (char));
	switch(seqFormat) {

/************************************/
		case EMBLSWISS:
			while( !linetype(line,"ID") )
				fgets(line,MAXLINE+1,fin);
			
                        for(i=5;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
			sname[MAXNAMES]=EOS;
			rtrim(sname);
                        blank_to_(sname);
			fprintf ( stderr, "\n%s", sname);
			while( !linetype(line,"SQ") )
				fgets(line,MAXLINE+1,fin);
			
			*len=0;
			while(fgets(line,MAXLINE+1,fin)) {
				for(i=0;*len < SEQ_MAX_LEN;i++) {
					c=line[i];
				if(c == '\n' || c == EOS || c == '/')
					break;			/* EOL */
				if( (c=chartab[c]))
					seq[++(*len)]=c;
				}
			if(*len == SEQ_MAX_LEN || c == '/') break;
			}
		break;
		
/************************************/
		case PIR:
			while(*line != '>')
				{
				
				fgets(line,MAXLINE+1,fin);
				 }			
                        for(i=4;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
			sname[MAXNAMES]=EOS;
			rtrim(sname);
                        blank_to_(sname);

			fgets(line,MAXLINE+1,fin);
			strncpy(tit,line,MAXTITLES);
			tit[MAXTITLES]=EOS;
			i=strlen(tit);
			if(tit[i-1]=='\n') tit[i-1]=EOS;
			
			*len=0;
			while(fgets(line,MAXLINE+1,fin)) {
				
				for(i=0;*len < SEQ_MAX_LEN;i++) 
				   {
				   c=line[i];
				   if(c == '\n' || c == EOS || c == '*')
					break;			/* EOL */
			
				   if( (c=chartab[c]))
					seq[++(*len)]=c;

				   }
			if(*len == SEQ_MAX_LEN || c == '*') break;
			}
		break;
/***********************************************/
	case (PEARSON):
			
			while(*line != '>')
				{
				fgets(line,MAXLINE+1,fin);
				}
                        for(i=1;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
			sname[MAXNAMES]=EOS;
			rtrim(sname);
                        blank_to_(sname);

			*tit=EOS;
			
			*len=0;
			while(fgets(line,MAXLINE+1,fin)) 
				{
				for(i=0;*len < SEQ_MAX_LEN;i++) 
					{
					c=line[i];
					if(c == '\n' || c == EOS || c == '>')
						break;			/* EOL */
						
					if( (c=chartab[c]))
						{seq[++(*len)]=c;
					         }
					}
				if(*len == SEQ_MAX_LEN || c == '>') break;
				}
			break;
/**********************************************/
		case GDE:
		
			while(*line != '#' ||*line != '%' )
					fgets(line,MAXLINE+1,fin);
			
			
			
			
			for (i=1;i<=MAXNAMES;i++) {
				if (line[i] == '(' || line[i] == '\n')
                                  {
                                    i--;
                                    break;
                                  }
				sname[i-1] = line[i];
			}
			sname[i]=EOS;
			offset=0;
			if (sname[i-1] == '(') sscanf(&line[i],"%d",&offset);
			else offset = 0;
			for(i=MAXNAMES-1;i > 0;i--) 
				if(isspace(sname[i])) {
					sname[i]=EOS;	
					break;
				}		
                        blank_to_(sname);


			*tit=EOS;
			
			*len=0;
			for (i=0;i<offset;i++) seq[++(*len)] = '-';
			while(fgets(line,MAXLINE+1,fin)) {
			if(*line == '%' || *line == '#' || *line == '"') break;
				for(i=0;*len < SEQ_MAX_LEN;i++) {
					c=line[i];
				if(c == '\n' || c == EOS) 
					break;			/* EOL */
			
				if( (c=chartab[c]))
					seq[++(*len)]=c;
				}
			if(*len == SEQ_MAX_LEN) break;
			}
		break;
/***********************************************/
	}
	
	if(*len == SEQ_MAX_LEN)
		warning("Sequence %s truncated to %d residues",
				sname,(pint)SEQ_MAX_LEN);
				
	seq[*len+1]=EOS;
}


int readseqs(char *saga_file,char ***SAGA_SEQ, char*** SAGA_NAMES, int ***SAGA_LEN) /*first_seq is the #no. of the first seq. to read */
{
	
	static char *line;
	static char *seq1;
	static char *sname1;
	static char *title;
	int i,j,no_seqs;
	static int l1;
	int a;
	int b, l;
	int first_seq=0;
	
	if ( !line) line=vcalloc   ( (FILENAMELEN+1), sizeof (char));
	if ( !seq1) seq1=vcalloc   ( (SEQ_MAX_LEN+2), sizeof (char));
	if ( !sname1)sname1=vcalloc ( (MAXNAMES+1), sizeof (char));
	if ( !title)title=vcalloc  ( (MAXTITLES+1), sizeof (char));
	
	fill_chartab();
	
	fin=vfopen(saga_file,"r");
	
	
	no_seqs=0;
	check_infile(&no_seqs);

	if(no_seqs == 0)
		return 0;       /* return the number of seqs. (zero here)*/
	
	SAGA_SEQ[0]= vcalloc ( no_seqs, sizeof ( char*));
	SAGA_NAMES[0]= vcalloc ( no_seqs, sizeof ( char*));
	SAGA_LEN[0]= declare_int ( no_seqs,3);
	
	
	for(i=first_seq;i<=first_seq+no_seqs-1;i++) 
		{    /* get the seqs now*/
		if(seqFormat == CLUSTAL) 
			get_clustal_seq(sname1,seq1,&l1,title,i-first_seq+1);
		if(seqFormat == MSF)
			    get_msf_seq(sname1,seq1,&l1,title,i-first_seq+1);
		else
			get_seq(sname1,seq1,&l1,title);
		
		if(l1 > SEQ_MAX_LEN) 
			{
			error("Sequence too long. Maximum is %d",(pint)SEQ_MAX_LEN);
			return 0;       /* also return zero if too many */
			}
		
		
		
		for ( a=0; a<l1; a++)
			seq1[a]=seq1[a+1];	
		seq1[l1]='\0';
		
		(SAGA_SEQ[0])[i]=vcalloc ( strlen ( seq1)+1, sizeof ( char));
		(SAGA_NAMES[0])[i]=vcalloc ( strlen ( sname1)+1, sizeof ( char));
		
		l1=strlen (seq1);
		(SAGA_LEN[0])[i][0]=l1;
		(SAGA_SEQ[0])[i][l1]='\0';
		
		sprintf ( (SAGA_SEQ[0])[i], "%s", seq1);
		sprintf ( (SAGA_NAMES[0])[i], "%s", sname1);
		
		l=strlen (SAGA_NAMES[0][i]);
		for ( b=0; b< l; b++)
		  {
		    if ( isspace(SAGA_NAMES[0][i][b])){SAGA_NAMES[0][i][b]='\0';break;}
		  }
			  
		
	}

	fclose(fin);
/*
   JULIE
   check sequence names are all different - otherwise phylip tree is 
   confused.
*/
	for(i=first_seq;i<=first_seq+no_seqs-1;i++) {
		for(j=i+1;j<=first_seq+no_seqs-1;j++) {
			if (strncmp((SAGA_NAMES[0])[i],(SAGA_NAMES[0])[j],MAXNAMES) == 0) {
				error("Multiple sequences found with same name (first %d chars are significant)", MAXNAMES);
				return -1;
			}
		}
	}
			
	return no_seqs;    /* return the number of seqs. read in this call */
}




static void check_infile(int *nseqs)
{
	static char *line;
	int i;	

	if ( !line)line=vcalloc ( (MAXLINE+1), sizeof (char));
	
	for ( i=0; i<=MAXLINE; i++)line[i]='a';
	*nseqs=0;
	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if ((*line != '\n') && (*line != ' ') && (*line != '\t'))
			break;
	}
        
	for(i=0;i<=6;i++) line[i] = toupper(line[i]);

	if( linetype(line,"ID") ) {					/* EMBL/Swiss-Prot format ? */
		seqFormat=EMBLSWISS;
		(*nseqs)++;
	}
        else if( linetype(line,"CLUSTAL") ) {
		seqFormat=CLUSTAL;
	}
        else if( linetype(line,"PILEUP") ) {
		seqFormat = MSF;
	}
	else if(*line == '>') {						/* no */
		seqFormat=(line[3] == ';')?PIR:PEARSON; /* distinguish PIR and Pearson */
		(*nseqs)++;
	}
	else if((*line == '"') || (*line == '%') || (*line == '#')) {
		seqFormat=GDE; /* GDE format */
		if (*line == '%') {
                        (*nseqs)++;
			
		}
		else if (*line == '#') {
			(*nseqs)++;
			
		}
	}
	else {
		seqFormat=UNKNOWN;
		return;
	}

	while(fgets(line,MAXLINE+1,fin) != NULL) {
		switch(seqFormat) {
			case EMBLSWISS:
				if( linetype(line,"ID") )
					(*nseqs)++;
				break;
			case PIR:
			case PEARSON:
				if( *line == '>' )
					(*nseqs)++;
				break;
			case GDE:
				if(( *line == '%' ) )
					(*nseqs)++;
				else if (( *line == '#') )
					(*nseqs)++;
				break;
			case CLUSTAL:
				*nseqs = count_clustal_seqs();
/* DES */ 			/* fprintf(stdout,"\nnseqs = %d\n",(pint)*nseqs); */
				fseek(fin,0,0);
				return;
				break;
			case MSF:
				*nseqs = count_msf_seqs();
				fseek(fin,0,0);
				return;
				break;
			case USER:
			default:
				break;
		}
	}
	fseek(fin,0,0);
}


static int count_clustal_seqs(void)
/* count the number of sequences in a clustal alignment file */
{
	static char *line;
	int  nseqs;

	if ( !line)line=vcalloc ( (MAXLINE+1), sizeof (char));

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) break;		/* Look for next non- */
	}						/* blank line */
	nseqs = 1;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(blankline(line)) return nseqs;
		nseqs++;
	}

	return 0;	/* if you got to here-funny format/no seqs.*/
}

static int count_msf_seqs(void)
{
/* count the number of sequences in a PILEUP alignment file */

	static char *line;
	int  nseqs;

	if ( !line)line=vcalloc ( (MAXLINE+1), sizeof (char));

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(linetype(line,"/")) break;
	}

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) break;		/* Look for next non- */
	}						/* blank line */
	nseqs = 1;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(blankline(line)) return nseqs;
		nseqs++;
	}

	return 0;	/* if you got to here-funny format/no seqs.*/
}



