#FILE:    /home/suhre/bin/fugue_align.sh
#AUTHOR:  Karsten Suhre
#DATE:    Wed Dec 3 17:01:56 CET 2003
#PURPOSE: simulate fugue_pair
#BUGS:    
#MODIF:   

PDB=${2?"usage: $0 fastafile pdbfile"}
SEQ=$1
PDB=`basename $PDB .pdb`
PDBFILE="$PDB.pdb"

SEQFILE=`basename $SEQ`

if [ ! -f $SEQ ] ; then
  echo "ERROR: file $SEQ not found"
  exit 1
fi

if [ ! -f $PDBFILE ] ; then
  echo "ERROR: file $PDBFILE not found"
  exit 1
fi

if grep -q '^ATOM' $PDBFILE ; then
  echo pdbfile OK > /dev/null
else
  echo "ERROR: $PDBFILE does not contain any ATOM records"
  exit 1
fi

export FUGUE_ROOT=/home/igs/Tools/T-COFFEE_web/Struct
PATH=$FUGUE_ROOT/joy-5.10:$PATH
PATH=$FUGUE_ROOT/joy_related-1.03:$PATH
PATH=$FUGUE_ROOT/joy_related-1.03/sstruc/:$PATH
PATH=$FUGUE_ROOT/joy_related-1.03/hbond:$PATH
PATH=$FUGUE_ROOT/fugue:$PATH
PATH=$FUGUE_ROOT/run_fugue:$PATH

export FUGUE_LIB_LIST=.
export HOMSTRAD_PATH=.
export HOMS_PATH=.

export MELODY_CLASSDEF=$FUGUE_ROOT/fugue/classdef.dat
export MELODY_SUBST=$FUGUE_ROOT/fugue/allmat.dat

mkdir /tmp/$$.fugue
cp $SEQ /tmp/$$.fugue
cp $PDBFILE /tmp/$$.fugue
cd /tmp/$$.fugue

#/home/igs/Tools/T-COFFEE_web/Struct/magicPDB.pl -i $PDBFILE  -split -o toto > /dev/null 2> /dev/null
#if [ -s toto_A.pdb ] ; then
#    mv toto_A.pdb $PDBFILE
#fi

joy $PDBFILE > /dev/null 2> /dev/null
melody -t $PDB.tem > /dev/null 2> /dev/null
fugueali -print -seq $SEQFILE -prf $PDB.fug | sed '/^ *$/d'

cd ..
rm -r /tmp/$$.fugue
