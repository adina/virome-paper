!/bin/bash

PREFIX=$1
SEQS=${*:2}
SCRIPTPATH=~ubuntu/khmer/scripts
SCRIPTPATH2=~ubuntu/khmer/sandbox

K=20
HASHBITS_SIZE=50e9
N_TABLES=4
C=20

echo python $SCRIPTPATH/normalize-by-median.py -k $K -C $C -N $N_TABLES -x $HASHBITS_SIZE -s $PREFIX.diginorm1.kh -R $PREFIX.diginorm1.report $SEQS

python $SCRIPTPATH/normalize-by-median.py -k $K -C $C -N $N_TABLES -x $HASHBITS_SIZE -s $PREFIX.diginorm1.kh -R $PREFIX.diginorm1.report $SEQS

echo python $SCRIPTPATH/normalize-by-median.py -p -k $K -C $C -N $N_TABLES -x $HASHBITS_SIZE -s $PREFIX.diginorm2.kh --loadhash $PREFIX.diginorm1.kh -R $PREFIX.diginorm2.report $SEQS

python $SCRIPTPATH/normalize-by-median.py -p -k $K -C $C -N $N_TABLES -x $HASHBITS_SIZE -s $PREFIX.diginorm2.kh --loadhash $PREFIX.diginorm1.kh -R $PREFIX.diginorm2.report $SEQS

echo python $SCRIPTPATH/filter-abund.py -V $PREFIX.diginorm2.kh *keep

python $SCRIPTPATH/filter-abund.py -V $PREFIX.diginorm2.kh *keep

python $SCRIPTPATH/load-into-count.py partition.ht 

python $SCRIPTPATH/partition-graph.py --subset-size $SUBSET_SIZE $BASENAME
python $SCRIPTPATH/merge-partitions.py -k $K $BASENAME
python $SCRIPTPATH/annotate-partitions.py -k $K $BASENAME $SEQS
python  $SCRIPTPATH/extract-partitions.py $BASENAME `basename ${SEQS}`.part
