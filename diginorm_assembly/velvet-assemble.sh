#! /bin/bash

K=$1
filename=$2
scriptpath=`dirname $0`
workdir=`dirname $filename`
BASE=`basename $filename`


if [ \! -f $BASE.se -o \! -f $BASE.pe \]; then
   python $scriptpath/strip-and-split-for-assembly.py $filename $workdir/$BASE
fi

if [ \! -s $BASE.se -o \! -s $BASE.pe ]; then
   echo 'WARNING -- one or more sequences files are EMPTY; may fail'
fi

velveth $workdir/$BASE.ass.$K $K -fasta -short $workdir/${BASE}.se -shortPaired $workdir/${BASE}.pe && \
velvetg $workdir/$BASE.ass.$K -read_trkg yes -exp_cov auto -cov_cutoff auto -scaffolding no

rm $workdir/${BASE}.ass.$K/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
