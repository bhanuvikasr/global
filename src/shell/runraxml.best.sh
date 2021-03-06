#!/bin/bash
set -x

id=$1 # This is gene ID
C=$2  # number of CPUs
H=$3  # Directory where genes are located
in=$4 # name of input alignment
N=$5  # Number of replicates
model="GTRGAMMA" # DNA Model to use
T=`echo $model|sed -e "s/\(.*\)/\L\1/g"`
dirn=raxmlboot.$T # output raxml directory

part=""
pp="unpart"
[ $# -gt 5 ] && [ $6 != "-" ] && part="-M -q $6"
[ $# -gt 5 ] && [ $6 != "-" ] && pp="part"
dirn=$dirn.$pp

out=""
echo Rooting at $out

cd $H/$id/

mkdir logs

mkdir -p $dirn
cd $dirn

if [ ! -s $in.phylip ] || [ "`head -n1 $in.phylip`" == "0 0" ] ; then
 $WS_HOME/global/src/shell/convert_to_phylip.sh ../$in $in.phylip
 rm $in.phylip.reduced
fi
[ $# == 6 ] && ln -s ../$6 .

# In case this is a re-run, rename outputs from the old runs
tar cfz `mktemp ./BEST.backup.XXXXX.tgz` --remove-files *best*

# Estimate the RAxML best tree
if [ $C -gt 1 ]; then
 raxmlHPC-8.0.19-PTHREADS-SSE3-modified $out -m $model -T $C -n best -s $in.phylip -N $N -p $RANDOM $part &>$H/$id/logs/best
else
 raxmlHPC-8.0.19-SSE3-modified $out -m $model -n best -s $in.phylip -N $N -p $RANDOM $part &>$H/$id/logs/best
fi

#/share/home/01721/smirarab/bin/mapsequences.py raxml/RAxML_bipartitions.ml namemap ml.mapped -rev &>logs/map

cd ..
if [ -s $H/$id/$dirn/RAxML_bestTree.best ]; then 
  echo "Done">.done.raxml.$T.$pp.1
  cd $dirn
  tar cfj best.run.logs.tar.bz --remove-files RAxML_log* RAxML_pars* RAxML_res*
  rm *reduced
else 
  echo Sorry, Failed!
fi
