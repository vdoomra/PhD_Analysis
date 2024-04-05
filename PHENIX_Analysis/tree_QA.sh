#!/bin/sh

LIST=`ls -lhtr /phenix/plhf/vdoomra/taxi/Run15pp200CAERTP108/18940/data/4*.root | awk '{printf("%s\n",$9)}'`

NUM=0

export MYINSTALL=/direct/phenix+u/vdoomra/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

INPUT=$(( $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR

cp /phenix/plhf/vdoomra/tree_QA_C.so .
cp /phenix/plhf/vdoomra/field_map.root .
for file in $LIST
do
  if (( $NUM == $1 ))
  then

    NAME=`echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    output=QA_$NAME
    echo $file $output

    root -l -b <<EOF
    gSystem->Load("libDileptonAnalysisEvent");
    gSystem->Load("tree_QA_C.so");
    tree_QA("$file","$output")
EOF
 
mv $output /phenix/plhf/vdoomra/QA_output/
  fi
  NUM=$(( $NUM + 1 ))
done

popd
rm -r $DIR
