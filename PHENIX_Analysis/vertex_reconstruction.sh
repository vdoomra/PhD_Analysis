#!/bin/sh

LIST=`ls -lhtr /phenix/plhf/vdoomra/taxi/Run15pp200CAMBP108/18981/data/4*.root | awk '{printf("%s\n",$9)}'`

NUM=0
APPLY_DEADMAP=1
MODE=0

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

INPUT=$(( $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR

cp /phenix/plhf/vdoomra/vertex_reconstruction_C.so .
cp /phenix/plhf/vdoomra/field_map.root .

for file in $LIST
do
  if (( $NUM == $1 ))
  then

    NAME=`echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    RUNNO=`echo $NAME | awk -F \. '{printf("%s\n",$1)}'`
    output=tree_$NAME
    output_text=/phenix/plhf/vdoomra/VTXDeadMap/deadmap_files/out_$RUNNO.txt
    echo $file $output $output_text

    root -l -b <<EOF
    gSystem->Load("/phenix/u/vdoomra/install/lib/libDileptonAnalysisEvent");
    gSystem->Load("vertex_reconstruction_C.so")
    vertex_reconstruction("$file","$output", "$output_text",$APPLY_DEADMAP,$MODE)
EOF
 
mv $output /phenix/plhf/vdoomra/vertex_reconstruction_output/pp200_mb_onlyL1/
  fi
  NUM=$(( $NUM + 1 ))
done

popd
rm -r $DIR
