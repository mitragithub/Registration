#! /opt/hpc/bin/bash

if (($# < 1))
then
    echo "No input brainname argument... exiting"
    exit 1
fi
#define input
BrainNO=$1
# thick_dz=$2
# thick_sections=$3 #manully input, should be changed to database later, !FIXME

#define environment
source /sonas-hs/it/hpc/home/easybuild/lmod-setup.sh
module load foss/2016a
module load IntelPython/2.7.12
MATLABCMD="/opt/hpc/pkg/MATLAB/R2018a/bin/matlab -nodisplay -nodesktop -nosplash -r "

#define path
PIPELINE_DIR=/sonas-hs/mitra/hpc/home/xli/RegistrationPipelineV3/
OUTPUT_DIR=$PIPELINE_DIR/Data/$BrainNO/Registration_OUTPUT/
ATLAS_DIR=$PIPELINE_DIR/ATLAS/
CODE_DIR=$PIPELINE_DIR/Bin/2.Registration/

#start process
mkdir -p $OUTPUT_DIR
$MATLABCMD "cd ${CODE_DIR}; atlas_free_rigid_alignment('${PIPELINE_DIR}/Data/${BrainNO}/INPUT_DATA/', '*.tif*', '${OUTPUT_DIR}', 10, [32,16,8,4],500); done('$BrainNO'); exit"

python $CODE_DIR/update_database.py $BrainNO