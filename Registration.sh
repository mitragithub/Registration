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
$MATLABCMD "cd ${CODE_DIR}; registration_pipeline_rnaseq_atlas_free_align('${ATLAS_DIR}/ara_nissl_50.vtk', '${PIPELINE_DIR}/Data/${BrainNO}/INPUT_DATA/', 'rnaseq_config.ini', '${OUTPUT_DIR}', '${OUTPUT_DIR}', '${ATLAS_DIR}/annotation_50.vtk'); exit"
#_atlas_free_align
#qsub -pe threads 16 -l m_mem_free=2G Registration.sh MD718