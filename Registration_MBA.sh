#!/bin/bash

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
#source /sonas-hs/it/hpc/home/easybuild/lmod-setup.sh
#module load foss/2016a
#module load IntelPython/2.7.12
MATLABCMD="matlab -nodisplay -nodesktop -nosplash -r "

#define path
PIPELINE_DIR=/grid/mitra/home/xli/RegistrationPipelineV3/
#define path and run step1, data collection.
OUTPUT_DIR=$PIPELINE_DIR/Data/$BrainNO/INPUT_DATA/
CODE_DIR=$PIPELINE_DIR/Bin/1.DataPreparation/

mkdir -p $OUTPUT_DIR
cp /mitranfs/*/main/mba_converted_imaging_data/*/${BrainNO}/*.tif $OUTPUT_DIR/
$MATLABCMD "cd ${CODE_DIR}; create_location_csv_MBA('${OUTPUT_DIR}',14.72, 14.72, 20,'${OUTPUT_DIR}'); exit"


OUTPUT_DIR=$PIPELINE_DIR/Data/$BrainNO/Registration_OUTPUT/
ATLAS_DIR=$PIPELINE_DIR/ATLAS/
CODE_DIR=$PIPELINE_DIR/Bin/2.Registration/

#start process
mkdir -p $OUTPUT_DIR
$MATLABCMD "cd ${CODE_DIR}; maxNumCompThreads(16);registration_pipeline_MBA('${ATLAS_DIR}/ara_nissl_50.vtk', '${PIPELINE_DIR}/Data/${BrainNO}/INPUT_DATA/', 'mba_nissl_fluoro_config.ini', '${OUTPUT_DIR}', '${OUTPUT_DIR}', '${ATLAS_DIR}/annotation_50.vtk'); done('$BrainNO'); exit"

# python $CODE_DIR/update_database.py $BrainNO
python $CODE_DIR/update_database_qcpipeline.py $BrainNO
