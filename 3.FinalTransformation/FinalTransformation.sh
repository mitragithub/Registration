#! /opt/hpc/bin/bash

if (($# < 1))
then
    echo "No input brainname argument... exiting"
    exit 1
fi
#define input
BrainNO=$1

#define environment
source /sonas-hs/it/hpc/home/easybuild/lmod-setup.sh
module load foss/2016a
module load IntelPython/2.7.12
MATLABCMD="/opt/hpc/pkg/MATLAB/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -r "

#define path
PIPELINE_DIR=/sonas-hs/mitra/hpc/home/xli/RegistrationPipelineV3/
OUTPUT_DIR=$PIPELINE_DIR/Data/$BrainNO/Transformation_OUTPUT/
ATLAS_DIR=$PIPELINE_DIR/ATLAS/
CODE_DIR=$PIPELINE_DIR/Bin/3.FinalTransformation/

#################################start process: target to registered space##############################

#1 transform high resolution images
img_path=/nfs/mitraweb2/mnt/disk127/main/mba_converted_imaging_data/*/$BrainNO/
# img_path=/nfs/mitraweb2/mnt/disk127/main/mba_converted_imaging_data/MD721\&720/MD720/
recon_path=$PIPELINE_DIR/Data/$BrainNO/Registration_OUTPUT/
csv_path=$PIPELINE_DIR/Data/$BrainNO/INPUT_DATA/
mkdir -p $OUTPUT_DIR/reg_high_tif/
$MATLABCMD "cd('$CODE_DIR'); maxNumCompThreads(2); transform('$img_path', '$recon_path', '$OUTPUT_DIR/reg_high_tif/'); exit"

#padding the tif
mkdir -p $OUTPUT_DIR/reg_high_tif_pad/
mkdir -p $OUTPUT_DIR/reg_high_tif_pad_jp2/
$MATLABCMD "cd('$CODE_DIR'); maxNumCompThreads(2); padtif('$OUTPUT_DIR/reg_high_tif/', '$OUTPUT_DIR/reg_high_tif_pad/', '$OUTPUT_DIR/reg_high_tif_pad_jp2/');exit"

# thumbnail
sh $CODE_DIR/kdujpg.sh $OUTPUT_DIR/reg_high_tif_pad_jp2/
rm -rf $OUTPUT_DIR/reg_high_tif_pad_jp2/*.tif

# ################################start process: atlas to registered space################################
#1 transform atlas to registered space in 5um space
mkdir -p $OUTPUT_DIR/low_seg/
$MATLABCMD "cd('$CODE_DIR'); transform_seg('$ATLAS_DIR/annotation_50.vtk', '$PIPELINE_DIR/Data/$BrainNO/INPUT_DATA/', '$recon_path','$OUTPUT_DIR/low_seg/', 5); exit"

mkdir -p $OUTPUT_DIR/reg_high_seg/
$MATLABCMD "cd('$CODE_DIR'); maxNumCompThreads(2); segresize('$OUTPUT_DIR/low_seg/', '$OUTPUT_DIR/reg_high_tif/', '$OUTPUT_DIR/reg_high_seg/'); exit"


mkdir -p $OUTPUT_DIR/reg_high_seg_pad/
$MATLABCMD "cd('$CODE_DIR'); maxNumCompThreads(2); padseg('$OUTPUT_DIR/reg_high_seg/', '$OUTPUT_DIR/reg_high_seg_pad/'); done('$BrainNO'); exit"

cd $CODE_DIR/makejson/
python brain_region.py $BrainNO

cd $CODE_DIR
python update_database.py $BrainNO