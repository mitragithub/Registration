#


PIPELINE_DIR=/sonas-hs/mitra/hpc/home/xli/snRNA_registration_pipeline/

brainid=$1
img_path=/nfs/mitraweb2/mnt/disk125/main/mba_converted_imaging_data/*/${brainid}/
recon_path=/sonas-hs/mitra/hpc/home/xli/snRNA_registration_pipeline/RegisteredData/cis.jhu.edu/~dtward/share/for_xu/${brainid}*/

OUTPUT_DIR=$PIPELINE_DIR/FinalTransformation/$brainid
mkdir $OUTPUT_DIR

MATLABCMD="/opt/hpc/pkg/MATLAB/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -r "



mkdir -p $OUTPUT_DIR/reg_high_tif/
$MATLABCMD "cd('$PIPELINE_DIR/FinalTransformation'); transform('$img_path', '$recon_path', '$OUTPUT_DIR/reg_high_tif/');exit"

mkdir -p $OUTPUT_DIR/reg_high_tif_pad/
mkdir -p $OUTPUT_DIR/reg_high_tif_pad_jp2/
$MATLABCMD "cd('$PIPELINE_DIR/FinalTransformation'); padtif('$OUTPUT_DIR/reg_high_tif/', '$OUTPUT_DIR/reg_high_tif_pad/', '$OUTPUT_DIR/reg_high_tif_pad_jp2/');exit"



    # sh $PIPELINE_DIR/FinalTransformation/kducomp.sh $OUTPUT_DIR/reg_high_tif_pad/ $OUTPUT_DIR/reg_high_tif_pad_jp2/


sh $PIPELINE_DIR/FinalTransformation/kdujpg.sh $OUTPUT_DIR/reg_high_tif_pad_jp2/

rm -rf $OUTPUT_DIR/reg_high_tif_pad/*.tif

mkdir -p $OUTPUT_DIR/lowseg/
$MATLABCMD "cd('$PIPELINE_DIR/FinalTransformation'); transform_seg('$recon_path', '$recon_path', '$OUTPUT_DIR/lowseg/');exit"

mkdir -p $OUTPUT_DIR/reg_high_seg/
$MATLABCMD "cd('$PIPELINE_DIR/FinalTransformation'); segresize('$OUTPUT_DIR/lowseg/', '$OUTPUT_DIR/reg_high_tif/', '$OUTPUT_DIR/reg_high_seg/');exit"

mkdir -p $OUTPUT_DIR/reg_high_seg_pad/
$MATLABCMD "cd('$PIPELINE_DIR/FinalTransformation'); padseg('$OUTPUT_DIR/reg_high_seg/', '$OUTPUT_DIR/reg_high_seg_pad/');exit"

cd $PIPELINE_DIR/FinalTransformation/makejson/
python brain_region.py $brainid




