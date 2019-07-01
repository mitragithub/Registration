input_dir=$1
output_dir=$2

# IMGDIR="/scratch/783UPDATE/reg_high_tif_pad_jp2/"
for filename in ${input_dir}*.jp2; do
    outfile=${filename/.jp2/}
    outfile=${outfile/"$input_dir"/}
    echo $outfile
    IMGFILE=$outfile
    kdu_expand -i "${input_dir}${IMGFILE}.jp2" -o "${input_dir}${IMGFILE}.tif" -reduce 6 -num_threads 8
    convert "${input_dir}${IMGFILE}.tif" -depth 8 "${input_dir}${IMGFILE}.jpg"
done
