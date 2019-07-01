filename=$1
output_dir=$2

    outfile="$(basename $filename)"
    outfile=${outfile/.tif/.jp2}
    # outfile=${outfile/"$filename"/}
    echo $outfile
    kdu_compress -i $filename -o "$output_dir/${outfile}" -rate 1 Creversible=yes Clevels=7 Clayers=8 Stiles=\{1024,1024\} Corder=RPCL Cuse_sop=yes ORGgen_plt=yes ORGtparts=R Cblk=\{32,32\} -num_threads 1

rm -rf $filename