#!/usr/bin/sh


script="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/sshic/pipeline.py"

fragments="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt"
oligos="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/inputs/capture_oligo_positions.csv"
centromeres="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna//data/samples/inputs/S288c_chr_centro_coordinates.tsv"
binning="1000 2000 3000 5000 10000 20000 40000 50000 80000 100000"
wt4h="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/inputs/wt4h_pcrdupkept.tsv"
wt2h="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/inputs/wt2h_pcrdupkept.tsv"
ws_centros=150000
ws_telos=150000
excluded_chr="chr2 chr3 chr5 2_micron mitochondrion chr_artificial"

# AD403
sample="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD403_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
reference=$wt2h
python3 $script -s $sample \
               -f $fragments \
               -o $oligos \
               -c $centromeres \
               -r $reference \
               -b $binning \
               --window-size-centros $ws_centros \
               --window-size-telos $ws_telos \
               --excluded-chr $excluded_chr \
               --inter-norm

## AD404
sample="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD404_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
reference=$wt4h
python3 $script -s $sample \
               -f $fragments \
               -o $oligos \
               -c $centromeres \
               -r $reference \
               -b $binning \
               --window-size-centros $ws_centros \
               --window-size-telos $ws_telos \
               --excluded-chr $excluded_chr \
               --inter-norm
# AD401
sample="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD401_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
reference="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD403/AD403_global_statistics.tsv"
python3 $script -s $sample \
               -f $fragments \
               -o $oligos \
               -c $centromeres \
               -r $reference \
               -b $binning \
               --window-size-centros $ws_centros \
               --window-size-telos $ws_telos \
               --excluded-chr $excluded_chr \
               --inter-norm

# AD402
sample="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD402_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
reference="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD404/AD404_global_statistics.tsv"
python3 $script -s $sample \
               -f $fragments \
               -o $oligos \
               -c $centromeres \
               -r $reference \
               -b $binning \
               --window-size-centros $ws_centros \
               --window-size-telos $ws_telos \
               --excluded-chr $excluded_chr \
               --inter-norm
# AD405
sample="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD405_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
reference="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD404/AD404_global_statistics.tsv"
python3 $script -s $sample \
               -f $fragments \
               -o $oligos \
               -c $centromeres \
               -r $reference \
               -b $binning \
               --window-size-centros $ws_centros \
               --window-size-telos $ws_telos \
               --excluded-chr $excluded_chr \
               --inter-norm
# AD406
sample="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD406_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
reference="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD403/AD403_global_statistics.tsv"
python3 $script -s $sample \
               -f $fragments \
               -o $oligos \
               -c $centromeres \
               -r $reference \
               -b $binning \
               --window-size-centros $ws_centros \
               --window-size-telos $ws_telos \
               --excluded-chr $excluded_chr \
               --inter-norm


# AD407
sample="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD407_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q20.txt"
reference="/home/nicolas/Documents/Projects/ssHiC/hic_ssdna/data/samples/pcrdupkept/AD404/AD404_global_statistics.tsv"
python3 $script -s $sample \
               -f $fragments \
               -o $oligos \
               -c $centromeres \
               -r $reference \
               -b $binning \
               --window-size-centros $ws_centros \
               --window-size-telos $ws_telos \
               --excluded-chr $excluded_chr \
               --inter-norm



