#!/bin/bash
# conda activate pyscenic

## environment
ref_dir=SCENIC_REF/cisTarget_database 
tfs=$ref_dir/allTFs_mm.txt
feather1=$ref_dir/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
feather2=$ref_dir/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=$ref_dir/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl

outdir=out/SCENIC_withCTL
input_loom=$outdir/sample.loom
grn_output=$outdir/grn_output.tsv
ctx_output=$outdir/ctx_output.csv
aucell_output=$outdir/aucell_output.loom
ls $tfs  $feather1 $feather2 $tbl  


## transform count matrix to loom format
python scripts\utils\csv2loom.py \
    $outdir/for.scenic.count.ctl-tall.csv \
    $input_loom

## step1 grn
arboreto_with_multiprocessing.py \
    $input_loom \
    $tfs \
    --method grnboost2 \
    --output $grn_output \
    --num_workers 80 \
    --seed 777

## step2  ctx
pyscenic ctx \
  $grn_output \
  $feather1 \
  --annotations_fname $tbl \
  --expression_mtx_fname $input_loom \
  --mode custom_multiprocessing \
  --output $ctx_output \
  --num_workers 80 \
  --mask_dropouts


## AUCell
pyscenic aucell \
  $input_loom  \
  $ctx_output \
  --output $aucell_output \
  --num_workers 80