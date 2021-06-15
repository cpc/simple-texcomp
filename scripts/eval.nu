#!/usr/bin/env nu

# meson compile
./transcoder ~/data/combined_set/* ../test/combined

ls ~/git/cpc/simple-texcomp/test/combined/* | get name | each { let out = (echo $it | path parse | update parent ~/git/cpc/astc-encoder-dissection/test/combined_simple | update extension png | path join); astcenc-dissect -dl $it $out }

source ~/git/scripts/conda_eval.nu
conda-run common ~/git/scripts/calculate_psnr.py '~/data/combined_set' '--list' '~/git/scripts/combined_dir.txt' '-o' '~/git/scripts/results/psnr_new'
