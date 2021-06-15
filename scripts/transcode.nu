#!/usr/bin/env nu

# Encode/decode list of files
def transcode [
    ...rest  # files to be transcoded
] {
    ./transcoder $rest ../test
    echo $rest | each {
        let inp = (path dirname -r '../test' | path parse | update extension astc | path join)
        let out = (path dirname -r '~/git/cpc/astc-encoder-dissection/test' | path expand | path parse | update extension png | path join)
        # let inp = (path dirname -r '../test' | path extension -r astc)
        # let out = (path dirname -r ~/git/cpc/astc-encoder-dissection/test | path extension -r png)
        echo [[inp out]; [$inp $out]]
        astcenc-dissect -dl $inp $out
    }
}

# Encode/decode all files in a directory
def transcode-dir [
    inp-dir: path  # directory with files to be encoded
] {
    let dir_basename = ($inp-dir | path basename)

    let astc_out_dir = ([ '../test' $dir_basename ] | path join)
    let png_out_dir = ([ '~/git/cpc/astc-encoder-dissection/test' $dir_basename ] | path join | path expand)

    let inp_files = (ls $inp-dir).name

    mkdir $astc_out_dir
    ./transcoder $inp_files $astc_out_dir

    mkdir $png_out_dir
    for $inp_file in $inp_files {
        let astc_file = ($inp_file | path parse | update extension astc | update parent $astc_out_dir | path join)
        let png_file = ($astc_file | path parse | update extension png | update parent $png_out_dir | path join)

        astcenc -dl $astc_file $png_file
    }

    source ~/git/scripts/conda_eval.nu
    conda-run common ~/git/scripts/calculate_psnr.py ~/data/combined_set "--list" "../combined_dir.txt" "-o" "../test"
}

# Decode all .astc files in a directory
def decode-dir [
    inp-dir: path  # directory with files to be encoded
] {
    let dir_basename = ($inp-dir | path basename)

    let png_out_dir = ([ '~/git/cpc/astc-encoder-dissection/test' $dir_basename ] | path join | path expand)

    let inp_files = (ls $"($inp-dir)/*.astc").name

    mkdir $png_out_dir
    for $inp_file in $inp_files {
        let png_file = ($inp_file | path parse | update extension png | update parent $png_out_dir | path join)

        astcenc -dl $inp_file $png_file
    }

    source ~/git/scripts/conda_eval.nu
    conda-run common ~/git/scripts/calculate_psnr.py ~/data/combined_set "--list" "../combined_dir.txt" "-o" "../test"
}
