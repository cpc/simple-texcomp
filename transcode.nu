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
        astcenc -dl $inp $out
    }
}

# Encode all files in a directory
def encode-dir [
    inp-dir: path        # directory with files to be encoded
    --out-dir (-o): path # optional output directory (default is '../test')
] {
    let dir-basename = ($inp-dir | path basename)
    let inp-files = (ls $inp-dir).name

    let out-dir = (if ($out-dir | empty?) { '../test' } { $out-dir })
    let out-subdir = ([ $out-dir $dir-basename ] | path join)

    mkdir $out-subdir
    ./transcoder $inp-files $out-subdir
}

def calculate-psnr [
    ref-dir: path        # reference directory
    list-file: path  # file with directories to be evaluated
] {
    source ~/git/extern/nu_scripts/virtual_environments/conda.nu
    load-env (conda-env common)

    python ~/git/scripts/calculate_psnr.py $ref-dir '--list' $list-file '-o' '~/git/scripts/results/psnr_new'
}

# Encode/decode all files in a directory
def transcode-dir [
    inp-dir: path         # directory with files to be encoded
    --astc-dir (-a): path # optional astc output directory (default is '../test')
    --png-dir (-p): path  # optional png output directory (default is '~/git/cpc/astc-encoder-dissection/test')
    --decode              # whether to decode the astc files
] {
    let dir-basename = ($inp-dir | path basename)

    let astc-dir = (if ($astc-dir | empty?) { '../test' } { $astc-dir } | path expand)
    let astc-subdir = ([ $astc-dir $dir-basename ] | path join | path expand)

    let png-dir = (if ($png-dir | empty?) { '~/git/cpc/astc-encoder-dissection/test' } { $png-dir })
    let png-subdir = ([ $png-dir $dir-basename ] | path join)

    # let astc-out-dir = ([ '../test' $dir-basename ] | path join)
    # let png-out-dir = ([ '~/git/cpc/astc-encoder-dissection/test' $dir-basename ] | path join | path expand)

    let inp-files = (ls ($inp-dir | path expand)).name

    mkdir $astc-subdir
    ./transcoder $inp-files $astc-subdir

    if ($decode | empty?) {
        echo "--- Not decoding, finished ---"
    } {
        mkdir $png-subdir
        for $inp-file in $inp-files {
            let astc-file = ($inp-file | path parse | update extension astc | update parent $astc-subdir | path join)
            let png-file = ($astc-file | path parse | update extension png | update parent $png-subdir | path join)

            astcenc -dl $astc-file $png-file
        }

        if "combined_set" in ($inp-dir | into string) {
            source ~/git/extern/nu_scripts/virtual_environments/conda.nu
            load-env (conda-env common)
            ~/git/scripts/calculate_psnr.py ~/data/combined_set "--list" "../combined_set.txt" "-o" "../test"
        } { }
    }
}

# Decode all .astc files in a directory
def decode-dir [
    inp-dir: path  # directory with files to be encoded
] {
    let dir-basename = ($inp-dir | path basename)

    let png-out-dir = ([ '~/git/cpc/astc-encoder-dissection/test' $dir-basename ] | path join | path expand)

    let inp-files = (ls $"($inp-dir)/*.astc").name

    mkdir $png-out-dir
    for $inp-file in $inp-files {
        let png-file = ($inp-file | path parse | update extension png | update parent $png-out-dir | path join)

        astcenc -dl $inp-file $png-file
    }

    source ~/git/extern/nu_scripts/virtual_environments/conda.nu
    load-env (conda-env common)
    ~/git/scripts/calculate_psnr.py ~/data/combined_set "--list" "../combined_set.txt" "-o" "../test"
}