#!/bin/bash

[ -z "$CBLEN" ] && CBLEN=16
[ -z "$UMILEN" ] && UMILEN=12

setup_alevin() {
    [ -n "$alevin_setup" ] && return
    export AF_SAMPLE_DIR=$indices_dir/alevin
    export ALEVIN_FRY_HOME="$AF_SAMPLE_DIR/af_home"
    mkdir -p $AF_SAMPLE_DIR
    mkdir -p $ALEVIN_FRY_HOME
    simpleaf set-paths
    alevin_setup=1
}

run_star() {
    if [ -n "$starcmd" ]; then
        chmod +x $starcmd
        $starcmd
        return $?       
    fi
    for fastq_dir in "${fastq_array[@]}"; do
        [ -z "$overwrite" ] && [ -f "$outputdir/$(basename $fastq_dir)/star/Solo.out/GeneFull/Summary.csv" ] && echo "results exist - skipping $fastq_dir" && continue
        reads1=$(find -L ${fastq_dir} -name "*$R1pattern*.$fastq_suffix" -type f | sort | tr '\n' ',' | sed 's/,$//')
        reads2=$(find -L ${fastq_dir} -name "*$R2pattern*"$fastq_suffi -type f | sort | tr '\n' ',' | sed 's/,$//')
        echo "$reads1"
        echo "$reads2"
        tmpDir="$outputdir/$(basename $fastq_dir)/star/tmp"
        [ -d "$tmpDir" ] && rm -rf "$tmpDir"
        mkdir -p "$outputdir/$(basename $fastq_dir)/star"
       cmd="STAR --runThreadN $nthreads \
        --outTmpDir $tmpDir \
        --soloType CB_UMI_Simple \
        --soloCBlen $CBLEN \
        --soloUMIlen $UMILEN \
        --soloUMIstart $((CBLEN+1)) \
        --soloCBstart 1 \
        --soloBarcodeReadLength 0 \
        --soloCBwhitelist $whitelist \
        --genomeDir $indices_dir/star \
        --limitIObufferSize 50000000 50000000 \
        --outSJtype None \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingBinsN 500 \
        --outBAMsortingThreadN 6 \
        --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCellFilter EmptyDrops_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --soloFeatures Gene GeneFull Velocyto \
        --outFileNamePrefix $outputdir/$(basename $fastq_dir)/star/ \
        --readFilesIn $reads2 $reads1 \
        --outFilterMatchNmin 28 \
        --outFilterMatchNminOverLread 0 \
        --outFilterMultimapNmax 1 \
        --alignEndsType Local \
        --outFilterScoreMinOverLread 0"

       [ -n "$gzipped_fastq" ] && cmd="$cmd --readFilesCommand zcat"
        echo "$cmd"
        eval $cmd
    done

}

run_kallisto() {
    if [ -n "$kallistocmd" ]; then
        chmod +x $kallistocmd
        $kallistocmd
        return $?       
    fi
    for fastq_dir in "${fastq_array[@]}"; do
        [ -z "$overwrite" ] && [ -f "$outputdir/$(basename $fastq_dir)/kallisto/counts_unfiltered/cells_x_genes.total.mtx" ] && echo "results exist - skipping $fastq_dir" && continue
        echo "working on $fastq_dir"
        mkdir -p $outputdir/$(basename $fastq_dir)/kallisto
        [ -d "$outputdir/$(basename $fastq_dir)/kallisto/tmp" ] && rm -rf "$outputdir/$(basename $fastq_dir)/kallisto/tmp"
        reads1=( $(find -L ${fastq_dir} -name "*$R1pattern*.$fastq_suffix" -type f | sort | tr '\n' ' ') )
        reads2=( $(find -L ${fastq_dir} -name "*$R2pattern*.$fastq_suffix" -type f | sort | tr '\n' ' ') )
        reads=""
        for i in "${!reads1[@]}"; do
            reads="$reads ${reads1[$i]} ${reads2[$i]}" 
        done
        echo "reads: $reads"
        cmd="kb count -x 10XV3 \
        --workflow nac --sum=total -i $indices_dir/kallisto/nac_offlist_1/index.idx \
        -g $indices_dir/kallisto/nac_offlist_1/g \
        -c1 $indices_dir/kallisto/nac_offlist_1/c1 \
        -c2 $indices_dir/kallisto/nac_offlist_1/c2 \
        -o $outputdir/$(basename $fastq_dir)/kallisto/ --overwrite --verbose \
        -t $nthreads \
        $reads"
        echo "$cmd"
        eval $cmd
    done
}
run_piscem() {
   setup_alevin
   if [ -n "$piscemcmd" ]; then
        chmod +x $piscemcmd
        $piscemcmd
        return $?       
    fi
    for fastq_dir in "${fastq_array[@]}"; do
        echo "working on $fastq_dir"
        [ -n "useSplicei" ] && mkdir -p $outputdir/$(basename $fastq_dir)/piscem/splicei
        [ -n "useSpliceu" ] && mkdir -p $outputdir/$(basename $fastq_dir)/piscem/spliceu

        reads1="$(find -L $fastq_dir -name "*$R1pattern*.$fastq_suffix" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
        reads2="$(find -L $fastq_dir -name "*$R2pattern*.$fastq_suffix" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
        if [ -n "$useSplicei" ]; then
            [  -z "$overwrite" ] && [ -f "$outputdir/$(basename $fastq_dir)/piscem/splicei/af_quant/alevin/quants_mat.mtx" ] && continue
            cmd="simpleaf quant \
            --use-piscem \
            --reads1 $reads1 \
            --reads2 $reads2 \
            --threads $nthreads \
            --chemistry 10xv3 --resolution cr-like \
            --expected-ori fw --unfiltered-pl \
            --index $indices_dir/piscem/splicei/index \
            --t2g-map $indices_dir/piscem/splicei/index/t2g_3col.tsv \
            --output $outputdir/$(basename $fastq_dir)/piscem/splicei/ " 
            echo "$cmd"
            eval $cmd
        fi
        if [ -n "$useSpliceu" ]; then
            [  -z "$overwrite" ] && [ -f "$outputdir/$(basename $fastq_dir)/piscem/spliceu/af_quant/alevin/quants_mat.mtx" ] && echo "results exist - skipping $fastq_dir" && continue
            cmd="simpleaf quant \
            --use-piscem \
            --reads1 $reads1 \
            --reads2 $reads2 \
            --threads $nthreads \
            --chemistry 10xv3 --resolution cr-like \
            --expected-ori fw --unfiltered-pl \
            --index $indices_dir/piscem/spliceu/index \
            --t2g-map $indices_dir/piscem/spliceu/index/t2g_3col.tsv \
            --output $outputdir/$(basename $fastq_dir)/piscem/spliceu/ " 
            echo "$cmd"
            eval $cmd
        fi
    done 
}

run_salmon() {
    setup_alevin
    if [ -n "$salmoncmd" ]; then
        chmod +x $salmoncmd
        $salmoncmd
        return $?       
    fi
    for fastq_dir in "${fastq_array[@]}"; do
        echo "working on $fastq_dir"
        [ -n "useSplicei" ] && mkdir -p $outputdir/$(basename $fastq_dir)/salmon/splicei
        [ -n "useSpliceu" ] && mkdir -p $outputdir/$(basename $fastq_dir)/salmon/spliceu

        reads1="$(find -L $fastq_dir -name "*$R1pattern*.$fastq_suffix" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
        reads2="$(find -L $fastq_dir -name "*$R2pattern*.$fastq_suffix" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
        if [ -n "$useSplicei" ]; then
            [  -z "$overwrite" ] && [ -f "$outputdir/$(basename $fastq_dir)/salmon/splicei/af_quant/alevin/quants_mat.mtx" ] && echo "results exist - skipping $fastq_dir" && continue
            cmd="simpleaf quant \
            --no-piscem \
            --reads1 $reads1 \
            --reads2 $reads2 \
            --threads $nthreads \
            --chemistry 10xv3 --resolution cr-like \
            --expected-ori fw --unfiltered-pl \
            --index $indices_dir/salmon/splicei/index \
            --t2g-map $indices_dir/salmon/splicei/index/t2g_3col.tsv \
            --output $outputdir/$(basename $fastq_dir)/salmon/splicei/ " 
            echo "$cmd"
            eval $cmd
        fi
        if [ -n "$useSpliceu" ]; then
            [  -z "$overwrite" ] && [ -f "$outputdir/$(basename $fastq_dir)/salmon/spliceu/af_quant/alevin/quants_mat.mtx" ] && echo "results exist - skipping $fastq_dir" && continue
            cmd="simpleaf quant \
            --no-piscem \
            --reads1 $reads1 \
            --reads2 $reads2 \
            --threads $nthreads \
            --chemistry 10xv3 --resolution cr-like \
            --expected-ori fw --unfiltered-pl \
            --index $indices_dir/salmon/spliceu/index \
            --t2g-map $indices_dir/salmon/spliceu/index/t2g_3col.tsv \
            --output $outputdir/$(basename $fastq_dir)/salmon/spliceu/ " 
            echo "$cmd"
            eval $cmd
        fi
    done 
}
parse_string_into_array(){
    local string=$1
    local -n array=$2
    # Remove square brackets and initial and final double quotes and split on commas
    IFS=',' read -r -a array <<< $(echo "$string" | sed 's/[][]//g; s/\"//g ')
    IFS=' '
}
translate_cmd(){
    local cmd=$1
    local -n array=$2
    IFS=',' read -r -a array <<< $(echo "$cmd" | sed 's/[][]//g; s/\"//g ')
    IFS=' '
}
#print out key pairs of env variables
print_env(){
    for key in "${!env[@]}"; do
        echo "$key=${env[$key]}"
    done
}

[ -z "$useStar" ] && [ -z "$useKallisto" ] && [ -z "$useSalmon" ] && [ -z "$usePiscem" ] && echo "Must choose one aligner" && exit 1
if [ -n "$usePiscem" ] || [ -n "$useSalmon" ]; then
    [ -z "$useSplicei" ] && [ -z "$useSpliceu" ] && echo "Must choose one reference for salmon/piscem" && exit 1
fi
# check if fastq_suffix ends in .gz
[[ "$fastq_suffix" != *.gz ]] && [ -n "$gzipped_fastq" ] && fastq_suffix="$fastq_suffix.gz"


# parse fastqdirs into an array
parse_string_into_array "$fastqdirs" fastq_array

[ -n "$useStar" ] && run_star || true
[ -n "$useKallisto" ] && run_kallisto || true
[ -n "$useSalmon" ] && run_salmon || true
[ -n "$usePiscem" ] && run_piscem  || true 

