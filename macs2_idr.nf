#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run main.nf --samples samples.txt --outdir outdir

  """.stripIndent()
  }


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * VALIDATE
 */

if (params.samples)     { ch_samples = Channel.fromPath(params.samples, checkIfExists: true) } else { exit 1, 'Samples not specified' }
    ch_samples
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.sample_id,  file(row.bam1), file(row.bam1_ctrl), file(row.bam2), file(row.bam2_ctrl)] }
        .into { ch_samples_split_macs; ch_samples_split_pooled_pseudo; ch_samples_split_self_pseudo; ch_samples_split_bigwig}


println ("""
        ===========================================================================================
                                            Macs2 peakcallig & IDR
        ===========================================================================================
        Sample file: ${params.samples}
        Macs2 q-value: ${params.macs_q}
        IDR threshold: ${params.idr_threshold}
        ===========================================================================================
        """)



/*
 * 1. Peak calling for true replicates
 */
process PEAKCALLING {
    publishDir "${params.outdir}/${sample_id}/peaks", mode: 'copy', pattern: '*_sort_peaks.narrowPeak'

    input:
    set val(sample_id), path(bam1), path(bam1_ctrl), path(bam2), path(bam2_ctrl) from ch_samples_split_macs

    output:
    tuple val(sample_id), path("${sample_id}_rep1_sort_peaks.narrowPeak"), path("${sample_id}_rep2_sort_peaks.narrowPeak")  into ch_peak_reps_true, ch_peak_reps_self_pseudo


    script:
    """
    # peak calling for replicates
    macs2 callpeak -t ${bam1} -c ${bam1_ctrl} -f BAM -g ${params.genome_size} -n ${sample_id}_rep1 -B -q ${params.macs_q}  2> ${sample_id}_rep1_macs2.log
    macs2 callpeak -t ${bam2} -c ${bam2_ctrl} -f BAM -g ${params.genome_size} -n ${sample_id}_rep2 -B -q ${params.macs_q}  2> ${sample_id}_rep2_macs2.log

    #Sort peak by -log10(p-value)
    sort -k8,8nr ${sample_id}_rep1_peaks.narrowPeak > ${sample_id}_rep1_sort_peaks.narrowPeak
    sort -k8,8nr ${sample_id}_rep2_peaks.narrowPeak > ${sample_id}_rep2_sort_peaks.narrowPeak
    """
}

/*
 * 2. Create self-pseudoreps
 */
process CREATE_SELF_PSEUDOREPS {

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(bam1), path(bam1_ctrl), path(bam2), path(bam2_ctrl) from ch_samples_split_self_pseudo

    output:
    tuple val(sample_id), path("${sample_id}_rep1_self_pseudorep1.bam"), path("${sample_id}_ctrl_rep1_self_pseudorep1.bam"), path("${sample_id}_rep1_self_pseudorep2.bam"),path("${sample_id}_ctrl_rep1_self_pseudorep2.bam"), path("${sample_id}_rep2_self_pseudorep1.bam"), path("${sample_id}_ctrl_rep2_self_pseudorep1.bam"), path("${sample_id}_rep2_self_pseudorep2.bam"),path("${sample_id}_ctrl_rep2_self_pseudorep2.bam")  into ch_bam_self_pseudoreps

    script:
    """
    ## Rep1
    samtools view -H ${bam1} > ${sample_id}_rep1_header.sam

    # Split merged treatments
    nlines=\$(samtools view ${bam1} | wc -l ) # Number of reads in the BAM file
    nlines=\$(( (nlines + 1) / 2 )) # half that number
    samtools view ${bam1} | shuf - | split -d -l \${nlines} - "${sample_id}_rep1" # Shuffle lines in file and split into two SAM files
    cat ${sample_id}_rep1_header.sam ${sample_id}_rep100 | samtools view -bS - > ${sample_id}_rep1_self_pseudorep1.bam
    cat ${sample_id}_rep1_header.sam ${sample_id}_rep101 | samtools view -bS - > ${sample_id}_rep1_self_pseudorep2.bam

    #Split merged treatment BAM
    nlines=\$(samtools view ${bam1_ctrl} | wc -l ) # Number of reads in the BAM file
    nlines=\$(( (nlines + 1) / 2 )) # half that number
    samtools view ${bam1_ctrl} | shuf - | split -d -l \${nlines} - "${sample_id}_ctrl_rep1" # This will shuffle the lines in the file and split in two
    cat ${sample_id}_rep1_header.sam ${sample_id}_ctrl_rep100 | samtools view -bS - > ${sample_id}_ctrl_rep1_self_pseudorep1.bam
    cat ${sample_id}_rep1_header.sam ${sample_id}_ctrl_rep101 | samtools view -bS - > ${sample_id}_ctrl_rep1_self_pseudorep2.bam


    ## Rep2
    samtools view -H ${bam2} > ${sample_id}_rep2_header.sam

    # Split merged treatments
    nlines=\$(samtools view ${bam2} | wc -l ) # Number of reads in the BAM file
    nlines=\$(( (nlines + 1) / 2 )) # half that number
    samtools view ${bam2} | shuf - | split -d -l \${nlines} - "${sample_id}_rep2" # Shuffle lines in file and split into two SAM files
    cat ${sample_id}_rep2_header.sam ${sample_id}_rep200 | samtools view -bS - > ${sample_id}_rep2_self_pseudorep1.bam
    cat ${sample_id}_rep2_header.sam ${sample_id}_rep201 | samtools view -bS - > ${sample_id}_rep2_self_pseudorep2.bam

    #Split merged treatment BAM
    nlines=\$(samtools view ${bam2_ctrl} | wc -l ) # Number of reads in the BAM file
    nlines=\$(( (nlines + 1) / 2 )) # half that number
    samtools view ${bam2_ctrl} | shuf - | split -d -l \${nlines} - "${sample_id}_ctrl_rep2" # This will shuffle the lines in the file and split in two
    cat ${sample_id}_rep2_header.sam ${sample_id}_ctrl_rep200 | samtools view -bS - > ${sample_id}_ctrl_rep2_self_pseudorep1.bam
    cat ${sample_id}_rep2_header.sam ${sample_id}_ctrl_rep201 | samtools view -bS - > ${sample_id}_ctrl_rep2_self_pseudorep2.bam
    """
}

/*
 * 3. Create pooled & pooled pseudo repliates
 */
process CREATE_POOLED_PSEUDOREPS {

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(bam1), path(bam1_ctrl), path(bam2), path(bam2_ctrl) from ch_samples_split_pooled_pseudo

    output:
    tuple val(sample_id), path("${sample_id}_pooled_merged.bam"), path("${sample_id}_pooled_ctrl_merged.bam") into ch_bam_pooled_reps, ch_bam_pooled_reps_for_bigwig
    tuple val(sample_id), path("${sample_id}_pooled_pseudorep1.bam"), path("${sample_id}_ctrl_pooled_pseudorep1.bam"), path("${sample_id}_pooled_pseudorep2.bam"),path("${sample_id}_ctrl_pooled_pseudorep2.bam")  into ch_bam_pooled_pseudoreps

    script:
    """
    # Merge bam
    samtools merge -u ${sample_id}_pooled_merged.bam ${bam1} ${bam2}
    samtools view -H ${sample_id}_pooled_merged.bam > ${sample_id}_header.sam

    #Split merged treatments
    nlines=\$(samtools view ${sample_id}_pooled_merged.bam | wc -l ) # Number of reads in the BAM file
    nlines=\$(( (nlines + 1) / 2 )) # half that number
    samtools view ${sample_id}_pooled_merged.bam | shuf - | split -d -l \${nlines} - "${sample_id}" # Shuffle lines in file and split into two SAM files
    cat ${sample_id}_header.sam ${sample_id}00 | samtools view -bS - > ${sample_id}_pooled_pseudorep1.bam
    cat ${sample_id}_header.sam ${sample_id}01 | samtools view -bS - > ${sample_id}_pooled_pseudorep2.bam

    # Merge ctrl bams
    samtools merge -u ${sample_id}_pooled_ctrl_merged.bam ${bam1_ctrl} ${bam2_ctrl}

    #Split merged treatment BAM
    nlines=\$(samtools view ${sample_id}_pooled_ctrl_merged.bam | wc -l ) # Number of reads in the BAM file
    nlines=\$(( (nlines + 1) / 2 )) # half that number
    samtools view ${sample_id}_pooled_ctrl_merged.bam | shuf - | split -d -l \${nlines} - "${sample_id}_ctrl" # This will shuffle the lines in the file and split in two
    cat ${sample_id}_header.sam ${sample_id}_ctrl00 | samtools view -bS - > ${sample_id}_ctrl_pooled_pseudorep1.bam
    cat ${sample_id}_header.sam ${sample_id}_ctrl01 | samtools view -bS - > ${sample_id}_ctrl_pooled_pseudorep2.bam
    """
}

/*
 * 4. Peak calling for pooled replicates
 */
process PEAKCALLING_POOLED_REPS {

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(bam_pr), path(bam_pr_ctrl)  from ch_bam_pooled_reps

    output:
    tuple val(sample_id), path("${sample_id}_pooled_sort_peaks.narrowPeak")  into ch_peak_pooled_true, ch_peak_pooled_pseudo

    script:
    """
    # peak calling for replicates
    macs2 callpeak -t ${bam_pr} -c ${bam_pr_ctrl} -f BAM -g ${params.genome_size} -n ${sample_id}_pooled -B -q ${params.macs_q}  2> ${sample_id}_pooled_macs2.log

    #Sort peak by -log10(p-value)
    sort -k8,8nr ${sample_id}_pooled_peaks.narrowPeak > ${sample_id}_pooled_sort_peaks.narrowPeak
    """
}

/*
 * 5. Peak calling for self replicates
 */
process PEAKCALLING_SELF_PSEUDOREPS {

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(bam_rep1_pr1), path(bam_ctrl_rep1_pr1), path(bam_rep1_pr2), path(bam_ctrl_rep1_pr2), path(bam_rep2_pr1), path(bam_ctrl_rep2_pr1), path(bam_rep2_pr2), path(bam_ctrl_rep2_pr2)  from ch_bam_self_pseudoreps

    output:
    tuple val(sample_id), path("${sample_id}_rep1_self_pseudorep1_sort_peaks.narrowPeak"), path("${sample_id}_rep1_self_pseudorep2_sort_peaks.narrowPeak"), path("${sample_id}_rep2_self_pseudorep1_sort_peaks.narrowPeak"), path("${sample_id}_rep2_self_pseudorep2_sort_peaks.narrowPeak")  into ch_peak_self_pseudoreps

    script:
    """
    ## Rep1
    # peak calling for replicates
    macs2 callpeak -t ${bam_rep1_pr1} -c ${bam_ctrl_rep1_pr1} -f BAM -g ${params.genome_size} -n ${sample_id}_rep1_self_pseudorep1 -B -q ${params.macs_q}  2> ${sample_id}_rep1_self_pseudorep1_macs2.log
    macs2 callpeak -t ${bam_rep1_pr2} -c ${bam_ctrl_rep1_pr2} -f BAM -g ${params.genome_size} -n ${sample_id}_rep1_self_pseudorep2 -B -q ${params.macs_q}  2> ${sample_id}_rep1_self_pseudorep2_macs2.log

    #Sort peak by -log10(p-value)
    sort -k8,8nr ${sample_id}_rep1_self_pseudorep1_peaks.narrowPeak > ${sample_id}_rep1_self_pseudorep1_sort_peaks.narrowPeak
    sort -k8,8nr ${sample_id}_rep1_self_pseudorep2_peaks.narrowPeak > ${sample_id}_rep1_self_pseudorep2_sort_peaks.narrowPeak

    ## Rep2
    # peak calling for replicates
    macs2 callpeak -t ${bam_rep2_pr1} -c ${bam_ctrl_rep2_pr1} -f BAM -g ${params.genome_size} -n ${sample_id}_rep2_self_pseudorep1 -B -q ${params.macs_q}  2> ${sample_id}_rep2_self_pseudorep1_macs2.log
    macs2 callpeak -t ${bam_rep2_pr2} -c ${bam_ctrl_rep2_pr2} -f BAM -g ${params.genome_size} -n ${sample_id}_rep2_self_pseudorep2 -B -q ${params.macs_q}  2> ${sample_id}_rep2_self_pseudorep2_macs2.log

    #Sort peak by -log10(p-value)
    sort -k8,8nr ${sample_id}_rep2_self_pseudorep1_peaks.narrowPeak > ${sample_id}_rep2_self_pseudorep1_sort_peaks.narrowPeak
    sort -k8,8nr ${sample_id}_rep2_self_pseudorep2_peaks.narrowPeak > ${sample_id}_rep2_self_pseudorep2_sort_peaks.narrowPeak
    """

}

/*
 * 6. Peak calling for pooled pseudo replicates
 */
process peakCALLING_POOLED_PSEUDOREPS {

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(bam_ppr1), path(bam_ppr1_ctrl), path(bam_ppr2),path(bam_ppr2_ctrl)  from ch_bam_pooled_pseudoreps

    output:
    tuple val(sample_id), path("${sample_id}_pooled_pseudorep1_sort_peaks.narrowPeak"), path("${sample_id}_pooled_pseudorep2_sort_peaks.narrowPeak")  into ch_peak_pooled_pseudoreps


    script:
    """
    # peak calling for replicates
    macs2 callpeak -t ${bam_ppr1} -c ${bam_ppr1_ctrl} -f BAM -g ${params.genome_size} -n ${sample_id}_pooled_pseudorep1 -B -q ${params.macs_q}  2> ${sample_id}_pooled_pseudorep1_macs2.log
    macs2 callpeak -t ${bam_ppr2} -c ${bam_ppr2_ctrl} -f BAM -g ${params.genome_size} -n ${sample_id}_pooled_pseudorep2 -B -q ${params.macs_q}  2> ${sample_id}_pooled_pseudorep2_macs2.log

    #Sort peak by -log10(p-value)
    sort -k8,8nr ${sample_id}_pooled_pseudorep1_peaks.narrowPeak > ${sample_id}_pooled_pseudorep1_sort_peaks.narrowPeak
    sort -k8,8nr ${sample_id}_pooled_pseudorep2_peaks.narrowPeak > ${sample_id}_pooled_pseudorep2_sort_peaks.narrowPeak
    """
}

/*
 * 7. IDR of true replicates
 */
process IDR_TRUE_REPS {
    publishDir "${params.outdir}/${sample_id}/idr/true_reps", mode: 'copy', pattern: '*.narrowPeak.gz'

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(peaks_rep1), path(peaks_rep2) from ch_peak_reps_true
    set val(sample_id), path(peaks_pooled)  from ch_peak_pooled_true

    output:
    tuple val(sample_id), path("${sample_id}_true_reps.IDR0.05.narrowPeak.gz")  into ch_idr_truereps


    script:
    """
    idr --samples ${peaks_rep1} ${peaks_rep1} \
    --peak-list ${peaks_pooled} \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file ${sample_id}_true_reps_idr \
    --plot \
    --soft-idr-threshold ${params.idr_threshold} \
    --use-best-multisummit-IDR \
    --log-output-file ${sample_id}_true_reps.idr.log

    idr_treshold_transformed=\$(awk -v p=${params.idr_threshold} 'BEGIN{print -log(p)/log(10)}')

    awk 'BEGIN{OFS="\t"} \$12>='"\${idr_treshold_transformed}"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${sample_id}_true_reps_idr | sort | uniq | sort -k7n,7n | gzip -nc > ${sample_id}_true_reps.IDR0.05.narrowPeak.gz
    """
}

/*
 * 8. IDR of self pseudo replicates
 */
process IDR_SELF_PSEUDOREPS {
    publishDir "${params.outdir}/${sample_id}/idr/self_pseudoreps", mode: 'copy', pattern: '*.narrowPeak.gz'

    when:
    !params.skip_idr

    input:
    set val(sample_id),  path(peaks_rep1_spr1), path(peaks_rep1_spr1),path(peaks_rep2_spr1), path(peaks_rep2_spr2) from ch_peak_self_pseudoreps
    set val(sample_id), path(peaks_rep1), path(peaks_rep2) from ch_peak_reps_self_pseudo

    output:
    tuple val(sample_id), path("${sample_id}_rep1_self_pseudoreps.IDR0.05.narrowPeak.gz"), path("${sample_id}_rep2_self_pseudoreps.IDR0.05.narrowPeak.gz")  into ch_idr_self_pseudoreps


    script:
    """
    ## Rep1
    idr --samples ${peaks_rep1_spr1} ${peaks_rep1_spr1} \
    --peak-list ${peaks_rep1} \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file ${sample_id}_rep1_self_pseudoreps_idr \
    --plot \
    --soft-idr-threshold ${params.idr_threshold} \
    --use-best-multisummit-IDR \
    --log-output-file ${sample_id}_rep1_self_pseudoreps.idr.log

    idr_treshold_transformed=\$(awk -v p=${params.idr_threshold} 'BEGIN{print -log(p)/log(10)}')

    awk 'BEGIN{OFS="\t"} \$12>='"\${idr_treshold_transformed}"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${sample_id}_rep1_self_pseudoreps_idr | sort | uniq | sort -k7n,7n | gzip -nc > ${sample_id}_rep1_self_pseudoreps.IDR0.05.narrowPeak.gz


    ## Rep1
    idr --samples ${peaks_rep2_spr1} ${peaks_rep2_spr1} \
    --peak-list ${peaks_rep2} \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file ${sample_id}_rep2_self_pseudoreps_idr \
    --plot \
    --soft-idr-threshold ${params.idr_threshold} \
    --use-best-multisummit-IDR \
    --log-output-file ${sample_id}_rep2_self_pseudoreps.idr.log

    awk 'BEGIN{OFS="\t"} \$12>='"\${idr_treshold_transformed}"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${sample_id}_rep2_self_pseudoreps_idr | sort | uniq | sort -k7n,7n | gzip -nc > ${sample_id}_rep2_self_pseudoreps.IDR0.05.narrowPeak.gz
    """
}

/*
 * 9. IDR of pooled pseudo replicates
 */
process IDR_POOLED_PSEUDOREPS {
    publishDir "${params.outdir}/${sample_id}/idr/pooled_pseudoreps", mode: 'copy', pattern: '*.narrowPeak.gz'

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(peaks_ppr1), path(peaks_ppr2) from ch_peak_pooled_pseudoreps
    set val(sample_id), path(peaks_pooled)  from ch_peak_pooled_pseudo

    output:
    tuple val(sample_id), path("${sample_id}_pooled_pseudoreps.IDR0.05.narrowPeak.gz")  into ch_idr_pooled_pseudoreps


    script:
    """
    idr --samples ${peaks_ppr1} ${peaks_ppr2} \
    --peak-list ${peaks_pooled} \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file ${sample_id}_pooled_pseudoreps_idr \
    --plot \
    --soft-idr-threshold ${params.idr_threshold} \
    --use-best-multisummit-IDR \
    --log-output-file ${sample_id}_pooled_pseudoreps.idr.log

    idr_treshold_transformed=\$(awk -v p=${params.idr_threshold} 'BEGIN{print -log(p)/log(10)}')

    awk 'BEGIN{OFS="\t"} \$12>='"\${idr_treshold_transformed}"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${sample_id}_pooled_pseudoreps_idr | sort | uniq | sort -k7n,7n | gzip -nc > ${sample_id}_pooled_pseudoreps.IDR0.05.narrowPeak.gz
    """
}


/*
 * 10. Optimal set & QC
 */
process QC {
    publishDir "${params.outdir}/${sample_id}/idr/optimal_set", mode: 'copy', pattern: '*.narrowPeak.gz'
    publishDir "${params.outdir}/${sample_id}/idr/QC", mode: 'copy', pattern: '*.txt'

    when:
    !params.skip_idr

    input:
    set val(sample_id), path(idr_tr) from ch_idr_truereps
    set val(sample_id), path(idr_ppr) from ch_idr_pooled_pseudoreps
    set val(sample_id), path(idr_spr1),path(idr_spr2) from ch_idr_self_pseudoreps


    output:
    tuple val(sample_id), path("${sample_id}_optimal_set.IDR0.05.narrowPeak.gz")  into ch_optimal_set
    path("${sample_id}_idr_QC.txt")  into ch_qc


    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Read peak files
    idr_tr = pd.read_csv('${idr_tr}', compression='gzip', header=None,sep='\t')
    idr_ppr = pd.read_csv('${idr_ppr}', compression='gzip', header=None,sep='\t')
    idr_spr1 = pd.read_csv('${idr_spr1}', compression='gzip', header=None,sep='\t')
    idr_spr2 = pd.read_csv('${idr_spr2}', compression='gzip', header=None,sep='\t')

    # Number of peaks
    Nt = idr_tr.shape[0]
    Np = idr_ppr.shape[0]
    N1 = idr_spr1.shape[0]
    N2 = idr_spr2.shape[0]

    # QC
    Rescue_ratio = max(Np,Nt) / min(Np,Nt)
    Self_consistency_ratio = max(N1,N2) / min(N1,N2)
    qc = pd.DataFrame([['Nt', round(Nt)], ['Np', Np], ['N1', N1], ['N2', N2],['Rescue Ratio', Rescue_ratio], ['Self-consistency Ratio', Self_consistency_ratio]] )
    qc.to_csv("${sample_id}_idr_QC.txt", index=False, sep='\t',header=False )

    if Nt > Np:
        Optimal_set = idr_tr
    else:
        Optimal_set = idr_ppr

    Optimal_set.to_csv("${sample_id}_optimal_set.IDR0.05.narrowPeak.gz", index=False, sep='\t', header=False, compression='gzip')
    """
}



/*
 * 11. Generate bigwigs
 */
process GENERATE_BIGWIGS {
    publishDir "${params.outdir}/${sample_id}/bigwigs", mode: 'copy', pattern: '*.bw'

    when:
    !params.skip_bigwig

    input:
    set val(sample_id), path(bam1), path(bam1_ctrl), path(bam2), path(bam2_ctrl) from ch_samples_split_bigwig
    set val(sample_id), path(bam_pr), path(bam_pr_ctrl)  from ch_bam_pooled_reps_for_bigwig

    output:
    tuple val(sample_id), path("${sample_id}_rep1.bw"), path("${sample_id}_Input_rep1.bw"), path("${sample_id}_rep2.bw"), path("${sample_id}_Input_rep2.bw")  into ch_bigwig_reps
    tuple val(sample_id), path("${sample_id}_pooled.bw"), path("${sample_id}_Input_pooled.bw")  into ch_bigwig_pooled

    script:
    """
    # Generate bigwigs for individual reps
    samtools index ${bam1}
    bamCoverage -b ${bam1} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample_id}_rep1.bw
    samtools index ${bam2}
    bamCoverage -b ${bam2} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample_id}_rep2.bw

    # Generate bigwigs for individual inputs
    samtools index ${bam1_ctrl}
    bamCoverage -b ${bam1_ctrl} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample_id}_Input_rep1.bw
    samtools index ${bam2_ctrl}
    bamCoverage -b ${bam2_ctrl} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample_id}_Input_rep2.bw

    # Generate bigwigs for pooled reps/inputs
    samtools index ${bam_pr}
    bamCoverage -b ${bam_pr} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample_id}_pooled.bw
    samtools index ${bam_pr_ctrl}
    bamCoverage -b ${bam_pr_ctrl} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample_id}_Input_pooled.bw
    """
}
