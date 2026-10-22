'''
MODULE: BSBolt alignment, BAM processing, QC and methylation calling
DESCRIPTION:
    Aligns bisulfite/EM-seq reads with BSBolt, fixes mates, sorts and marks
    duplicates with samtools, runs samtools + Qualimap QC (aggregated with MultiQC)
    and calls CpG methylation with BSBolt CallMethylation.
    For xenograft samples, run Xengsort first and use the graft reads as input.
REQUIRES (defined in the main Snakefile):
    IN_DIR, R1_BS, R2_BS, OUT_BS, INDEX_BS, SAMPLES, APP['bsbolt'], APP['qualimap']
'''

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule bsbolt_align:
    input:
        f1 = lambda wc: f"{IN_DIR}/{R1_BS.format(sample=wc.sample)}",
        f2 = lambda wc: f"{IN_DIR}/{R2_BS.format(sample=wc.sample)}",
        index_dir = INDEX_BS
    params:
        prefix = f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}"
    output:
        bam = temp(f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}.bam")
    log:
        f"{OUT_BS}/01_Alignment/{{sample}}.log"
    benchmark:
        f"{OUT_BS}/01_Alignment/{{sample}}.benchmark"
    threads:
        config['bsbolt_call']['thread_align']
    resources:
        mem_mb = config['bsbolt_call']['mem_align']
    container:
        APP['bsbolt']
    shell:
        """
        bsbolt Align \
            -F1 {input.f1} \
            -F2 {input.f2} \
            -DB {input.index_dir} \
            -O {params.prefix} \
            -t {threads} \
            -OT 4 &> {log}
        """

# Fix mate information, sort by coordinate, mark duplicates and index
rule process_bam:
    input:
        bam = f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}.bam"
    params:
        tmp_prefix = f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}.sort_tmp"
    output:
        bam = f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}_dedup.bam",
        bai = f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}_dedup.bam.bai"
    log:
        f"{OUT_BS}/01_Alignment/{{sample}}_process_bam.log"
    threads:
        config['bsbolt_call']['thread_qc']
    resources:
        mem_mb = config['bsbolt_call']['mem_qc']
    envmodules:
        config['modules']['samtools']
    shell:
        """
        set -o pipefail
        (
            samtools fixmate -p -m {input.bam} - \
            | samtools sort -@ {threads} -T {params.tmp_prefix} - \
            | samtools markdup -@ {threads} - {output.bam}
            samtools index -@ {threads} {output.bam} {output.bai}
        ) &> {log}
        """

rule samtools_qc:
    input:
        f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}_dedup.bam"
    output:
        f"{OUT_BS}/02_Align_QC/{{sample}}/Samtools_stats.txt"
    threads:
        config['bsbolt_call']['thread_qc']
    resources:
        mem_mb = config['bsbolt_call']['mem_qc']
    envmodules:
        config['modules']['samtools']
    shell:
        """
        samtools quickcheck {input}
        samtools idxstats -@ {threads} {input} > {output}
        samtools flagstat -@ {threads} {input} >> {output}
        samtools stats -@ {threads} {input} >> {output}
        """

rule qualimap_qc:
    input:
        f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}_dedup.bam"
    output:
        directory(f"{OUT_BS}/02_Align_QC/{{sample}}/Qualimap")
    log:
        f"{OUT_BS}/02_Align_QC/{{sample}}/Qualimap.log"
    container:
        APP['qualimap']
    threads:
        config['bsbolt_call']['thread_qc']
    resources:
        mem_mb = config['bsbolt_call']['mem_qc']
    shell:
        """
        qualimap bamqc \
            -bam {input} \
            -outdir {output} \
            -nt {threads} \
            --java-mem-size={resources.mem_mb}M &> {log}
        """

rule multiqc:
    input:
        samtools = expand(f"{OUT_BS}/02_Align_QC/{{sample}}/Samtools_stats.txt", sample=SAMPLES),
        qualimap = expand(f"{OUT_BS}/02_Align_QC/{{sample}}/Qualimap", sample=SAMPLES)
    params:
        extra = "--verbose"
    output:
        report = f"{OUT_BS}/02_Align_QC/multiqc_report.html"
    log:
        f"{OUT_BS}/02_Align_QC/multiqc.log"
    wrapper:
        "v9.15.0/bio/multiqc"

rule bsbolt_call:
    input:
        bam = f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}_dedup.bam",
        bai = f"{OUT_BS}/01_Alignment/{{sample}}/{{sample}}_dedup.bam.bai"
    params:
        index_dir = INDEX_BS,
        out_prefix = f"{OUT_BS}/03_CallMethylation/{{sample}}/{{sample}}_meth"
    output:
        cgmap = f"{OUT_BS}/03_CallMethylation/{{sample}}/{{sample}}_meth.CGmap.gz"
    log:
        f"{OUT_BS}/03_CallMethylation/{{sample}}/{{sample}}_meth.log"
    benchmark:
        f"{OUT_BS}/03_CallMethylation/{{sample}}/{{sample}}_meth.benchmark"
    threads:
        config['bsbolt_call']['thread_call']
    resources:
        mem_mb = config['bsbolt_call']['mem_call']
    container:
        APP['bsbolt']
    shell:
        """
        bsbolt CallMethylation \
            -I {input.bam} \
            -DB {params.index_dir} \
            -O {params.out_prefix} \
            -CG \
            -t {threads} \
            -verbose &> {log}
        """
