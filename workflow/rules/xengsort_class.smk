'''
MODULE: Xengsort human/mouse read classification
DESCRIPTION:
    Classifies paired-end reads into host / graft / both / neither / ambiguous.
    Only the graft (human) reads are passed on to BSBolt; the other categories
    are kept for QC and can be deleted once the MultiQC report has been reviewed.
REQUIRES (defined in the main Snakefile):
    RAW_DIR, R1_PATTERN, R2_PATTERN, OUT_XENG, INDEX_XENG, SAMPLES, APP['xengsort']
'''

XENG_CATEGORIES = ["host", "graft", "both", "neither", "ambiguous"]

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule xengsort_classify:
    input:
        r1 = lambda wc: f"{RAW_DIR}/{R1_PATTERN.format(sample=wc.sample)}",
        r2 = lambda wc: f"{RAW_DIR}/{R2_PATTERN.format(sample=wc.sample)}",
        index_info = f"{INDEX_XENG}.info",
        index_hash = f"{INDEX_XENG}.hash"
    params:
        index = INDEX_XENG,
        outprefix = f"{OUT_XENG}/{{sample}}/{{sample}}",
        chunksize = config['xeng_class']['chunksize'],
        prefetch = config['xeng_class']['prefetch']
    output:
        graft1 = f"{OUT_XENG}/{{sample}}/{{sample}}-graft.1.fq.gz",
        graft2 = f"{OUT_XENG}/{{sample}}/{{sample}}-graft.2.fq.gz",
        others = expand(f"{OUT_XENG}/{{{{sample}}}}/{{{{sample}}}}-{{cat}}.{{read}}.fq.gz",
                        cat=[c for c in XENG_CATEGORIES if c != "graft"], read=[1, 2])
    log:
        f"{OUT_XENG}/{{sample}}/{{sample}}.log"
    benchmark:
        f"{OUT_XENG}/{{sample}}/{{sample}}.benchmark"
    threads: config['xeng_class']['threads']
    resources:
        mem_mb = config['xeng_class']['mem_mb']
    container:
        APP['xengsort']
    shell:
        """
        xengsort -DD classify \
            --fastq {input.r1} \
            --pairs {input.r2} \
            --index {params.index} \
            --threads {threads} \
            --chunksize {params.chunksize} \
            --out {params.outprefix} \
            --prefetchlevel {params.prefetch} &> {log}
        """

# MultiQC config: strip the .log extension and make the xengsort module
# search only the first lines of each log (avoids scanning huge files).
rule xengsort_multiqc_config:
    output:
        config = temp(f"{OUT_XENG}/multiqc_config.yaml")
    run:
        with open(output.config, "w") as f:
            f.write('fn_clean_exts:\n    - ".log"\n')
            f.write('sp:\n    xengsort:\n        contents: "# Xengsort classify"\n        num_lines: 10\n')

rule xengsort_multiqc:
    input:
        expand(f"{OUT_XENG}/{{sample}}/{{sample}}.log", sample=SAMPLES),
        config = f"{OUT_XENG}/multiqc_config.yaml"
    output:
        f"{OUT_XENG}/multiqc_report.html"
    params:
        extra = "--verbose"
    log:
        f"{OUT_XENG}/multiqc.log"
    wrapper:
        "v9.15.0/bio/multiqc"
