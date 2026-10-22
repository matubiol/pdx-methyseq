'''
MODULE: Xengsort index building (RNA-seq / WES / WGBS)
DESCRIPTION:
    Builds a human (graft) / mouse (host) k-mer index for xenograft read sorting.
    - wgbs:   genomic DNA, bisulfite-aware index (--bisulfite)
    - wes:    genomic DNA
    - rnaseq: genomic DNA + cDNA
REQUIRES (defined in the main Snakefile):
    REF_DIR, DATA_TYPE, INDEX_XENG, K_XENG, APP['xengsort']
'''

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
BISULFITE = DATA_TYPE == "wgbs"
USE_CDNA = DATA_TYPE == "rnaseq"

DNA_HUMAN = str(REF_DIR / config['xeng_ref']['dna_human'])
DNA_MOUSE = str(REF_DIR / config['xeng_ref']['dna_mouse'])
HOST_REFS = [DNA_MOUSE]
GRAFT_REFS = [DNA_HUMAN]
if USE_CDNA:
    HOST_REFS.append(str(REF_DIR / config['xeng_ref']['cdna_mouse']))
    GRAFT_REFS.append(str(REF_DIR / config['xeng_ref']['cdna_human']))

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule xengsort_index:
    input:
        host = HOST_REFS,
        graft = GRAFT_REFS
    params:
        index = INDEX_XENG,
        size = config['xeng_index']['size'],
        fill = config['xeng_index']['fill'],
        weakthreads = config['xeng_index']['weakthreads'],
        ksize = K_XENG,
        bisulfite_flag = "--bisulfite" if BISULFITE else ""
    output:
        info = f"{INDEX_XENG}.info",
        hash = f"{INDEX_XENG}.hash"
    log:
        f"{INDEX_XENG}.log"
    benchmark:
        f"{INDEX_XENG}.benchmark"
    threads: config['xeng_index']['threads']
    resources:
        mem_mb = config['xeng_index']['mem_mb']
    container:
        APP['xengsort']
    shell:
        """
        xengsort -DD index \
            --index {params.index} \
            -H {input.host} \
            -G {input.graft} \
            -n {params.size} \
            -k {params.ksize} \
            --fill {params.fill} \
            {params.bisulfite_flag} \
            -W {params.weakthreads} &> {log}
        """
