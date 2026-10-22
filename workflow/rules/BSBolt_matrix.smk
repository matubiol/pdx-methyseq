'''
MODULE: BSBolt methylation matrix aggregation and imputation
DESCRIPTION:
    Aggregates per-sample CGmap files into a CpG x sample methylation matrix
    (proportions in [0, 1]) and imputes missing values with BSBolt's kNN imputation.
REQUIRES (defined in the main Snakefile):
    CALL_DIR, SAMPLES, APP['bsbolt']
'''

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule bsbolt_matrix:
    input:
        cgmaps = expand(f"{CALL_DIR}/{{sample}}/{{sample}}_meth.CGmap.gz", sample=SAMPLES)
    params:
        min_cov = config['bsbolt_matrix']['min_cov'],
        min_sample = config['bsbolt_matrix']['min_sample'],
        # Explicit sample labels so matrix columns are exactly the sample IDs
        labels = ",".join(SAMPLES)
    output:
        matrix = f"{CALL_DIR}/methyl_matrix.tsv",
        tmp_list = temp(f"{CALL_DIR}/cgmaps_list.txt")
    log:
        f"{CALL_DIR}/methyl_matrix.log"
    container:
        APP['bsbolt']
    threads:
        config['bsbolt_matrix']['thread_agg']
    resources:
        mem_mb = config['bsbolt_matrix']['mem_agg']
    shell:
        """
        printf '%s\\n' {input.cgmaps} > {output.tmp_list}

        bsbolt AggregateMatrix \
            -F {output.tmp_list} \
            -S {params.labels} \
            -CG \
            -min-coverage {params.min_cov} \
            -min-sample {params.min_sample} \
            -t {threads} \
            -O {output.matrix} \
            -verbose &> {log}
        """

rule impute_matrix:
    input:
        matrix = f"{CALL_DIR}/methyl_matrix.tsv"
    params:
        window_size = config['bsbolt_matrix']['window_size'],
        kneighbors = config['bsbolt_matrix']['kneighbors']
    output:
        imputed = f"{CALL_DIR}/methyl_matrix_imputed.tsv"
    log:
        f"{CALL_DIR}/methyl_matrix_imputed.log"
    container:
        APP['bsbolt']
    threads:
        config['bsbolt_matrix']['thread_impute']
    resources:
        mem_mb = config['bsbolt_matrix']['mem_impute']
    shell:
        """
        bsbolt Impute \
            -M {input.matrix} \
            -O {output.imputed} \
            -t {threads} \
            -W {params.window_size} \
            -k {params.kneighbors} \
            -verbose &> {log}
        """
