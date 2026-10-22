'''
MODULE: Downstream methylation analysis
DESCRIPTION:
    1. Unsupervised clustering (PCA + k-means / hclust / PAM, PERMANOVA validation)
    2. Differential methylation per chromosome with DSS (DMLs and DMRs)
    3. Genome-wide merge of DSS results
    4. DMR filtering, annotation (annotatr) and GO/KEGG enrichment (clusterProfiler)
REQUIRES (defined in the main Snakefile):
    CALL_DIR, METADATA, OUT_DOWN, SAMPLES, K_CLUST, VAR, TEST_GROUP, CHROMOSOMES
'''

# =============================================================================
# PIPELINE RULES
# =============================================================================
rule clustering_analysis:
    input:
        matrix = f"{CALL_DIR}/methyl_matrix_imputed.tsv",
        metadata = METADATA
    params:
        k = K_CLUST,
        batch_var = config['downstream'].get('batch_var') or "",
        samples_col = config['downstream']['samples_col'],
        outdir = f"{OUT_DOWN}/01_Clustering"
    output:
        expand(f"{OUT_DOWN}/01_Clustering/k_{K_CLUST}/{{f}}",
               f=["Scree_plot.svg", "Silhouette_k_selection.svg",
                  "Silhouette_kmeans.svg", "Silhouette_hclust.svg", "Silhouette_PAM.svg",
                  "adonis.txt", "pairwiseAdonis.tsv",
                  "PCA_2D_plot.svg", "PCA_3D_plot.html"])
    log:
        f"{OUT_DOWN}/01_Clustering/k_{K_CLUST}/Clustering_analysis.log"
    resources:
        mem_mb = config['downstream'].get('mem_cluster', 34000)
    envmodules:
        *config['modules']['r']
    script:
        "../scripts/Clustering_analysis.R"

rule differential_meth:
    input:
        cgmaps = expand(f"{CALL_DIR}/{{sample}}/{{sample}}_meth.CGmap.gz", sample=SAMPLES),
        metadata = METADATA
    params:
        samples = SAMPLES,
        samples_col = config['downstream']['samples_col'],
        indir = CALL_DIR,
        outdir = f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}",
        var = VAR,
        test_group = TEST_GROUP,
        min_cov = config['downstream'].get('dss_min_cov', 5)
    output:
        bsseq = f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}/BSseq_object_{{chrom}}.rds",
        dmls = f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}/DML_results_{VAR}_{{chrom}}.csv",
        dmrs = f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}/DMR_results_{VAR}_{{chrom}}.csv"
    log:
        f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}/differential_meth_{{chrom}}.log"
    benchmark:
        f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}/differential_meth_{{chrom}}.benchmark"
    wildcard_constraints:
        chrom = "|".join(re.escape(c) for c in CHROMOSOMES)
    resources:
        mem_mb = config['downstream'].get('mem_dss', 240000)
    envmodules:
        *config['modules']['r']
    script:
        "../scripts/Differential_DSS.R"

rule merge_differential_results:
    input:
        dmls = expand(f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}/DML_results_{VAR}_{{chrom}}.csv", chrom=CHROMOSOMES),
        dmrs = expand(f"{OUT_DOWN}/02_Diff_Methyl/{{chrom}}/DMR_results_{VAR}_{{chrom}}.csv", chrom=CHROMOSOMES)
    output:
        dmls_merged = f"{OUT_DOWN}/02_Diff_Methyl/DML_results_{VAR}_all.csv",
        dmrs_merged = f"{OUT_DOWN}/02_Diff_Methyl/DMR_results_{VAR}_all.csv"
    params:
        chrom_order = CHROMOSOMES
    run:
        import pandas as pd

        def merge(files, pos_col, out):
            # Skip chromosomes without results (header-only or empty files)
            dfs = []
            for f in files:
                try:
                    df = pd.read_csv(f)
                except pd.errors.EmptyDataError:
                    continue
                if not df.empty:
                    dfs.append(df)
            if not dfs:
                pd.DataFrame(columns=["chr", pos_col]).to_csv(out, index=False)
                return
            merged = pd.concat(dfs, ignore_index=True)
            # Natural chromosome order (chr1, chr2, ..., chrX, chrY)
            merged["chr"] = pd.Categorical(merged["chr"], categories=params.chrom_order, ordered=True)
            merged.sort_values(["chr", pos_col]).to_csv(out, index=False)

        merge(input.dmls, "pos", output.dmls_merged)    # DSS callDML: chr, pos, ...
        merge(input.dmrs, "start", output.dmrs_merged)  # DSS callDMR: chr, start, end, ...

rule dmr_annotation:
    input:
        dmrs = f"{OUT_DOWN}/02_Diff_Methyl/DMR_results_{VAR}_all.csv"
    params:
        gen_ver = config['downstream']['gen_ver'],
        test = TEST_GROUP,
        outdir = f"{OUT_DOWN}/03_DMR_annotation"
    output:
        dmrs_ann = f"{OUT_DOWN}/03_DMR_annotation/DMR_results_{VAR}_all_ann.csv"
    log:
        f"{OUT_DOWN}/03_DMR_annotation/DMR_annotation.log"
    resources:
        mem_mb = config['downstream'].get('mem_annot', 34000)
    envmodules:
        *config['modules']['r']
    script:
        "../scripts/DMR_annotation.R"
