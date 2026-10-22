'''
MODULE: BSBolt index building
DESCRIPTION:
    Builds the BSBolt bisulfite index from a reference FASTA. The index is written
    to `bsbolt_call.index_dir` so that the alignment step picks it up directly.
REQUIRES (defined in the main Snakefile):
    REF_DIR, INDEX_BS, APP['bsbolt']
'''

REF_FILE = f"{REF_DIR}/{config['bsbolt_index']['ref_genome']}"

rule bsbolt_index:
    input:
        db = REF_FILE
    output:
        index_dir = directory(INDEX_BS)
    log:
        f"{INDEX_BS}.log"
    resources:
        mem_mb = config['bsbolt_index'].get('mem_mb', 48000)
    container:
        APP['bsbolt']
    shell:
        """
        bsbolt Index \
            -G {input.db} \
            -DB {output.index_dir} &> {log}
        """
