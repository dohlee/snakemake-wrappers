rule fetch_pdb:
    output: 'pdb_dir/{pdb_id}.pdb'
    threads: 1
    resources: network=1
    shell: 'wget https://files.rcsb.org/download/{pdb_id}.pdb -O {output}'
