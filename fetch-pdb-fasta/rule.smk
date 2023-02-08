rule fetch_pdb_fasta:
    output:
        # Uncomment according to your use case.

        # 1. Fetch FASTA sequences per PDB entry
        # e.g., 4HHB.fasta
        # Uncomment below:
        # 'pdb_fasta_dir/{pdb_id}.fasta'
        
        # 2. Fetch FASTA sequence per polymer entity.
        # (identified by <pdb_id>_<entity_id>) e.g., 4HHB_1.fasta
        # Uncomment below:
        # 'pdb_fasta_dir/{pdb_id}_{entity_id}.fasta'

        # 3. Fetch FASTA sequence per polymer entity instance (chain)
        # (identified by <pdb_id>_<asym_id>) e.g., 4HHB_A.fasta
        # Uncomment below:
        'pdb_fasta_dir/{pdb_id}_{chain_id}.fasta'
    threads: 1
    resources: network=1
    shell:
        # Uncomment according to your use case.

        # 1. Fetch FASTA sequences per PDB entry. Uncomment below:
        # 'wget https://www.rcsb.org/fasta/entry/{pdb_id}/download -O {output}'

        # 2. Fetch FASTA sequence per polymer entity
        # (identified by <pdb_id>_<entity_id>). Uncomment below:
        # 'wget https://www.rcsb.org/fasta/entity/{pdb_id}_{entity_id}/download -O {output}'

        # 3. Fetch FASTA sequence per polymer entity instance (chain)
        # (identified by <pdb_id>_<asym_id>). Uncomment below:
        'wget https://www.rcsb.org/fasta/chain/{pdb_id}.{chain_id}/download -O {output}'
