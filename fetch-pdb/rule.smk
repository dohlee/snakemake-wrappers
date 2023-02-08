rule fetch_pdb:
    output:
        # Required output.
        '{pdb_id}.pdb'
    threads: 1
    resources: network=1
    wrapper:
        'http://dohlee-bio.info:9193/fetch-pdb'
