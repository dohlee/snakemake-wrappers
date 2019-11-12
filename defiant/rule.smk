GROUP2FILES = {
    # SPECIFY YOUR METHYLATION FILES HERE.
}

def defiant_input(wildcards):
    return {
        'a': GROUP2FILES[wildcards.a],
        'b': GROUP2FILES[wildcards.b],
    }

rule defiant:
    input:
        unpack(defiant_input)
    output:
        '{a}_vs{b}.defiant.tsv'  
    params:
        extra = '',
        # Specify annotation file.
        a = False,
        # Output DMRs in bed file. This option does not take an argument.
        b = False,
        # Minimum coverage, e.g. "-c 10". This option accepts positive integers
        # and can be parallelized to test multiple options.
        c = False,
        # Minimum CpN/CpG/CH/CHH in a DMR, e.g. "-CpN 10". This option accepts
        # positive integers and can be parallelized. "CpN" is case insensitive.
        CpN = False,
        # Minimum differential nucleotide count in a DMR, e.g. "-d 3".
        # This option can be parallelized.
        d = False,
        # Turn on debugging mode. This slows down the execution significantly,
        # but can help diagnosse problem if they arise. This option does not
        # accept any arguments.
        debug = False,
        # Maximum non-default options in a parallel run, e.g. "-D 4"
        D = False,
        # Print statistics for every CpN. This option does not take an argument.
        # This slows Defiant down significantly. This file will start with
        # "nucleotide" and is used as input for the "regions_of_interest" program roi.c.
        E = False,
        # Make EPS figures for each DMR. Warning: requires R installation.
        # This option does not take an argument, and will slow defiant's execution.
        f = False,
        # Calculate FDR-adjusted q-value for each CpN. 'FDR' is case insensitive. 
        # This option can take case-insensitive arguments 'fdr' or 'bh' for Benjamini-
        # Hochberg method, 'Bonferroni', 'Hochberg', 'Hommel', 'Holm', or 'BY' for 
        # Benjamini & Yekutieli. If no argument is given, 'Holm' is assumed. 
        # This function is a translation of R's 'p.adjust'. I recommend against using
        # this as for genome-scale CpG measurements, almost everything will be
        # q = 1 and no DMRs will be obtained in any case. This option will
        # substantially increase RAM use and slow execution. 'Hommel' is so slow
        # I strongly recommend against it.
        fdr = False,
        # Maximum allowed gap between CpN, e.g. "-G 1000".
        G = False,
        # Set output file(s) label, e.g. "-l new"
        l = False,
        # List CpG Nucleotides in the DMR in output file. This option does not take
        # an argument.
        N = False,
        # Maximum p-value, which is 0<=p<=1. THis option can be parallelized to test
        # multiple options. Default 0.05.
        p = 0.05,
        # Minimum Percent methylation difference (0<=P<=100). This option can be parallelized
        # to test multiple options. Default 10%.
        P = 10,
        # Promoter cutoff for gene assignment of intergenic DMRs (default 10,000 nucleotides).
        # This option accepts positive integers, e.g. "-q 15000".
        q = False,
        # Minimum nucleotide range, which accepts a non-negative integer.
        # Default range is 0 nucleotides.
        r = False,
        # Include "Random" chromosomes. This option does not accept an argument.
        R = False,
        # Maximum allowed consecutive similar CpN, default is 5 CpN. This accepts non-negative
        # integers, e.g. "-s 3".
        s = 5,
        # Allow some number of consecutive skips of low coverage, default is 0. 
        # This accepts positive integers, e.g. "-S 1".
        S = 0,
        # Include "Un" chromsomes (default is to ignore them). This option does not accept an
        # argument.
        U = False,
        # Print a p-value for each DMR. This option accepts the same arguments that
        # the '-FDR' option does.
        v = False,
    threads: 1
    log: 'logs/defiant/{a}_vs{b}.log'
    benchmark: 'benchmarks/defiant/{a}_vs_{b}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/defiant'
