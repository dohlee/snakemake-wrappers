rule multiqc:
    input:
        expand('FASTQC_OUTDIR/{run}_fastqc.zip', run=SOME_RUNS)
    output:
        'FASTQC_OUTDIR/multiqc_report.html',
    params:
        extra = '',
        # Overwrite any existing reports.
        # Default: True
        force = True,
        # Prepend directory to sample names.
        # Default: False
        dirs = False,
        # Prepend [INT] directories to sample names.
        # Negative number to take from start of path.
        # Default: False
        dirs_depth = False,
        # Do not clean the sample names (leave as full file name).
        # Default: False
        fullnames = False,
        # Report title. Printed as page header, used for filename
        # if not otherwise specified.
        # Default: False,
        title = False,
        # Custom comment, will be printed at the top of the report.
        # Default: False
        comment = False,
        # Create report in the specified output directory.
        # Default: False
        outdir = False,
        # Report template to use.
        # Options: [default|default_dev|geo|sections|simple]
        # Default: False
        template = False,
        # Use only modules which tagged with this keyword, eg. RNA.
        # Default: False
        tag = False,
        # Ignore analysis files (glob expression).
        # Default: False
        ignore = False, # Ignore symlink directories and files.
        # Default: False
        ignore_symlinks = False,
        # File containing alternative sample names.
        # Default: False
        sample_names = False,
        # Supply a file containing a list of file paths to be searched,
        # one per row.
        # Default: False
        file_list = False,
        # Do not use this module. Can specify multiple times.
        # Default: False
        exclude = False,
        # Use only this module. Can specify multiple times.
        # Default: False
        module = False,
        # Force the parsed data directory from being created.
        # Default: False
        data_dir = False,
        # Prevent the parsed data directory from being created.
        # Default: False
        no_data_dir = False,
        # Output parsed data in a different format.
        # Options: [tsv|json|yaml]
        # Default: tsv
        data_format = 'tsv',
        # Compress the data directory.
        # Default: False
        zip_data_dir = False,
        # Export plots as static images in addition to the report.
        # Default: False
        export = False,
        # Use only flat plots (static images)
        # Default: False
        flat = False,
        # Use only interactive plots (HighCharts Javascript).
        # Default: False
        interactive = False,
        # Use strict linting (validation) to help code development.
        # Default: False
        lint = False,
        # Creates PDF report with 'simple' template.
        # Requires Pandoc to be installed.
        # Default: False
        pdf = False,
        # Don't upload generated report to MegaQC, even if MegaQC options are found.
        # Default: False
        no_megaqc_upload = False,
        # Specific config file to load, after those in MultiQC dir / home dir / working dir.
        # Default: False
        config = False,
        # Specify MultiQC config YAML on the commandline.
        # Default: False
        cl_config = False,
        # Increase output verbosity.
        # Default: False
        verbose = False,
        # Only show log warnings.
        # Default: False
        quiet = False,
    threads: 1
    log: 'logs/multiqc.log'
    wrapper:
        'http://dohlee-bio.info:9193/multiqc'
