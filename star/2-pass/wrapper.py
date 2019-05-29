__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.input[parameter])
        if param:
            return option + ' ' + str(snakemake.input[parameter])
        else:
            return ''
    except AttributeError:
        return ''

def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        if str(snakemake.params[parameter]) == '':
            return ''
        if type(snakemake.params[parameter]) == bool:
            if snakemake.params[parameter]:
                return option
            else:
                return ''
        else:
            return option + ' ' + str(snakemake.params[parameter])
    except AttributeError:
        return '' 

# Extract required inputs.
reads = snakemake.input.reads
if isinstance(reads, str):
    reads = [reads]
assert len(reads) in [1, 2], "Currently star/2-pass wrapper only supports single-sample analysis."
star_index = snakemake.input.star_index

# Extract required outputs.
genome_alignment = snakemake.output.genome_alignment
transcriptome_alignment = snakemake.output.get('transcriptome_alignment', None)
output_prefix = genome_alignment[:-10]  # Strip tailing 'sorted.bam'.

# Extract parameters.
# Extract optional parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('run_rng_seed', '--runRNGseed'))
user_parameters.append(optionify_params('quant_mode', '--quantMode'))
user_parameters.append(optionify_params('genome_load', '--genomeLoad'))
user_parameters.append(optionify_params('sjdb_file_chr_start_end', '--sjdbFileChrStartEnd'))
user_parameters.append(optionify_params('sjdb_gtf_file', '--sjdbGTFfile'))
user_parameters.append(optionify_params('sjdb_gtf_chr_prefix', '--sjdbGTFchrPrefix'))
user_parameters.append(optionify_params('sjdb_overhang', '--sjdbOverhang'))
user_parameters.append(optionify_params('sjdb_gtf_tag_exon_parent_transcript', '--sjdbGTFtagExonParentTranscript'))
user_parameters.append(optionify_params('out_filter_type', '--outFilterType'))
user_parameters.append(optionify_params('out_filter_multimap_score_range', '--outFilterMultimapScoreRange'))
user_parameters.append(optionify_params('out_filter_multimap_n_max', '--outFilterMultimapNmax'))
user_parameters.append(optionify_params('out_filter_mismatch_n_max', '--outFilterMismatchNmax'))
user_parameters.append(optionify_params('out_filter_mismatch_n_over_l_max', '--outFilterMismatchNoverLmax'))
user_parameters.append(optionify_params('out_filter_mismatch_n_over_read_l_max', '--outFilterMismatchNoverReadLmax'))
user_parameters.append(optionify_params('out_filter_score_min', '--outFilterScoreMin'))
user_parameters.append(optionify_params('out_filter_score_min_over_l_read', '--outFilterScoreMinOverLread'))
user_parameters.append(optionify_params('out_filter_match_n_min', '--outFilterMatchNmin'))
user_parameters.append(optionify_params('out_filter_match_n_min_over_l_read', '--outFilterMatchNminOverLread'))
user_parameters.append(optionify_params('out_filter_intron_motifs', '--outFilterIntronMotifs'))
user_parameters.append(optionify_params('out_filter_intron_strands', '--outFilterIntronStrands'))
user_parameters.append(optionify_params('out_sj_filter_reads', '--outSJfilterReads'))
user_parameters.append(optionify_params('out_sj_filter_overhang_min', '--outSJfilterOverhangMin'))
user_parameters.append(optionify_params('out_sj_filter_count_unique_min', '--outSJfilterCountUniqueMin'))
user_parameters.append(optionify_params('out_sj_filter_count_total_min', '--outSJfilterCountTotalMin'))
user_parameters.append(optionify_params('out_sj_filter_dist_to_other_sj_min', '--outSJfilterDistToOtherSJmin'))
user_parameters.append(optionify_params('out_sj_filter_intron_max_vs_read_n', '--outSJfilterIntronMaxVsReadN'))
user_parameters.append(optionify_params('align_intron_min', '--alignIntronMin'))
user_parameters.append(optionify_params('align_intron_max', '--alignIntronMax'))
user_parameters.append(optionify_params('align_sj_overhang_min', '--alignSJoverhangMin'))

# If gzipped reads are given, but user did not specify readFilesCommand option,
# kindly add the option.
read_files_command = '--readFilesCommand cat'
if reads[0].endswith('.gz') and '--readFilesCommand' not in extra:
    read_files_command = '--readFilesCommand zcat'
user_parameters.append(read_files_command)

user_parameters = ' '.join([p for p in user_parameters if p != ''])

# rename {output_prefix}Aligned.sortedByCoord.out.bam into {output_prefix}.sorted.bam
rename_command = '&& mv {prefix}Aligned.sortedByCoord.out.bam {prefix}sorted.bam'.format(prefix=output_prefix)

transcriptome_rename_command = ''
if snakemake.params.quant_mode == 'TranscriptomeSAM':
    transcriptome_rename_command = '&& mv {prefix}Aligned.toTranscriptome.sortedByCoord.out.bam {prefix}.transcriptome.sorted.bam'.format(prefix=output_prefix)

# Execute shell command.
shell(
    "("
    "STAR "
    "--runMode alignReads "
    "--twopassMode Basic "
    "--runThreadN {snakemake.threads} "
    "--readFilesIn {reads} "
    "--genomeDir {star_index} "
    "--outFileNamePrefix {output_prefix} "
    "--outSAMtype BAM SortedByCoordinate "
    "{extra} "
    "{user_parameters} "
    "{rename_command} "
    "{transcriptome_rename_command} "
    ") "
    "{log}"
)
