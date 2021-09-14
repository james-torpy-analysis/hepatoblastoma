# Run command:
# snakemake --reason --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N ewfus.smk -wd '/share/ScratchGeneral/jamtor/projects/hepatoblastoma/logs' -b y -j y -V -P DSGClinicalGenomics' -j 23

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

### This script remaps bams to hg19 and uses SvABA and Manta to identify 
# breakpoints in genomic data ###


# define variables:
project_name = 'hepatoblastoma'
capture_id = 'CDHS-33412Z-324'

# define/create directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/'
results_dir = project_dir + 'results/'
script_dir = project_dir + 'scripts/'

fq_dir = 'raw_files/'
smcounter_dir = 'results/smcounter2/'
bam_dir = 'results/picard/'
VAF_dir = 'results/detection_and_VAF/'

R_dir = "/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/envs/snkenv/bin/"

#SAMPLES = list([
#  '324_072_DBV4V_ATCTCAGG-CTCTCTAT_L001'
#])

SAMPLES = list([
    '324_001_DB674_TAAGGCGA-CTCTCTAT_L001', '324_003_DB674_AGGCAGAA-CTCTCTAT_L001', 
    '324_021_1_D9HGF_CTCTCTAC-CTCTCTAT_L001', '324_022_D9HGF_CGAGGCTG-CTCTCTAT_L001', 
    '324_025_D9HGF_GCTCATGA-CTCTCTAT_L001', 
    '324_026_DCR52_ATCTCAGG-CTCTCTAT_L001', '324_027_D9HGF_TAAGGCGA-CTCTCTAT_L001', 
    '324_028_D9HGF_CGTACTAG-CTCTCTAT_L001', '324_029_D9YW9_AGGCAGAA-CTCTCTAT_L001', 
    '324_030_D9HGF_TCCTGAGC-CTCTCTAT_L001', '324_031_D9HGF_GGACTCCT-CTCTCTAT_L001', 
    '324_032_D9HGF_TAGGCATG-CTCTCTAT_L001',  
    '324_033_DCR52_CTCTCTAC-CTCTCTAT_L001', '324_034_DB62M_CGAGGCTG-CTCTCTAT_L001', 
    '324_035_D9HGF_AAGAGGCA-CTCTCTAT_L001', '324_036_D9HGF_GTAGAGGA-CTCTCTAT_L001', 
    '324_037_D9HGF_AGGCAGAA-CTCTCTAT_L001', '324_039_D9YW9_CGTACTAG-CTCTCTAT_L001', 
    '324_040_D9YW9_TAGGCATG-CTCTCTAT_L001', '324_041_D9YW9_CTCTCTAC-CTCTCTAT_L001', 
    '324_042_D9YW9_CGAGGCTG-CTCTCTAT_L001', '324_043_D9YW9_AAGAGGCA-CTCTCTAT_L001', 
    '324_044_D9YW9_GCTCATGA-CTCTCTAT_L001', 
    '324_045_DCR52_TAAGGCGA-CTCTCTAT_L001', '324_046_D9YW9_GTAGAGGA-CTCTCTAT_L001', 
    '324_047_D9YW9_ATCTCAGG-CTCTCTAT_L001', '324_048_D9YWF_CGTACTAG-CTCTCTAT_L001', 
    '324_049_D9YWF_TCCTGAGC-CTCTCTAT_L001', '324_050_D9YWF_GGACTCCT-CTCTCTAT_L001', 
    '324_051_D9YWF_TAGGCATG-CTCTCTAT_L001', '324_052_D9YWF_CTCTCTAC-CTCTCTAT_L001', 
    '324_053_DCR52_CGAGGCTG-CTCTCTAT_L001', 
    '324_054_D9YWF_AAGAGGCA-CTCTCTAT_L001', '324_055_D9YWF_GCTCATGA-CTCTCTAT_L001', 
    '324_056_D9YWF_TAAGGCGA-CTCTCTAT_L001', '324_057_DB674_CGTACTAG-CTCTCTAT_L001', 
    '324_058_DB674_TCCTGAGC-CTCTCTAT_L001', '324_059_DB674_GGACTCCT-CTCTCTAT_L001', 
    '324_060_DB674_TAGGCATG-CTCTCTAT_L001', '324_061_DB674_CTCTCTAC-CTCTCTAT_L001', 
    '324_062_DB674_CGAGGCTG-CTCTCTAT_L001', '324_063_DB674_AAGAGGCA-CTCTCTAT_L001', 
    '324_064_DB674_GTAGAGGA-CTCTCTAT_L001', '324_065_DB674_GCTCATGA-CTCTCTAT_L001', 
    '324_066_DB674_ATCTCAGG-CTCTCTAT_L001', '324_067_DB62M_TAAGGCGA-CTCTCTAT_L001', 
    '324_068_DB62M_CGTACTAG-CTCTCTAT_L001', '324_069_DB62M_AGGCAGAA-CTCTCTAT_L001', 
    '324_070_DB62M_TCCTGAGC-CTCTCTAT_L001', '324_071_DB62M_AAGAGGCA-CTCTCTAT_L001', 
    '324_072_DBV4V_ATCTCAGG-CTCTCTAT_L001', '409_006_DB62M_GGACTCCT-CTCTCTAT_L001', 
    '409_007_DB62M_TAGGCATG-CTCTCTAT_L001', '409_008_DB62M_GTAGAGGA-CTCTCTAT_L001', 
    '409_009_DB62M_GCTCATGA-CTCTCTAT_L001', '409_010_DB62M_ATCTCAGG-CTCTCTAT_L001'
])

rule all:
    input:
        expand(
            VAF_dir + '{sample}/Rdata/selected_VAF.rds',
            sample = SAMPLES
        )


######################################################################################################
### 1. Trim reads using smCounter2 ###
######################################################################################################

rule trim_reads:
    input:
        fq1 = fq_dir + '{sample}/{sample}_R1.fastq.gz',
        fq2 = fq_dir + '{sample}/{sample}_R2.fastq.gz'
    output:
        smcounter_dir + '{sample}/{sample}.align.bam'
    threads: 8
    shell:
        "mkdir -p logs/trim_reads/{wildcards.sample}/; " + 
        "cd logs/trim_reads/{wildcards.sample}/; " +
        " ../../../scripts/1.smcounter2.sh " +
        "{wildcards.sample} " +
        "{capture_id} 2> {wildcards.sample}.smcounter2.errors"


######################################################################################################
### 2. Deduplicate UMIs ###
######################################################################################################

rule dedup:
    input:
        smcounter_dir + '{sample}/{sample}.align.bam'
    output:
        bam_dir + '{sample}/{sample}.dedup.sorted.by.coord.bam'
    threads: 8
    shell:
        "mkdir -p logs/dedup/{wildcards.sample}/; " + 
        "cd logs/dedup/{wildcards.sample}/; " +
        " ../../../scripts/2.UMI_dedup.sh " +
        "{wildcards.sample} 2> {wildcards.sample}.dedup.errors"


######################################################################################################
### 3. Detection and VAFs ###
######################################################################################################

rule detect_and_VAF:
    input:
        bam_dir + '{sample}/{sample}.dedup.sorted.by.coord.bam'
    output:
        VAF_dir + '{sample}/Rdata/selected_VAF.rds'
    threads: 8
    shell:
        "mkdir -p logs/detect_and_VAF/{wildcards.sample}/; " + 
        "cd logs/detect_and_VAF/{wildcards.sample}/; " +
        "{R_dir}/R CMD BATCH  --no-save '--args" + 
        " {project_name}" + 
        " {wildcards.sample}" + 
        " 10" + # min total reads needed for breakpoint to be considered
        " 19" + # min overlap required for split reads vs breakpoints
        " 2" + # min supporting reads required for VAF
        "' ../../../scripts/3.detect_and_vaf.R"



