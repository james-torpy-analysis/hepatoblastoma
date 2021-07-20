# Run command:
# snakemake --reason --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N ewfus.smk -wd '/share/ScratchGeneral/jamtor/projects/hepatoblastoma/logs' -b y -j y -V -P DSGClinicalGenomics' -j 23

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

### This script remaps bams to hg19 and uses SvABA and Manta to identify 
# breakpoints in genomic data ###


# define variables:
project_name = 'hepatoblastoma'
chromosomes = 'chr3'
ROI = "41265552_41266717" # region covered by primers
supps_allowed = '2'
SV_type = 'deletion'
min_overlap = '19'
disc_read_window = '200'
same_SV_window = '10'

# define/create directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/'
results_dir = project_dir + 'results/'
genome_dir = project_dir + 'genome/'
script_dir = project_dir + 'scripts/'

conda_dir = '/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/'
env_dir = conda_dir + 'envs/snkenv/bin/'

fq_dir = 'raw_files/'
align_dir = 'results/BWA_and_picard/bams/'
svaba_dir = 'results/svaba/BWA_and_picard/'
int_dir = 'results/BWA_and_picard/int_bams/'
SV_dir = 'results/SVs/'
VAF_dir = 'results/VAF_calculation/'

R_dir = "/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/envs/snkenv/bin/"

#SAMPLES = list([
#    '324_003_DB674_AGGCAGAA-CTCTCTAT_L001'
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

#rule all:
#    input:
#        expand(
#            align_dir + '{sample}/{sample}.consensus.bam.bai',
#            sample = SAMPLES
#        )

#rule all:
#    input:
#        expand(
#            svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
#            sample=SAMPLES
#        )

#rule all:
#    input:
#        expand(
#            'logs/completed_jobs/{sample}_complete',
#            sample = SAMPLES
#        )

rule all:
    input:
        expand(
            SV_dir + '{sample}/detected_SVs.Rdata',
            sample = SAMPLES
        )

#rule all:
#    input:
#        expand(
#            VVAF_dir + '{sample}/non_specific_SV_supporting_reads.tsv',
#            sample = SAMPLES
#        )


######################################################################################################
### 1. UMI collapse and BWA ###
######################################################################################################

rule BWA_and_umi_collapse:
    input:
        fq1 = fq_dir + '{sample}/{sample}_R1.fastq.gz',
        fq2 = fq_dir + '{sample}/{sample}_R2.fastq.gz'
    output:
        bam = align_dir + '{sample}/{sample}.consensus.bam',
        bai = align_dir + '{sample}/{sample}.consensus.bam.bai',
    threads: 8
    shell:
        'mkdir -p logs/BWA_and_picard; ' +
        'cd logs/BWA_and_picard; ' + 
        'touch {wildcards.sample}; ' + 
        script_dir + '1.UMI_collapse.sh' +
            ' {project_name}' + 
            ' {wildcards.sample}' +
            ' 2>&1 {wildcards.sample}.alignment.log'


######################################################################################################
### 2. SvABA ###
######################################################################################################

rule svaba:
   input:
       bam = align_dir + '{sample}/{sample}.consensus.bam',
       bai = align_dir + '{sample}/{sample}.consensus.bam.bai'
   output:
       filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       unfilt = svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
   threads: 8
   shell:
       'mkdir -p logs/svaba; ' + 
        'cd logs/svaba; ' + 
        'mkdir -p ../../' + svaba_dir + '{wildcards.sample}/; ' +
        'svaba run -t ' + project_dir + '{input.bam} -G ' + 
            genome_dir + 'GRCh37.p13.genome.fa -a ../../' + 
            svaba_dir + 
            '{wildcards.sample}/{wildcards.sample}' + 
            ' -p 6 --override-reference-check' +
            ' --min-overlap 0.13' +
            ' 2> {wildcards.sample}.svaba.errors'

rule format_vcf:
    input:
        svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
    output:
        svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
    threads: 1
    shell:
        'scripts/fix_broken_svaba_vcf.sh {wildcards.sample} {input}'

rule remove_duds:
    input:
        svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
    output:
        svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf'
    threads: 1
    shell:
        'grep -v NODISC {input} > {output}'

rule vcf_index:
   input:
       filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       semifilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf'
   output:
        filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        semifilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf.idx'
   threads: 2
   shell:
        env_dir + 'igvtools index {input.filt}; ' + 
        env_dir + 'igvtools index {input.semifilt}'

    
######################################################################################################
### 3. Clean up ###
######################################################################################################

rule cleanup:
    input:
        filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf.idx'
    output:
        'logs/completed_jobs/{sample}_complete'
    threads: 1
    shell:
        'rm -fr ' + int_dir + '{wildcards.sample}; '
        'touch {output}'


######################################################################################################
### 4. Find SVs ###
######################################################################################################

rule find_specific_SVs:
    input:
        filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf.idx'
    output:
        SV_dir + '{sample}/detected_SVs.Rdata'
    threads: 8
    shell:
        "mkdir -p logs/find_SVs/{wildcards.sample}/; " + 
        "cd logs/find_SVs/{wildcards.sample}/; " +
        "{R_dir}/R CMD BATCH  --no-save '--args" + 
        " {project_name}" + 
        " {wildcards.sample}" + 
        " {chromosomes}" + 
        "' ../../../scripts/2.find_specific_SVs.R"


######################################################################################################
### 5. Calculate VAFs ###
######################################################################################################

rule calc_VAFs:
    input:
        SV_dir + '{sample}/detected_SVs.Rdata'
    output:
        VAF_dir + '{sample}/Rdata/VAFs.Rdata'
    threads: 8
    shell:
        "mkdir -p logs/VAF_calculation/{wildcards.sample}/; " + 
        "cd logs/VAF_calculation/{wildcards.sample}/; " +
        "{R_dir}/R CMD BATCH  --no-save '--args" + 
        " {project_name}" + 
        " {wildcards.sample}" + 
        "' ../../../scripts/3.vaf.R"


######################################################################################################
### 6. Find supporting reads ###
######################################################################################################

rule find_supp:
    input:
        VAF_dir + '{sample}/Rdata/VAFs.Rdata'
    output:
        VAF_dir + '{sample}/non_specific_SV_supporting_reads.tsv'
    threads: 8
    shell:
        "mkdir -p logs/find_supp/{wildcards.sample}/; " + 
        "cd logs/find_supp/{wildcards.sample}/; " +
        "{R_dir}/R CMD BATCH  --no-save '--args" + 
        " {project_name}" + 
        " {wildcards.sample}" + 
        " 19" + # min overlapp required for read 1 to be counted as supporting
        " 19" + # min overlapp required for read 2 to be counted as supporting
        "' ../../../scripts/4.find_supporting_reads.R"

