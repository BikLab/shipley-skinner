############### Qiime Scripts used for 16S/18S rRNA metabarcoding datasets ###############
### Input/output files should be changed based on the gene being analyzed
### These scripts should be used on a QIIME2 artifact that contained the all demultiplexed files per gene
### Demultiplexing of MiSeq files were done using the Perl script "Demul_trim_prep_flipedmerge2.1.pl"

##### 00_Importing_Silva_Taxonomy_18S
#!/bin/bash
#SBATCH --job-name="Importing_Silva_Taxonomy_18S"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=24:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e Importing_Silva_Taxonomy_18S.err-%N
#SBATCH -o Importing_Silva_Taxonomy_18S.out-%N

#Path Variables
INPUTSEQ=/scratch/tjp99569/shsk/database/CustomDatabase_SILVA_and_SingleNem_GenBank_Sequences_Nov_06_2021.fna
OUTPUTSEQ=/scratch/tjp99569/shsk/database/CustomDatabase_SILVA_and_SingleNem_GenBank_Sequences_Nov_06_2021.qza
INPUTAX=/scratch/tjp99569/shsk/database/CustomDatabase_SILVA_and_SingleNem_GenBank_Taxonomy_Nov_06_2021.tsv
OUTPUTAX=/scratch/tjp99569/shsk/database/CustomDatabase_SILVA_and_SingleNem_GenBank_Taxonomy_Nov_06_2021.qza

#Activate QIIME2/2020.11
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.11

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

# Qiime import commands

# The sequences
qiime tools import \
 --type 'FeatureData[Sequence]' \
 --input-path $INPUTSEQ \
 --output-path $OUTPUTSEQ \

# The taxonomy strings
qiime tools import \
 --type 'FeatureData[Taxonomy]' \
 --input-path $INPUTAX \
 --output-path $OUTPUTAX \
 --input-format HeaderlessTSVTaxonomyFormat
 
##### 01-Import_Manifest_SHSK.sh
#!/bin/bash
#SBATCH --job-name="Import_18S_Fastq_SHSK"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e Import_18S_Fastq_SHSK.err-%N
#SBATCH -o Import_18S_Fastq_SHSK.out-%N

#Path Variables
MANIFEST=/scratch/tjp99569/shsk/qiime-files/18S-Soil-and-SingleNematodes-manifest-01-22-2021.txt
ARTIFACT=/home/tjp99569/shsk/dada2-results-18S/01-SHSK_18S_Samples_ArtifactFile_notrimmed.qza
INPUTQZA=/home/tjp99569/shsk/dada2-results-18S/01-SHSK_18S_Samples_ArtifactFile_notrimmed.qza
INPUTQZV=/home/tjp99569/shsk/dada2-results-18S/01-SHSK_18S_Samples_ArtifactFile_notrimmed.qzv

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Create Qiime Artifact
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path $ARTIFACT \
  --input-format PairedEndFastqManifestPhred33V2 \

#Visualize Quality Scores
qiime demux summarize \
  --i-data $INPUTQZA \
  --o-visualization $INPUTQZV

##### 02-Remove_Primers_SHSK.sh
#!/bin/bash
#SBATCH --job-name="Removing_Primers_18S_SHSK"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e Removing_Primers_18S_SHSK.err-%N
#SBATCH -o Removing_Primers_18S_SHSK.out-%N

#Path Variables
UNTRIMMED=/home/tjp99569/shsk/dada2-results-18S/01-SHSK_18S_Samples_ArtifactFile_notrimmed.qza
TRIMMED=/home/tjp99569/shsk/dada2-results-18S/01-SHSK_18S_Samples_ArtifactFile_trimmed.qza
TRIMMEDQVZ=/home/tjp99569/shsk/dada2-results-18S/01-SHSK_18S_Samples_ArtifactFile_trimmed.qzv

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#cutadapt parameters
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $UNTRIMMED \
  --p-cores 5 \
  --p-adapter-f GGCTTAATCTTTGAGACAAGCCGAATTACCATA \
  --p-adapter-r TCCAAGGAAGGCAGCAGGCTACTGGCTGACT \
  --p-error-rate 0.11 \
  --p-discard-untrimmed True \
  --o-trimmed-sequences $TRIMMED \

#Visualize Quality Scores
qiime demux summarize \
  --i-data $TRIMMED \
  --o-visualization $TRIMMEDQVZ

#### 03-AssignASVs_DADA2_SHSK.sh
#!/bin/bash
#SBATCH --job-name="AssignASVs_DADA2_18S_SHSK"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e AssignASVs_DADA2_18S_SHSK.err-%N
#SBATCH -o AssignASVs_DADA2_18S_SHSK.out-%N

#Path Variables
INPUT=/home/tjp99569/shsk/dada2-results-18S/01-SHSK_18S_Samples_ArtifactFile_notrimmed.qza
TABLE=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureTable_DADA2.qza
SEQS=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2.qza
STATS=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_Denoise_Stats_DADA2.qza

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc


#Dada2 parameters
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $INPUT \
  --p-trim-left-f 0 \
  --p-trunc-len-f 245 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 287 \
  --p-n-threads 5 \
  --o-table $TABLE \
  --o-representative-sequences $SEQS \
  --o-denoising-stats $STATS \
  --verbose

#### 04-Visualize_Stats_SHSK.sh
#!/bin/bash
#SBATCH --job-name="Visuzalize_Stats_18S_SHSK"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e Visuzalize_Stats_18S_SHSK.err-%N
#SBATCH -o Visuzalize_Stats_18S_SHSK.out-%N

#Path Variables step #1
INPUTABLE=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureTable_DADA2.qza
OUTPUTABLE=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureTable_DADA2.qzv
MAPPING=/scratch/tjp99569/shsk/qiime-files/18S-Soil-and-SingleNematodes-mapping-01-22-2021.txt

#Path Variables step #2
INPUTSEQ=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2.qza
OUTPUTSEQ=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2.qzv

#Path Variables step #3
INPUTSTATS=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_Denoise_Stats_DADA2.qza
OUTPUTSTATS=/home/tjp99569/shsk/dada2-results-18S/

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Visualize Feature Table Summary, step #1
qiime feature-table summarize \
  --i-table $INPUTABLE \
  --o-visualization $OUTPUTABLE \
  --m-sample-metadata-file $MAPPING \

#Visualize Rep-Seqs, step #2
qiime feature-table tabulate-seqs \
  --i-data $INPUTSEQ \
  --o-visualization $OUTPUTSEQ \

#Export Stats, step #3
  qiime tools export \
  --input-path $INPUTSTATS \
  --output-path $OUTPUTSTATS \

#### 05-Export_FeatureTable_SHSK.sh
#!/bin/bash
#SBATCH --job-name="Export_FeatureTable_18S_SHSK"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e Export_FeatureTable_18S_SHSK.err-%N
#SBATCH -o Export_FeatureTable_18S_SHSK.out-%N

#Path Variables
INPUTABLE=/home/tjp99569/shsk/dada2-results/SHSK_trimmed_FeatureTable_DADA2.qza
OUTPUTABLE=/home/tjp99569/shsk/dada2-results/

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Export Feature table
  qiime tools export \
  --input-path $INPUTABLE \
  --output-path $OUTPUTABLE

#### 06-Assign_taxonomy_SHSK.sh
#!/bin/bash
#SBATCH --job-name="Assign_Taxonomy_18S_SHSK"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=24:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e Assign_Taxonomy_18S_SHSK.err-%N
#SBATCH -o Assign_Taxonomy_18S_SHSK.out-%N

#Path Variables
INPUT=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2.qza
REFSEQ=/scratch/tjp99569/shsk/database/CustomDatabase_SILVA_and_SingleNem_GenBank_Sequences_Nov_06_2021.qza
REFTAX=/scratch/tjp99569/shsk/database/CustomDatabase_SILVA_and_SingleNem_GenBank_Taxonomy_Nov_06_2021.qza
OUTPUT=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2_BLAST_Taxonomy.qza

#Activate QIIME2/2020.11
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.11

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Blast+ parameters
qiime feature-classifier classify-consensus-blast \
--i-query $INPUT \
--i-reference-reads $REFSEQ \
--i-reference-taxonomy $REFTAX \
--p-perc-identity 0.9 \
--o-classification $OUTPUT \
--verbose

##### 07-Taxonomy_Boxplots_SHSK.sh
#!/bin/bash
#SBATCH --job-name="Taxonomy_Boxplots_18S_SHSK"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=24:00:00
#SBATCH --mail-user=tiago.pereira@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e Taxonomy_Boxplots_18S_SHSK.err-%N
#SBATCH -o Taxonomy_Boxplots_18S_SHSK.out-%N

#Path Variables
INPUT_TAX=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2_BLAST_Taxonomy.qza
OUTPUT_TAX=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2_BLAST_Taxonomy.qzv
TABLE=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureTable_DADA2.qza
MAPPING=/scratch/tjp99569/shsk/qiime-files/18S-Soil-and-SingleNematodes-mapping-01-22-2021.txt
BOXPLOT=/home/tjp99569/shsk/dada2-results-18S/SHSK_18S_notrimmed_FeatureSeq_DADA2_BLAST_Taxonomy_BoxPlots.qzv

#Activate QIIME2/2020.11
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.11

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Visualize Taxonomy
qiime metadata tabulate \
  --m-input-file $INPUT_TAX \
  --o-visualization $OUTPUT_TAX \

#Generate boxplots
qiime taxa barplot \
  --i-table $TABLE \
  --i-taxonomy $INPUT_TAX \
  --m-metadata-file $MAPPING \
  --o-visualization $BOXPLOT \
  --verbose