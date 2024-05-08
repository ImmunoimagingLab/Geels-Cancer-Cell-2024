#!/bin/bash
#SBATCH --job-name=Salmon_Quant      ## Name of the job.
#SBATCH -A FMARANGO_LAB     ## account to charge (all labs have a separate account, if not request HPC)
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=8    ## number of cores the job needs
#SBATCH --mem=48G    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH -t 7-00:00:00           ## 7 day run time
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=rssousa@uci.edu


INPATH=/share/crsp/lab/fmarango/rssousa/PD1RNAseq/EGAData    # directory path
BASENAME=$1               # Basename of fastq file
OUTDIR=/share/crsp/lab/fmarango/rssousa/PD1RNAseq/Gide2019/ENAfastq/POST

IN_FQ1=${BASENAME}/*_R1.fastq.gz  # first read of pair fastQ
IN_FQ2=${BASENAME}/*_R2.fastq.gz  # second read of pair fastQ


SALMON_BIN=${INPATH}/scripts/salmon-latest_linux_x86_64/bin/salmon
SALMON_IDX=${INPATH}/data/Annotation/Salmon_Index/GRCh38_Gencode32                            # path to salmon index directory
GTF_FILE=${INPATH}/data/Annotation/Salmon_Index/GRCh38_Gencode32/gencode.v32.annotation.gtf   # gtf file used with the salmon index
THREADS=8 # number of threads


###### Salmon quantification

OUT_FILE_DONE=${BASENAME}.SalmonQuant.done
cd $OUTDIR

echo -e "Starting Salmon Quantification\n"


if [ ! -f $OUT_FILE_DONE ]; then

    ${SALMON_BIN} quant \
       -i         $SALMON_IDX \
       -o         ${BASENAME}_SalmonOut \
       -l         A \
       -p         $THREADS \
       -1         $IN_FQ1 \
       -2         $IN_FQ2 \
       -g         $GTF_FILE \
       --validateMappings \
       && touch $OUT_FILE_DONE

fi

echo -e "Script Complete \n\n"
