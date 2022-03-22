#!/bin/bash
#SBATCH -J CRCount		# Can be changed
#SBATCH -o CRCount-%j.out	# MODIFY WITH CAUTION
#SBATCH -e CRCount-%j.err	# MODIFY WITH CAUTION
#SBATCH -n 20

# Change the following for your particular experiment
SCRATCHDIR="/home/chinye.nwokolo/scRNA1"
INPDIR="/AIG_11-02/"
OUTDIR="AIG_out"
GENOME="mm10-2020-A"

# EDITS BELOW THIS LINE ONLY WITH CAUTION

ncores="20"
GENOMEDIR="/data/projects/cellranger/"

module load cellranger

cellranger count --id=${OUTDIR} \
		 --transcriptome=${GENOMEDIR}/refdata-gex-${GENOME} \
		 --fastqs=${SCRATCHDIR}/${INPDIR} \
		 --localcores=${ncores} # DO NOT MODIFY

exit 0

