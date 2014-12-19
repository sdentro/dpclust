LIBPATH=$1
DATASET=$2
ITERS=$3
BURNIN=$4
DATPATH=$5
PURITYFILE=$6
ANALYSIS=$7
PARALLEL=$8
THREADS=$9
BINSIZE=${10}
NBLOCKS=${11}
CMD="Rscript ${LIBPATH}/RunDP_pipeline.R"
${CMD} ${LIBPATH} ${DATASET} ${ITERS} ${BURNIN} ${DATPATH} ${PURITYFILE} "cons" ${PARALLEL} ${THREADS} ${BINSIZE} "1" ${NBLOCKS}
