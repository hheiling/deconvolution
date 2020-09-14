sbatch -t 2-12 -c 6 --mem=22528 --partition=largenode --wrap="R CMD BATCH --no-save --no-restore \"--args dataset='EGAD00001002671' sam_name='EGAF00001331297'\" step1_get_counts.R step1_get_counts_EGAD00001002671_EGAF00001331297.Rout"


