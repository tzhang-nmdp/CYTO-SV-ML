###############################################################################
# interface_data.R
#
# Prepares SV prediction data for the R Shiny web portal.
# Reads the combined prediction file (.all_pred) and cytoband coordinates,
# recodes numeric labels to human-readable strings, and saves as RData.
#
# Usage:
#   Rscript interface_data.R -i <input_prefix> -o <output_rdata> -c <cytoband_file>
#
# Arguments:
#   -i : Input file prefix (reads {prefix}.all_pred)
#   -o : Output RData file path
#   -c : Cytoband dictionary CSV (chromosome centromere coordinates)
#
# Output:
#   RData file containing:
#     all_sv      - data.frame with SV predictions and metadata
#     genome_coor - data.frame with chromosome centromere coordinates
###############################################################################

# ---- Load required packages ----
requiredPackages <- c('optparse')
for (p in requiredPackages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p)
    }
}
library(optparse)

# ---- Parse command-line arguments ----
option_list <- list(
    make_option(c("-i", "--input_file"),
        type = "character", help = "input file prefix",
        default = NA, metavar = "filename"
    ),
    make_option(c("-o", "--output_file"),
        type = "character", help = "output RData file",
        default = NA, metavar = "filename"
    ),
    make_option(c("-c", "--cyto_band_file"),
        type = "character", help = "cytoband dictionary file",
        default = NA, metavar = "filename"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input <- as.character(opt$input_file)
output <- as.character(opt$output_file)
cyto_band_file <- as.character(opt$cyto_band_file)

# ---- Read input data ----
# Combined TRS + nonTRS prediction file
all_sv <- read.table(
    paste(input, '.all_pred', sep = ''),
    sep = '\t', header = TRUE
)

# Chromosome centromere coordinates (for genome visualization)
genome_coor <- read.table(cyto_band_file, sep = '\t', header = TRUE)

# ---- Recode numeric labels to human-readable strings ----
# Database benchmark labels:
#   0  -> UC (Unclassified)
#   -1 -> TA (True Artifact)
#   1  -> TG (True Germline)
#   2  -> TS (True Somatic)
all_sv$label <- as.character(all_sv$label)
all_sv[all_sv['label'] == '0', 'label'] <- 'UC'
all_sv[all_sv['label'] == '-1', 'label'] <- 'TA'
all_sv[all_sv['label'] == '1', 'label'] <- 'TG'
all_sv[all_sv['label'] == '2', 'label'] <- 'TS'

# Model prediction labels (same mapping)
all_sv$predict_label <- as.character(all_sv$predict_label)
all_sv[all_sv['predict_label'] == '-1', 'predict_label'] <- 'TA'
all_sv[all_sv['predict_label'] == '1', 'predict_label'] <- 'TG'
all_sv[all_sv['predict_label'] == '2', 'predict_label'] <- 'TS'

# ---- Save as RData for Shiny app ----
save(all_sv, genome_coor, file = output)
