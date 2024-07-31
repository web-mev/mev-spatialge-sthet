suppressMessages(library(spatialGE))
#suppressMessages(library(reshape2))
suppressMessages(library(rjson))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
source('/usr/local/bin/prep_stlist.R')

# args from command line:
args <- commandArgs(TRUE)

option_list <- list(
    make_option(
        c('-f', '--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-c', '--coordinates_file'),
        help='Path to the barcode spatial coordinates input.'
    ),
    make_option(
        c('-s', '--sample_name'),
        help='Sample name'
    ),
    make_option(
        c('-n', '--normalization'),
        help='Normalization method of `log` or `sct`'
    ),
    make_option(
        c('-m', '--method'),
        help='The spatial statistic to compute. One of `moran` or `geary`'
    ),
    make_option(
        c('-g', '--gene_set'),
        help='The set of genes to test under SThet.'
    ),
    make_option(
        c('-x', '--xpos_col'),
        help='The column header for the x-position coordinate metadata'
    ),
    make_option(
        c('-y', '--ypos_col'),
        help='The column header for the y-position coordinate metadata'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Sanity checks on the inputs:
# Check that the file was provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$coordinates_file)){
    message('Need to provide a count matrix with the -c/--coordinates_file arg.')
    quit(status=1)
}

genes_to_test <- strsplit(opt$gene_set, ',')[[1]]

# transform the name of the normalization scheme:
if (is.null(opt$normalization)){
    message('Need to provide a normalization scheme with the -n/--normalization arg.')
    quit(status=1)
} else if(tolower(opt$normalization) == 'sctransform'){
    norm_scheme <- 'sct'
} else if(tolower(opt$normalization) == 'log'){
    norm_scheme <- 'log'
} else {
    message('We only accept `log` or `SCTransform` for the normalization scheme.')
    quit(status=1)
}

# transform the name of the test statistic:
if (is.null(opt$method)){
    message('Need to provide a test statistic with the -m/--method arg.')
    quit(status=1)
} else if(opt$method == "Moran's I"){
    stat_method <- 'moran'
} else if(opt$method == "Geary's C"){
    stat_method <- 'geary'
} else {
    message("We only accept `Moran's I` or `Geary's C` for the normalization scheme.")
    quit(status=1)
}

working_dir <- dirname(opt$input_file)
setwd(working_dir)

# prepare an STList instance. 
spat_list <- prep_stlist(opt$input_file, 
                         opt$coordinates_file,
                         opt$sample_name,
                         opt$xpos_col,
                         opt$ypos_col)
spat <- spat_list$spat

# normalize
spat <- transform_data(spat, method=norm_scheme)

# check requested genes against this object:
raw_counts <- spat@counts[[opt$sample_name]]
diff_set <- setdiff(genes_to_test, rownames(raw_counts))
num_invalid_genes <- length(diff_set)
if (num_invalid_genes > 0) {
    max_print <- 5
    if (num_invalid_genes < max_print) {
        invalid_genes <- paste(diff_set[1:num_invalid_genes], collapse=', ')
    } else {
        invalid_genes <- sprintf('%s, and %d others.', paste(diff_set[1:max_print], collapse=', '), num_invalid_genes-max_print)
    }
    message(sprintf('The following genes to test were not found in your count matrix: %s', invalid_genes))
    quit(status=1)
}

# Run SThet
spat <- SThet(
    spat,
    genes = genes_to_test,
    method = stat_method
)

# Pull dataframe from STList metadata and write to table
# skips sample column (the first)
stat_col <- if(stat_method == 'moran') 'moran_i' else 'geary_c'

out_table <- get_gene_meta(
    spat, 
    sthet_only=T
)[, c('gene', 'gene_mean', 'gene_stdevs', stat_col)] %>% drop_na()

if(dim(out_table)[1] > 0){
    output_filename <- paste(working_dir, 'sthet_output.tsv', sep='/')
    write.table(
        out_table,
        output_filename,
        sep="\t",
        quote=F,
        row.names=F
    )
    json_str = paste0('{"SThet_results":"', output_filename, '"}')
    output_json <- paste(working_dir, 'outputs.json', sep='/')
    write(json_str, output_json)
} else {
    message("Note that the results table was empty. The input gene was indeed found in your abundance matrix, but the analysis was unable to successfully complete and report the spatial statistic. Sometimes, this is due to the gene having zero counts across all sampling areas.")
    quit(status=1)
}
