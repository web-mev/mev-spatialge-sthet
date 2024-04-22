suppressMessages(library(spatialGE))
suppressMessages(library(reshape2))
suppressMessages(library(rjson))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
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
        c('-g', '--gene_set'),
        help='The set of genes to test under SThet.'
    ),
    make_option(
        c('-o', '--output'),
        help='The name of the output file'
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

# ToDo: exception checking for gene set input


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

# prepare an STList instance. Note that we are potentially re-mapping
# the original gene identifiers to symbols such that they work with
# the MSigDB files:
spat_list <- prep_stlist(opt$input_file, 
                         opt$coordinates_file,
                         opt$sample_name,
                         gene_mapping_df,
                         gene_ids,
                         'SYMBOL')
spat <- spat_list$spat

# normalize
spat <- transform_data(spat, method=opt$normalization)

# Run SThet
spat <- SThet(
    spat,
    genes = opt$gene_set,
    method = 'moran'
)

# Pull dataframe from STList metadata and write to table
# skips gene column
# note: there are two column for data in sthet output, moran and geary;
#   and both are input params to SThet under method
#   data is placed into appropriate column
#   will need logic for switching column if we expose moran vs geary method choice
write.table(
    get_gene_meta(
        spat, 
        sthet_only=T
    )[, c(2,3,4,5)],
    opt$output,
    sep="\t",
    quote=F,
    row.names=F
)