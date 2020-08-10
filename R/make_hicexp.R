#Edits from Monash for joining data.table instead of dlpyr
#' Make Hi-C experiment object from data
#' 
#' @export
#' @importFrom dplyr full_join left_join right_join
#' @import data.table
#' 
#' @param ... Hi-C data. Data must in sparse upper triangular 
#'     format with 4 columns: chr, region1, region2, IF or
#'     in 7 column BEDPE format with columns chr, start1, 
#'     end1, chr, start2, end2, IF.
#' @param data_list Alternate way to enter data. If you have
#'     your Hi-C data in the form of a list already with each
#'     entry of the list representing a sample use this option.
#' @param groups A vector of the experimental groups 
#'     corresponding to each Hi-C data object entered.
#'     If it is not in factor form when entered it will
#'     be converted to a factor.
#' @param covariates Optional data.frame containing 
#'     covariate information for your Hi-C experiment.
#'     Some examples are enzyme used, batch number, etc.
#'     Should have the same number of rows as the number
#'     of Hi-C data objects entered and columns corresponding
#'     to covariates. 
#' @param remove_zeros Logical, should rows with 1 or more
#'     zero IF values be removed?
#' @param zero.p The proportion of zeros in a row to filter by. 
#'     If the proportion of zeros in a row is <= zero.p
#'     the row will be filtered out, i.e. zero.p = 1 means
#'     nothing is filtered based on zeros and zero.p = 0 
#'     will filter rows that have any zeros.
#' @param A.min The minimum average expression value
#'     (row mean) for an interaction pair. If the 
#'     interaction pair has an average expression 
#'     value less than A.min the row will be filtered
#'     out.
#' @param filter Logical, should filtering be performed?
#'     Defaults to TRUE. If TRUE it will filter out
#'     the interactions that have low average IFs
#'     or large numbers of 0 IF values. As these
#'     interactions are not very interesting and
#'     are commonly false positives during difference
#'     detection it is better to remove them from
#'     the dataset. Additionally, filtering will
#'     help speed up the run time of multiHiCcompare. 
#'     Filtering can be performed before or after 
#'     normalization, however the best computational
#'     speed gain will occur when filtering is done
#'     before normalization. Filtering parameters
#'     are controlled by the zero.p and A.min options.
#' @param remove.regions A GenomicRanges object indicating
#'     specific regions to be filtered out. By default
#'     this is the hg19 centromeric, gvar, and stalk
#'     regions. Also included in the package is
#'     hg38_cyto. If your data is not hg19 you will 
#'     need to substitute this file. To choose not 
#'     to filter any regions set regions = NULL. NOTE:
#'     if you set filter = FALSE these regions will NOT
#'     be removed. This occurs in conjuction with the 
#'     filtering step.
#' @details Use this function to create a hicexp object for
#'     analysis in multiHiCcompare. Filtering can also be 
#'     performed in this step if the filter option is 
#'     set to TRUE. Filtering parameters are controlled
#'     by the zero.p and A.min options.
#'     
#' @return A hicexp object.
#' @examples 
#' # load data in sparse upper triangular format
#' data("HCT116_r1", "HCT116_r2", "HCT116_r3", "HCT116_r4", 
#'     "HCT116_r5", "HCT116_r6")
#' # make groups & covariate input
#' groups <- factor(c(1, 1, 1, 2, 2, 2))
#' covariates <- data.frame(enzyme = factor(c('mobi', 'mboi', 'mboi',
#'  'dpnii', 'dpnii', 'dpnii')), batch = c(1, 2, 1, 2, 1, 2))
#' # make the hicexp object
#' hicexp <- make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4,
#'      HCT116_r5, HCT116_r6, groups = groups, 
#'      covariates = covariates)


make_hicexp <- function(..., data_list = NA, groups, covariates = NULL, 
                        zero.p = 0.8, A.min = 5, filter = TRUE, remove.regions = hg19_cyto) {
  if (!is.na(data_list[1])) {
    tabs <- data_list
  } else {
    tabs <- list(...)
  }
  # check number of columns of input data
  if(ncol(tabs[[1]]) != 4 &  ncol(tabs[[1]]) != 7) {
    stop("You must enter data in 4 column sparse matrix format
         or in 7 column BEDPE format.")
  }
  # if BEDPE pull out only columns we need
  if (ncol(tabs[[1]]) == 7) {
    tabs <- lapply(tabs, function(x) {
      colnames(x) <- c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'IF')
      new.x <- x[, c('chr1', 'start1', 'start2', 'IF')]
      return(new.x)
    })
  }
  
  # convert groups to a factor if it is not already
  if (!is.factor(groups)) {
    groups <- as.factor(groups)
  }

  # initialize values
  grp_tbl <- table(groups)
  grp_singletons <- names(grp_tbl[grp_tbl <= 1])
  grp_logical <- groups %in% grp_singletons
  groups <- droplevels(groups[!grp_logical])

  if(length(grp_singletons) > 0) {
    grp_str <- paste(grp_singletons, collapse = "\n")
    warning(paste0("These groups have been dropped since they don't have replicates:\n", grp_str))
  }

  uniq_groups <- unique(groups)
  res <- vector(length = length(uniq_groups))
  
  #keep samples that have replicates
  tabs <- tabs[!grp_logical]

  # set column names of data and cast data.table
  tabs <- lapply(1:length(tabs), function(i) {
                     colnames(tabs[[i]]) <- c("chr", "region1", "region2", paste0("IF", i))
                     as.data.table(tabs[[i]]) })

  # check for correct input
  if (length(groups) != length(tabs)) {
    stop("Length of groups must equal the number of Hi-C data objects entered")
  }
  #NOTE this is similar to what I did above with grp_singletons I think..
  #if (sum(grp_tbl < 2) > 0) {
  #  warning("Each experimental condition should have at least 2 samples. 
  #          If you have less than 2 samples per group use HiCcompare instead")
  #}
  if (!is.null(covariates)) {
    # check for data.frame
    if (!is(covariates, "data.frame")) {
      stop("Enter a data.frame for covariates")
    }
    if (nrow(covariates) != length(tabs)) {
      stop("Number of rows in covariates should correspond to number 
           of Hi-C data objects entered")
    }
    }
  if (zero.p < 0 | zero.p > 1) {
    stop("zero.p must be in [0,1]")
  }
  if (A.min < 0) {
    stop("A.min must be >= 0")
  }
  
  # cycle through groups to create a table for each experimental condition
  print("CUSTOM")
  hic_tbl <- Reduce(function(x, y) { merge.data.table(x, y, by = c("chr", "region1", "region2"), no.dups = FALSE, all = TRUE) }, tabs)
  hic_tbl[is.na(hic_tbl)] <- 0
  
  # calculate resolution
  bins <- unique(c(hic_tbl$region1, hic_tbl$region2))
  bins <- bins[order(bins)]
  resolution <- min(diff(bins))
  hic_tbl[, `:=`(D, abs(region2 - region1) / resolution)]

  # rearrange columns
  setcolorder(hic_tbl, c('chr', 'region1', 'region2', 'D'))
  
  # check chr format, if "chr#" change to just the number
  if (!is.numeric(hic_tbl$chr)) {
    # replace any "chr"
    hic_tbl[, chr := sub("chr", "", chr)]
    # replace any X 
    hic_tbl[, chr := sub("X", "23", chr)]
    # replace any Y
    hic_tbl[, chr := sub("Y", "24", chr)]
    # convert to numeric
    hic_tbl[, chr := as.numeric(chr)]
  }
  # sort hic_tbl
  #TODO do we really need to do this ? I wonder if it would be faster to do it on individual samples
  #  and then outer join. I think outer join should preserve the sorting order..
  setorder(hic_tbl, chr, region1, region2)
  
  # check that all resolutions are equal
  if (length(unique(res)) > 1) {
    stop("Resolution of all datasets must be equal.")
  }
  
  # make metadata 
  metadata <- data.frame(group = groups)
  row.names(metadata) <- paste0("Sample", 1:length(groups))
  if (!is.null(covariates)) {
    metadata <- cbind(metadata, covariates)
  }
  
  # put into hicexp object
  experiment <- new("Hicexp", hic_table = hic_tbl, 
                    comparison = data.table::data.table(), 
                    metadata = metadata, resolution = resolution, 
                    normalized = FALSE)
  
  # filtering steps
  if(filter) {
    experiment <- hic_filter(experiment, zero.p = zero.p, A.min = A.min)
  }
  
  if(!is.null(remove.regions)) {
    experiment <- hic_filter_regions(experiment, remove.regions = remove.regions)
  }

  return(experiment)
}
