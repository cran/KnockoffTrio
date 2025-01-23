#' Example data for KnockoffTrio
#'
#' A toy example of haplotype and genotype data for trios and duos.
#'
#' @format KnockoffTrio.example contains the following items:
#' \describe{
#'   \item{trio}{A 9*5 numeric genotype matrix of 3 trios and 5 variants. Each trio contains 3 rows in the order of father, mother and offspring. Each column represents a variant.}
#'   \item{trio.hap}{A 18*5 numeric haplotype matrix of 3 trios and 5 variants. Each trio contains 6 rows in the order of father, mother and offspring. Each column represents a variant.}
#'   \item{duo.hap}{A 12*5 numeric haplotype matrix of 3 duos and 5 variants. Each duo contains 4 rows in the order of a single parent and offspring. Each column represents a variant.}
#'   \item{pos}{A numeric vector of length 5 for the position of 5 variants.}
#' }
"KnockoffTrio.example"
