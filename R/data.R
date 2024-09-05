#' Data for infoset function
#'
#' Contains daily prices of ETFs
#'
#' @format A data frame with 3174 rows and 44 columns
#'
#' @source {Created in-house to serve as an example}
#'
#' @examples
#' data(sample.data)
"sample.data"

#' Data with time points for portfolio construction using the LR_cp measure
#'
#' Contains daily prices of ETFs
#'
#' @format A data frame with 3175 rows and 45 columns
#'
#'
#' @source {Created in-house to serve as an example}
#'
#' @examples
#' data(sample.data.ts)
"sample.data.ts"

#' Data for clustering and labeling ETFs
#'
#' Contains asset class of ETFs
#'
#' @format A data frame with 44 observations (rows) on 3 variables (columns)
#' \describe{
#' \item{id}{name of ETF}
#' \item{label}{from 1 to 5 according to the specific asset class}
#' \item{class}{specific asset class (5 categories)}
#' }
#'
#' @source {Created in-house to serve as an example}
#'
#' @examples
#' data(asset.label)
"asset.label"
