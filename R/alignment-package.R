##################################################
#' Data Alignment Tools
##################################################
#'
#' Methods for the alignment and merging of partially or
#' inexactly aligned time-series.
#'
#' The main \code{alignment} functions are named \code{[type]_align} and
#' in intended to be used in the form:
#'
#' \code{ans <- _align(ts1, ts2, ...)}
#'
#' In the simplest cases \code{ts1} is a \code{vector} to be used as a
#' reference time-series and \code{ts2} is a second time-series vector to be aligned
#' with it. However, either can also be a \code{data.frame} in which case an
#' extra argument \code{by} should be included to identify the data.frame columns
#' to use, e.g. for \code{df1} containing \code{ts1} and \code{df2} containing
#' \code{ts2}:
#'
#' \code{ans <- _align(df1, df2, by=c("ts1", "ts2"), ...)}
#'
#' If the object to be aligned is a \code{data.frame}, all columns in that
#' \code{data.frame} are modified according to the \code{ts1/ts2} alignment.
#'
#' \code{alignment}s can also be applied within a single data, e.g. for
#' \code{df} that contains both \code{ts1} and \code{ts2}:
#'
#' \code{ans <- _align(df, by=c("ts1", "ts2"), ...)}
#'
#' (in which case only df$ts2 is aligned using df$ts1 as the reference...)
#'
#' The default output is typically a single \code{data.frame} containing the
#' aligned data, but other outputs may also be generated, e.g. plots and
#' summary reports. This behaviour can be changed using the common argument
#' \code{output}, and options include \code{"plot"}, \code{"summary"},
#' \code{"ans"}, and \code{"alignment"} (all outputs in the package object
#' class, see also \code{\link{alignment.generics}}). Multiple output may also
#' be requested but only the last is caught.
#'
#' Example \code{alignment}s include: \code{\link{n_align}},
#' \code{\link{cor_align}} and \code{\link{cow_align}}.
#'
#' Although the main \code{alignment} functions require two time-series,
#' some use sub-routines to reshape data, and these can also be applied
#' directly.
#'
#' Examples of these include: \code{\link{regularize}} and
#' \code{\link{warp}}.
#'
#' @references
#'
#' Carslaw, D.C., Ropkins, K., Laxen, D., Moorcroft, S., Marner, B. and
#' Williams, M.L., 2008. Near-field commercial aircraft contribution to
#' nitrogen oxides by engine, aircraft type, and airline by
#' individual plume sampling. Environmental science & technology, 42(6),
#' pp.1871-1876. \url{https://doi.org/10.1021/es071926a}.
#'
#' Ropkins, K., Carlsaw, D.C., Goodman, P.S. and Tate, J.E., 2009.
#' Application of non-linear time-alignment and integration methods to
#' environmental time series. TrAC Trends in Analytical Chemistry, 28(3),
#' pp.373-391. \url{https://doi.org/10.1016/j.trac.2008.11.013}.
#'
#' @name alignment-package
#' @aliases alignment alignment-package
"_PACKAGE"
