\name{dykstra_linealDF}
\alias{dykstra_linealDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dykstra_linealDF(DF, A, b, eq = rep("<=", length(b)), I = matrix(rep(1, dim(DF)[1] * dim(DF)[2]), ncol = dim(DF)[2], nrow = dim(DF)[1]), W = diag(length(DF[1, ])))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DF}{
%%     ~~Describe \code{DF} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{eq}{
%%     ~~Describe \code{eq} here~~
}
  \item{I}{
%%     ~~Describe \code{I} here~~
}
  \item{W}{
%%     ~~Describe \code{W} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (DF, A, b, eq = rep("<=", length(b)), I = matrix(rep(1, 
    dim(DF)[1] * dim(DF)[2]), ncol = dim(DF)[2], nrow = dim(DF)[1]), 
    W = diag(length(DF[1, ]))) 
{
    df_check = check_rows_lineal(DF, A, b, eq)
    df_final = data.frame()
    df_proyec = data.frame()
    for (i in 1:length(DF[, 1])) {
        xi = DF[i, ]
        ii = which(I[i, ][1:length(I[i, ])] == 1)
        xi_prima = xi[ii]
        A_prima = as.matrix(A[, ii])
        Wi = t(as.matrix(W[, ii]))[, ii]
        bi_prima = b - as.matrix(A[, -ii]) \%*\% as.numeric(xi[-ii])
        x = dykstra_lineal(xi_prima, A_prima, bi_prima, eq, Wi)
        xi[ii] = x
        df_final = rbind(df_final, xi)
        if (any(xi != DF[i, ])) {
            df_proyec = rbind(df_proyec, xi)
        }
    }
    return(list(df_proyec = df_proyec, df_check = df_check, df_original = DF, 
        df_final = df_final))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
