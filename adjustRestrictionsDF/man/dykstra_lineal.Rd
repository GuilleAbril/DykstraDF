\name{dykstra_lineal}
\alias{dykstra_lineal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dykstra_lineal(x0, A, b, eq, W = diag(length(x0)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
%%     ~~Describe \code{x0} here~~
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
function (x0, A, b, eq, W = diag(length(x0))) 
{
    M = length(as.vector(b))
    N = length(x0)
    E = matrix(rep(0, N * M), nrow = M, ncol = N)
    n = 1
    error_vec = c()
    error = 1
    while ((error > 1e-04 || n <= M) & n < 100) {
        i = 1 + (n - 1)\%\%M
        x0e = as.numeric(x0) + E[i, ]
        if (eq[i] == "==") {
            x = P_subspace(x0e, as.matrix(A)[i, ], as.vector(b)[i], 
                W)
        }
        else {
            x = P_Half_space(x0e, as.matrix(A)[i, ], as.vector(b)[i], 
                W)
        }
        e = as.numeric(x0) + E[i, ] - x
        error_vec[i] = calc_error(e, E[i, ], x, x0, W = W)
        error = sum(error_vec)
        E[i, ] = e
        x0 = x
        n = n + 1
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
