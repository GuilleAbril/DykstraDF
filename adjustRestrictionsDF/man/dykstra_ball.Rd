\name{dykstra_ball}
\alias{dykstra_ball}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dykstra_ball(x0, r = rep(1, 1), centers = matrix(rep(0, length(x0) * length(r)), nrow = length(r), ncol = length(x0)), eq = rep("<=", length(r)), W = diag(length(x0)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{centers}{
%%     ~~Describe \code{centers} here~~
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
function (x0, r = rep(1, 1), centers = matrix(rep(0, length(x0) * 
    length(r)), nrow = length(r), ncol = length(x0)), eq = rep("<=", 
    length(r)), W = diag(length(x0))) 
{
    M = length(r)
    N = length(x0)
    E = matrix(rep(0, N * M), nrow = M, ncol = N)
    n = 1
    error_vec = c()
    error = 1
    while ((error > 1e-04 || n <= M) & n < 100) {
        i = 1 + (n - 1)\%\%M
        x0e = as.numeric(x0) + E[i, ]
        x0e_prima = (x0e - centers[i, ])/r[i]
        if (eq[i] == "<=") {
            x = P_ball(x0e_prima, W)
        }
        else {
            x = P_circunference(x0_primae, W)
        }
        x = x * r[i] + centers[i, ]
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
