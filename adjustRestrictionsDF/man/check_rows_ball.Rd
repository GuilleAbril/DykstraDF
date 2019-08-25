\name{check_rows_ball}
\alias{check_rows_ball}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
check_rows_ball(DF, r = rep(1, 1), centers = matrix(rep(0, length(DF[1, ]) * length(r)), nrow = length(r), ncol = length(DF[1, ])), eq = rep("<=", length(r)), W = diag(length(DF[1, ])))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DF}{
%%     ~~Describe \code{DF} here~~
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
function (DF, r = rep(1, 1), centers = matrix(rep(0, length(DF[1, 
    ]) * length(r)), nrow = length(r), ncol = length(DF[1, ])), 
    eq = rep("<=", length(r)), W = diag(length(DF[1, ]))) 
{
    wrong_rows = data.frame()
    for (i in 1:length(DF[, 1])) {
        xi = DF[i, ]
        checks_vec = c()
        for (k in 1:length(r)) {
            r_prima = (as.numeric(xi) - centers[k, ]) \%*\% W \%*\% 
                (as.numeric(xi) - centers[k, ])
            if (eq[k] == "<=") {
                checks_vec[k] = as.vector(r_prima) <= r[k]
            }
            else {
                checks_vec[k] = as.vector(r_prima) == r[k]
            }
        }
        if (any(checks_vec == FALSE, na.rm = TRUE)) {
            wrong_rows = rbind(wrong_rows, xi)
        }
    }
    return(wrong_rows)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
