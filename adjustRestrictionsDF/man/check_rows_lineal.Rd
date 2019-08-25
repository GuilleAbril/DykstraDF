\name{check_rows_lineal}
\alias{check_rows_lineal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
check_rows_lineal(df, A, b, eq = rep("<=", length(b)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
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
function (df, A, b, eq = rep("<=", length(b))) 
{
    wrong_rows = data.frame()
    for (i in 1:length(df[, 1])) {
        xi = df[i, ]
        b_prima = A \%*\% as.numeric(xi)
        checks_vec = c()
        for (k in 1:length(b)) {
            if (eq[k] == "<=") {
                checks_vec[k] = as.vector(b_prima)[k] <= b[k]
            }
            else {
                checks_vec[k] = as.vector(b_prima)[k] == b[k]
            }
        }
        if (any(checks_vec == FALSE)) {
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
