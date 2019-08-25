\name{dykstra_ballDF}
\alias{dykstra_ballDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dykstra_ballDF(DF, r = rep(1, 1), centers = matrix(rep(0, length(DF[1, ]) * length(r)), nrow = length(r), ncol = length(DF[1, ])), eq = rep("<=", length(r)), I = matrix(rep(1, length(DF[1, ]) * length(DF[, 1])), ncol = length(DF[1, ]), nrow = length(DF[, 1])), W = diag(length(DF[1, ])))
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
function (DF, r = rep(1, 1), centers = matrix(rep(0, length(DF[1, 
    ]) * length(r)), nrow = length(r), ncol = length(DF[1, ])), 
    eq = rep("<=", length(r)), I = matrix(rep(1, length(DF[1, 
        ]) * length(DF[, 1])), ncol = length(DF[1, ]), nrow = length(DF[, 
        1])), W = diag(length(DF[1, ]))) 
{
    df_check = check_rows_ball(DF = DF, r = r, centers = centers, 
        eq = eq, W = W)
    df_final = data.frame()
    df_proyec = data.frame()
    for (i in 1:length(DF[, 1])) {
        xi = DF[i, ]
        ii = which(I[i, ][1:length(I[i, ])] == 1)
        if (any(I[i, ] == 1)) {
            if (all(I[i, ] == 1)) {
                x = dykstra_ball(xi, r, centers, eq, W)
            }
            else {
                xi_prima = xi[ii]
                centers_prima = t(centers[, ii])
                r_prima = c()
                for (k in 1:length(r)) {
                  r_prima[k] = sqrt(abs(r^2 - as.numeric(xi[-ii]) \%*\% 
                    as.matrix(W[, -ii])[-ii, ] \%*\% as.numeric(xi[-ii])))
                }
                Wi = t(as.matrix(W[, ii]))[, ii]
                x = dykstra_ball(xi_prima, r_prima, centers_prima, 
                  eq, Wi)
            }
            xi[ii] = x
            df_final = rbind(df_final, xi)
            if (any(xi != DF[i, ], na.rm = TRUE)) {
                df_proyec = rbind(df_proyec, xi)
            }
        }
        else {
            df_final = rbind(df_final, xi)
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
