\name{dykstra_proyecDF}
\alias{dykstra_proyecDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dykstra_proyecDF(DF, f_proyec = list(), arguments = list(), W = diag(length(arguments[[1]][[1]])))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DF}{
%%     ~~Describe \code{DF} here~~
}
  \item{f_proyec}{
%%     ~~Describe \code{f_proyec} here~~
}
  \item{arguments}{
%%     ~~Describe \code{arguments} here~~
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
function (DF, f_proyec = list(), arguments = list(), W = diag(length(arguments[[1]][[1]]))) 
{
    df_final = data.frame()
    df_proyec = data.frame()
    for (i in 1:length(DF[, 1])) {
        arguments[[1]][[1]] = DF[i, ]
        x = dykstra_proyec(f_proyec = f_proyec, arguments = arguments, 
            W = W)
        df_final = rbind(df_final, x)
        if (any(x != DF[i, ], na.rm = TRUE)) {
            df_proyec = rbind(df_proyec, x)
        }
    }
    colnames(df_final) = colnames(DF)
    colnames(df_proyec) = colnames(DF)
    return(list(df_proyec = df_proyec, df_original = DF, df_final = df_final))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
