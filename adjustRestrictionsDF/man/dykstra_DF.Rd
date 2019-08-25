\name{dykstra_DF}
\alias{dykstra_DF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dykstra_DF(DF, typeConstraints = "list() or c() of characters", f_proyec = "list of functions returning vec", arguments = "a list of arguments for the function f_proyec", A = "matrix", b = "vector", r = rep(1, 1), centers = matrix(rep(0, length(DF[1, ]) * length(r)), nrow = length(r), ncol = length(DF[1, ])), eq = c(), I = matrix(rep(1, dim(DF)[1] * dim(DF)[2]), ncol = dim(DF)[2], nrow = dim(DF)[1]), W = diag(length(DF[1, ])))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DF}{
%%     ~~Describe \code{DF} here~~
}
  \item{typeConstraints}{
%%     ~~Describe \code{typeConstraints} here~~
}
  \item{f_proyec}{
%%     ~~Describe \code{f_proyec} here~~
}
  \item{arguments}{
%%     ~~Describe \code{arguments} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
function (DF, typeConstraints = "list() or c() of characters", 
    f_proyec = "list of functions returning vec", arguments = "a list of arguments for the function f_proyec", 
    A = "matrix", b = "vector", r = rep(1, 1), centers = matrix(rep(0, 
        length(DF[1, ]) * length(r)), nrow = length(r), ncol = length(DF[1, 
        ])), eq = c(), I = matrix(rep(1, dim(DF)[1] * dim(DF)[2]), 
        ncol = dim(DF)[2], nrow = dim(DF)[1]), W = diag(length(DF[1, 
        ]))) 
{
    if (all(list("lineal", "ball", "other") \%in\% typeConstraints)) {
    }
    else if (all(list("lineal", "ball") \%in\% typeConstraints) & 
        !("other" \%in\% typeConstraints)) {
        if (length(eq) == 0) 
            eq = rep("<=", length(b) + length(r))
        results = dykstra_linealBallDF(DF, A = A, b = b, r = r, 
            centers = centers, eq = eq, I = I, W = W)
    }
    else if ("other" \%in\% typeConstraints) {
        results = dykstra_proyecDF(DF = DF, f_proyec = f_proyec, 
            arguments = arguments, W = W)
    }
    else if ("ball" \%in\% typeConstraints) {
        if (length(eq) == 0) 
            eq = rep("<=", length(r))
        results = dykstra_ballDF(DF = DF, r = r, centers = centers, 
            eq = eq, I = I, W = W)
    }
    else if ("lineal" \%in\% typeConstraints) {
        if (length(eq) == 0) 
            eq = rep("<=", length(b))
        results = dykstra_linealDF(DF = DF, A = A, b = b, eq = eq, 
            I = I, W = W)
    }
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
