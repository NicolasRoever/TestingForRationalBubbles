\name{bhargava_statistic}
\alias{bhargava_statistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bhargava_statstic
}
\description{
This function calculates the Bhargava Statistic from Homm, Ulrich, and Jörg Breitung. "Testing for speculative bubbles in stock markets: a comparison of alternative methods." Journal of Financial Econometrics 10.1 (2012): 198-231.
}
\usage{
bhargava_statistic(x, tau_0 = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x is a vector or a time series object containing the relevant time series values.
}
  \item{tau_0}{
  The test statistic is recursively calculated for subsamples of the time series. tau_0 gives the fraction for the smallest subsample to be considered. So if your sample is 100 and tau_0 is 0.2, the smalles subsample consists of 20 data points.
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
Homm, Ulrich, and Jörg Breitung. "Testing for speculative bubbles in stock markets: a comparison of alternative methods." Journal of Financial Econometrics 10.1 (2012): 198-231.
}
\author{
Nicolas Röver
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
