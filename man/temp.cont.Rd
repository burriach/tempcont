\name{temp.cont}
\alias{temp.cont}
\alias{TempCont}
\title{TempCont: Temporal contributions on trends using mixed models}
\usage{
  temp.cont(model, driver, random, timevar="year", plot=TRUE)
}
\arguments{
  \item{model}{Fitted "lme" model. I must be a mixed model with gaussian distribution, and it is advisable to control for temporal autocorrelation at lag 1 (e.g., corCAR1())}
  \item{driver}{Driver of interest. It must be a continous variable.}
  \item{random}{Grouping factor.}
  \item{plot}{if \code{TRUE}, plot the data and the predictions.}
  \item{timevar}{name of the time variable.}
}
\description{
When using temporal data, linear models provide parameter estimates that help to understand how, for instance, interannual variability in a response variable is affected by a predictor. However, if the response variable and the predictors present trends, things get a bit more tricky. This methodology is aimed at estimating the effect of the trend in the predictors on the observed trend of the response variable. To do so, our package first calculates the observed trend of the response (slope estimate ± standard error of the slope) in the data data using GLMMs with random slopes and temporal autocorrelation structure (ARMA, p = 1, q = 0). Then it calculates the trend predicted by the model provided by the user and the trend predicted by the same model but maintaining the predictor of interest constant (i.e., using the median values per grouping level, while all other predictors change according to the observations). The difference between the temporal predictions for the whole model and the temporal predictions of the model when one variable is controlled is the average temporal contribution of that predictor to the trend in the response variable. The difference between all individual temporal contributions and the observed trend (should be calculated separately) should considered to be unknown contributions. Finally, the package also calculates the average sensitivities of the response variable to predictor changes by dividing the temporal contributions by the trends of the predictor variables. All errors are calculated using the error-propagation method. See Fernández-Martínez et al., 2017 and 2019 for further information on the methodology.
}
\value{
A data.frame with the following variables for each driver variable. Row names correspon to the the name of the variable.
  \item{mod.slope}{Trend of the response variable predicted by the model.}
  \item{mod.slope.se}{Standard error of mod.slope.}
  \item{mod.slope.t}{t-value of mod.slope}
  \item{mod.slope.p}{P-value of mod.slope}
  \item{temp.contr}{Temporal contribution of the predictor}
  \item{temp.contr.se}{Standard error of temp.contr}
  \item{temp.contr.t}{t-value of temp.contr}
  \item{temp.contr.p}{P-value of temp.contr}
  \item{pred.trend}{Trend of the predictor}
  \item{pred.trend.se}{Standard error of pred.trend}
  \item{pred.trend.t}{t-value of pred.trend}
  \item{pred.trend.p}{P-value of pred.trend}
  \item{sensit}{Trend sensitivity of the response variable to the predictor}
  \item{sensit.se}{Standard error of sensit}
  \item{sensit.t}{t-value of sensit}
  \item{sensit.p}{P-value of sensit}
}
\references{
Fernández-Martínez et al., 2017. "Atmospheric deposition, CO2, and change in the land carbon sink".  Scientific Reports \bold{7}:9632.

Fernández-Martínez et al., 2019. "Global trends in carbon sinks and their relationships with CO2 and temperature".  Nature Climate Change \bold{9}:73–79.
}
\author{Marcos Fernández-Martínez, Joan Maspons}

\examples{
  library (nlme)
  data(tempcont_data)
  mod <- lme (trolls ~ bushes + temperature  + rainfall + drought + gremlins +
              orcs + warlocks, data=tempcont_data, random= ~1|cave,
              correlation = corCAR1(form=~year|cave), method="ML")
  summary (mod)

  temp.cont(model = mod, driver="temperature", random="cave")
  \donttest{temp.cont(model = mod, random="cave")} % Avoid long runtimes when checking the package

}
