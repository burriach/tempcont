\name{Temporal contributions on trends using mixed models}
\alias{TempCont}
\title{TempCont: Temporal contributions on trends using mixed models}
\usage{
  temp.cont(model, resp, driver, random, dataframe)
}
\arguments{
  \item{model}{Fitted "lme" model. I must be a mixed model with gaussian distribution, and it is advisable to control for temporal autocorrelation at lag 1 (e.g., corCAR1())}
  \item{resp}{Response continous variable.}
  \item{driver}{Driver of interest. It must be a continous variable.}
  \item{random}{Grouping factor.}
  \item{dataframe}{Dataset used in the analyses.}
}
\description{
When using temporal data, linear models provide parameter estimates that help to understand how, for instance, interannual variability in a response variable is affected by a predictor. However, if the response variable and the predictors present trends, things get a bit more tricky. This methodology is aimed at estimating the effect of the trend in the predictors on the observed trend of the response variable. To do so, our package first calculates the observed trend of the response (slope estimate ± standard error of the slope) in the data data using GLMMs with random slopes and temporal autocorrelation structure (ARMA, p = 1, q = 0). Then it calculates the trend predicted by the model provided by the user and the trend predicted by the same model but maintaining the predictor of interest constant (i.e., using the median values per grouping level, while all other predictors change according to the observations). The difference between the temporal predictions for the whole model and the temporal predictions of the model when one variable is controlled is the average temporal contribution of that predictor to the trend in the response variable. The difference between all individual temporal contributions and the observed trend should considered to be unknown contributions. Finally, the package also calculates the average sensitivities of the response variable to predictor changes by dividing the temporal contributions by the trends of the predictor variables. All errors are calculated using the error-propagation method. See Fernández-Martínez et al., 2017 and 2019 for further information on the methodology.
}
\value{}
\examples{
  library (nlme)
  library (MuMIn)
  library (tempcont)
  data(tempcont_data)
  mod <- lme (trolls ~ bushes + temperature  + rainfall + drought + gremlins + orcs + warlocks, data=tempcont_data, random= ~1|cave, correlation = corCAR1(form=~year|cave), method="ML")
  summary (mod)
  r.squaredGLMM(mod)
  temp.cont(model = mod, resp = "trolls", driver="temperature", random="cave", dataframe = tempcont_data)
  temp.cont(model = mod, resp = "trolls", driver="orcs", random="cave", dataframe = tempcont_data)
}
