\name{temp.cont}
\alias{temp.cont}
\title{Temporal contribution}
\usage{
  temp.cont(model, resp, driver, random, dataframe)
}
\arguments{
  \item{model}{}
  \item{resp}{}
  \item{driver}{}
  \item{random}{}
  \item{dataframe}{}
}
\description{
Prints 'temp.cont, world!'.
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
