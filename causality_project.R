# Load the libraries
# To install pcalg library you may first need to execute the following commands:
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
biocLite("RBGL")
biocLite("Rgraphviz")

require(vars)
require(urca)
require(pcalg)

# Read the input data

mydata <- read.csv("data.csv")

# Build a VAR model
criteria <- VARselect(mydata, lag.max = 10)
lagvalue <- criteria$selection["SC(n)"]
varmodel <- VAR(mydata, lagvalue)

# Select the lag order using the Schwarz Information Criterion with a maximum lag of 10
# see ?VARSelect to find the optimal number of lags and use it as input to VAR()

# Extract the residuals from the VAR model
# see ?residuals

res <- residuals(varmodel)

# Check for stationarity using the Augmented Dickey-Fuller test
# see ?ur.df

adf <- apply(res, 2, ur.df)

lapply(adf, summary)

# Check whether the variables follow a Gaussian distribution
# see ?ks.test

ks <- apply(res, 2, ks.test, "pnorm")

# Write the residuals to a csv file to build causal graphs using Tetrad software

write.csv(res, "residuals.csv", row.names = FALSE)

# OR Run the PC and LiNGAM algorithm in R as follows,
# see ?pc and ?LINGAM

# PC Algorithm
suffStat = list(C = cor(res), n = nrow(res))
pcmodel <-
  pc(
    suffStat = suffStat,
    indepTest = gaussCItest,
    alpha = 0.1,
    labels = colnames(res),
    verbose = TRUE,
    skel.method = "original"
  )
plot(pcmodel)

# LiNGAM Algorithm
lingam(res, verbose = TRUE)
