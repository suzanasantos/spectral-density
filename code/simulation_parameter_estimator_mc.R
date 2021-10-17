################################################################################
## Simulation - parameter estimator using Monte Carlo approximations of density
## function
################################################################################


# Required libraries and scripts -----------------------------------------------

library(snowfall)
library(igraph)
source('statGraph.R')

# Simulation setup -------------------------------------------------------------

directory <- "./" # results folder

model <- "BA" # Random graph model used to generate graphs.
# Possible values are:
# ER - Erdos-Renyi
# KR - d-regular graphs
# GRG - geometric random graphs
# BA - Barabasi-Albert random graph model
# WS - Watts Strogatz random graph model

eps <- 0.001 # search grid precision

maxSim <- 1000  # Number of repetitions of the experiment, that is, number of 
                # graphs generated for estimating the model parameter

maxGraphs <- 100 # Number of Monte Carlo instances of random graphs for 
                 # approximating the asymptotic density function.

n <- 1000  # Graph size


parameters <- seq(0, 1, eps) # search grid
parameters <- parameters[-c(1, length(parameters))]

if (model == "KR")
    parameters <- 1:(n/2)

if (model == "BA")
    parameters <- seq(0, 1, eps)


p <- 0.7

executionName <- paste("cum_vs_density_parameter_estimator2_", model, "_n=", n, "_p=", p, sep="")



# Precomputing the eigenvalues for the Monte Carlo approaximation of the 
# asymptotic density function of the random graph model ------------------------


set.seed(0)
spectra <- array(NA, c(length(parameters), n, maxGraphs))
rownames(spectra) <- as.character(parameters)
for (parameter in parameters) {
    m <- model
    if (laplacian)
        m <- modelLaplacian(model)
    spectra[as.character(parameter), , ] <- modelSpectra(m, n, parameter, ngraphs=maxGraphs)
}

# Parralel execution setup -----------------------------------------------------


sfInit(parallel = TRUE , type ="SOCK", socketHosts = c(rep("localhost", 10)))


# Exporting global variables ---------------------------------------------------

sfExportAll() # exporta todas as variáveis
sfLibrary(igraph) # exporta pacotes a serem utilizados


# Function that fits the raindom graph model for each random graph instance ----

run <- function(sim) {
    system(paste("renice 10 -p", Sys.getpid()))

    if (model == "ER") {
        g <- erdos.renyi.game(n, p)
        g <- as.matrix(get.adjacency(g))
    }
    if (model == "GRG") { 
        g <- grg.game(n, p)
        g <- as.matrix(get.adjacency(g))
    }
    if (model == "KR") {
    	g <- k.regular.game(n, p*10)
        g <- as.matrix(get.adjacency(g))
    }
    if (model == "WS") {
    	g <- watts.strogatz.game(1, n, 2, p)
        g <- as.matrix(get.adjacency(g))
    }
    if (model == "BA") {
    	g <- barabasi.game(n, power=p+1, directed = FALSE)
        g <- as.matrix(get.adjacency(g))
    }
  
    parameter <- c() 


    parameter["L1"] <- graph.param.estimator.cdf(g, spectra, dist="L1")$p

    parameter["L1density"] <- graph.param.estimator(g, spectra, dist="L1")$p
    return(parameter)
}


# Running parallel execution ---------------------------------------------------

results <- sfClusterApplySR(1:maxSim, run, name=executionName, restore=T, perUpdate=1)

sfStop()

save(results, file=paste(directory, executionName, ".RData", sep=""))


# Summarizing the results ------------------------------------------------------

methods <- names(results[[1]])
ps <- matrix(NA, length(results), length(methods))
colnames(ps) <- methods

for (i in 1:length(results)) {
    for (m in methods) {
        ps[i, m] <- as.numeric(results[[i]][m])
    }
}

colMeans(ps)

apply(ps, 2, sd)
