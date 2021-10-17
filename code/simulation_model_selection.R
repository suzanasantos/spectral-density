################################################################################
## Simulation - model selection
################################################################################


# Dependencies -----------------------------------------------------------------

library(snowfall)
library(igraph)
source('statGraph.R')

# Simulation setup -------------------------------------------------------------

directory <- "./" # results folder

models <- c("ER", "GRG", "KR", "WS", "BA")
# ER - Erdos-Renyi
# KR - d-regular graphs
# GRG - geometric random graphs
# BA - Barabasi-Albert random graph model
# WS - Watts Strogatz random graph model


maxSim <- 1000 # Number of repetitions of the experiment, that is, number of 
                # graphs generated for estimating the model parameter
maxGraphs <- 100 # Number of Monte Carlo instances of random graphs for 
                 # approximating the asymptotic density function.
n <- 50  # graph size (number of vertices)

p <- 0.3 # parameter used to generare graphs. Por the BA model, the parameter is
         # 1 + p


# Auxiliar function ------------------------------------------------------------

# Generate an instance of a random graph for a given model
generateGraph <- function(model) {
    if (model == "ER") {
        g <- erdos.renyi.game(n, p)
        g <- as.matrix(get.adjacency(g))
    }
    if (model == "GRG") {  ## o TORUS precisa ser FALSE, senao o estimador nao funciona
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
    return(g)
}

# Pre-computing the random graph model spectra in the search grid --------------

Spectra <- list()

for (model in models) {
    eps <- 0.001
    parameters <- seq(0, 1, eps)
    parameters <- parameters[-c(1, length(parameters))]
    
    if (model == "KR")
        parameters <- 1:(n/2)

    if (model == "BA")
        parameters <- seq(0, 4, 0.01)

    set.seed(0)
    spectra <- array(NA, c(length(parameters), n, maxGraphs))
    rownames(spectra) <- as.character(parameters)
    for (parameter in parameters) {
        spectra[as.character(parameter), , ] <- modelSpectra(model, n, parameter, ngraphs=maxGraphs)
    }

    Spectra[[model]] <- spectra

}


executionName <- paste("model_selection_", n, "_p=", p, sep="")

# Parallel execution setup -----------------------------------------------------

sfInit(parallel = TRUE , type ="SOCK", socketHosts = c(rep("localhost", 30)))


# Exporting global variables ---------------------------------------------------

sfExportAll() # exporta todas as variáveis
sfLibrary(igraph) # exporta pacotes a serem utilizados

# Function that runs each instance of the simulation ---------------------------
run <- function(sim) {
    if (!ime)
    system(paste("renice 10 -p", Sys.getpid()))

    result <- list() 

    for (model in models) {
        g <- generateGraph(model)
        result[[model]] <- list()
        result[[model]][["L1"]] <- graph.model.selection.cdf(g, Spectra, dist="L1")
        result[[model]][["KL"]] <- graph.model.selection(g, Spectra, dist="KL")
        result[[model]][["L1density"]] <- graph.model.selection(g, Spectra, dist="L1")
    }
    return(result)
}

# Paraleliza as N iterações em 40 núcleos (10 na máquina local, 12 na darwin e 18 na hulk3).
# Salva os resultados em "~/.sfCluster/restore/SAVE_DEFAULT_BAGULHO"
results <- sfClusterApplySR(1:maxSim, run, name=executionName, restore=TRUE, perUpdate=1)

sfStop()

save(results, file=paste(directory, executionName, ".RData", sep=""))

################################################################################
