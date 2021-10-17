################################################################################
## Simulations - parameter estimator based on the analytical density function
################################################################################


# Required libraries and scripts -----------------------------------------------

library(snowfall)
library(igraph)

source('statGraph.R')

# Simulation setup -------------------------------------------------------------

directory <- "./" # results folder

parallel <- TRUE # indicates whether simulatio will be run in parallel cores

model <- "ER" # Random graph model used to generate graphs.
# Possible values are:
# ER - Erdos-Renyi
# KR - d-regular graphs


interval <- c(0, 1)

ns <- c(20, 50, 100, 500, 1000, 10000)  

ps <- seq(0.1, 1, 0.2)

maxSim <- 1000

dists <- c("L1.ref", "L1.ref.cdf")

if (model == "ER") {
    ref <- p.semi.circle
    ref.cdf <- p.semi.circle.cdf
}

if (model == "KR") {
    ps <- c(3, 5, 10, -1)
    ref <- d.fixed.sd
    ref.cdf <- d.fixed.cdf
}



#p <- 0.4
#if (model == "SBM")
#    p <- "0.2_0.7"
#if (model == "PA")
#    p <- 1.4

executionName <- paste("parameter_estimator_theoretical_densities_", model, sep="")


set.seed(0)
seeds <- runif(maxSim, 0, 100000)

ncores <- 10 # number of cores for paralell execution

if (parallel)
    sfInit(parallel = TRUE , type ="SOCK", socketHosts = c(rep("localhost", ncores)))


# Exporting global variables ---------------------------------------------------

if (parallel) {
    sfExportAll() # exporta todas as variáveis
    sfLibrary(igraph) # exporta pacotes a serem utilizados
}


run <- function(sim) {
    system(paste("renice 10 -p", Sys.getpid()))
    result <- list()
    estimates <- matrix(NA, length(ns), length(ps)) 
    rownames(estimates) <- as.character(ns)
    colnames(estimates) <- as.character(ps)

    if (parallel) 
        set.seed(seeds[sim])

    for (dist in dists) {
        result[[dist]] <- estimates
    }

    for (p in ps) {
        for (n in ns) {
            parameter <- p
            if (p == -1) {
                parameter <- round(sqrt(n), 0)
            }
            if (model == "ER") {
                g <- erdos.renyi.game(n, parameter)
                g <- as.matrix(get.adjacency(g))
            }
            if (model == "KR") {
            	g <- k.regular.game(n, parameter)
                g <- as.matrix(get.adjacency(g))
            }            
            
            e <- as.numeric(eigen(g, only.values=TRUE, symmetric=TRUE)$values)[-1]
            if (model == "KR") {
               interval <- c(0, n)
            }
            else
                 e <- e/sqrt(n)
            cdf <- ecdf(e)

            f <- gaussianDensity(e, bandwidth="Silverman")

            fst <- gaussianDensity(e, bandwidth="Sturges")
            
            L1.ref <- function(x) {
                f2 <- list("x"=f$x, "y"=ref(f$x, x))
                return(L1(f, f2))
            }

            L1.ref.st <- function(x) {
                f2 <- list("x"=fst$x, "y"=ref(fst$x, x))
                return(L1(fst, f2))
            }

            L1.ref.cdf <- function(x) {
                cdf2 <- list("x"=cdf$x, "y"=ref.cdf(cdf$x, x))
                return(L1(cdf, cdf2))
            }
            
            matchFun <- function(name) {
                    return(match.fun(name))
            }
            for (dist in dists) {
                if (model == "ER") {
                    if (p < 0.5)
                        interval <- c(0, 0.5)
                    if (p > 0.5)
                        interval <- c(0.5, 1)
                    else
                        interval <- c(0, 1)
                }
                min <- optimize(matchFun(dist), interval=interval, maximum=FALSE,tol=1e-8)$minimum
                result[[dist]][as.character(n), as.character(p)] <- min
            }
        
        }
    }
    return(result)
}

# Paraleliza as N iterações em 40 núcleos (10 na máquina local, 12 na darwin e 18 na hulk3).
# Salva os resultados em "~/.sfCluster/restore/SAVE_DEFAULT_BAGULHO"

{if (parallel) {
    results <- sfClusterApplySR(1:maxSim, run, name=executionName, restore=F, perUpdate=1)

    sfStop()
}
else {
    results <- list()
    for (i in 1:maxSim) {
        results[[i]] <- run(i)
    }
}}

save(results, file=paste(directory, executionName, ".RData", sep=""))

################################################################################


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
