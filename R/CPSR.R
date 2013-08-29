##
## cPocockSimonRandomizer class
##
setClass("cPocockSimonRandomizer",
         representation(expt = "ClinicalExperiment",
                        seed = "integer",
                        stateTable = "matrix",
                        tr.assignments = "data.frame",
                        tr.ratios = "numeric",
                        max.counts = "numeric",
                        d.func = "function",
                        g.func = "function",
                        p.func = "function")
         contains = "PocockSimonRandomizer"
         )

setMethod("initialize", "cPocockSimonRandomizer",
          function (.Object, expt, seed, stateTable, tr.ratios, max.counts,
                    d.func=.default.d.func,
                    g.func=.default.g.func, p.func=.default.p.func) {
            ##print("I am called")
            if (missing(seed)) seed <- 12345
            if (missing(stateTable)) {
              m <- expt@number.of.treatments
              n <- sum(expt@number.of.factor.levels)
              stateTable <- matrix(0, nrow=m, ncol=n)
              rownames(stateTable) <- expt@treatment.names
              colnames(stateTable) <- unlist(sapply(1:expt@number.of.factors,
                                                    function(i) 
                                                    sapply(1:expt@number.of.factor.levels[i],
                                                           function(j) paste(expt@factor.names[i],
                                                                             expt@factor.level.names[[i]][j], sep=":"))))
            }
            if (missing(tr.ratios)) {
              m <- expt@number.of.treatments
              tr.ratios <- rep(1, m)
            }
            if (missing(max.counts)) {
              stop("Need max.counts for constrained Pocock Simon randomization.")
            }
            if (!length(max.counts) %in% c(1, expt@number.of.treatments)) {
              stop("max.counts must be vector of length 1 or number of treatments.")
            }
            if (length(max.counts) == 1) {
              max.counts <- rep(max.counts, expt@number.of.treatments)
            }

            .Object@expt <- expt
            .Object@seed <- as.integer(seed)
            .Object@stateTable <- stateTable
            .Object@tr.ratios <- tr.ratios / sum (tr.ratios)
            .Object@max.counts <- max.counts
            .Object@d.func <- d.func
            .Object@g.func <- g.func
            .Object@p.func <- p.func
            set.seed(seed)
            .Object
          })

setValidity("cPocockSimonRandomizer", .validPocockSimonRandomizer)

setMethod("randomize",
          signature(object="cPocockSimonRandomizer", subject.id="character", factor.values="character"),
          function (object, subject.id, factor.values) {
              if (missing(subject.id)) {
                  stop("Need subject id for randomization")
              }
              
              if (length(subject.id) > 1) {
                stop("Need a single subject id for randomization")
              }
              
              if (subject.id %in% rownames(object@tr.assignments)) {
                stop(paste("Subject ID", subject.id, "already randomized!"))
              }
                

              ##   factor.values <- as.integer(factor.values)
              if (missing(factor.values)) {
                  stop("Need factor values")
              }

              expt <- object@expt
              number.of.factors <- expt@number.of.factors
              if (length(factor.values) != number.of.factors) {
                stop("Not correct number of factors")
              }
              factor.level.names <- expt@factor.level.names
              factor.values.kosher <- sapply(1:number.of.factors,
                                             function(x) factor.values[x] %in% factor.level.names[[x]])
              if (! all(factor.values.kosher)) {
                stop("Incorrect factor values provided")
              }

              if (!"Counts" %in% expt@factor.names) {
                stop("'Counts' must be a factor for a constrained PS randomization")
              }
              
              number.of.treatments <- expt@number.of.treatments
              treatment.names <- expt@treatment.names
              factor.names <- expt@factor.names
              max.counts <- object@max.counts
              named.factors <- paste(factor.names, factor.values, sep=":")
              imbalances <- computeImbalances(object, factor.values)
              overallImbalance <- computeOverallImbalance(object, imbalances)
              treatmentCounts <- computeCounts(object, factor.values)
              ##print("Overall imbalance")
              ##print(overallImbalance)
              p.func <- object@p.func
              tr.ratios <- object@tr.ratios
              
              p.vec <- p.func(overallImbalance, treatmentCounts, max.counts)
              ##print("pvec")
              ##print(p.vec)
              tr.index <- sample(number.of.treatments, 1, prob=p.vec)
              tr.name <- treatment.names[tr.index]
              ##print(paste("Treatment is", tr.name))
              
              ## Update state
              state.matrix <- object@stateTable
              state.matrix[tr.name, named.factors] <- state.matrix[tr.name, named.factors] + 1
              stateTable(object) <- state.matrix
              
              ## Update assignment table
              current.assignment <- data.frame(as.list(factor.values), tr.name, stringsAsFactors=FALSE)
              rownames(current.assignment) <- subject.id
              colnames(current.assignment) <- c(factor.names, "Treatment")
              if (nrow(object@tr.assignments) == 0) {
                  assignments <- current.assignment
              } else {
                  assignments <- rbind(object@tr.assignments, current.assignment)
              }
              tr.assignments(object) <- assignments
              object
          })

##
## get the counts
##
if (!isGeneric("computeCounts")) {
    if (is.function("computeCounts")) {
        setGeneric("computeCounts", computeCounts)
    } else {
        setGeneric("computeCounts",
                   function(object, factor.values) standardGeneric("computeCounts"))
    }

}

setMethod("computeCounts",
          signature(object="cPocockSimonRandomizer"),
          function (object, factor.values) {
              ##   factor.values <- as.integer(factor.values)
              if (missing(factor.values)) {
                  stop("Need factor values")
              }
              expt <- object@expt
              number.of.factors <- expt@number.of.factors
              if (length(factor.values) != number.of.factors) {
                  stop("Not correct number of factors")
              }
              number.of.factor.levels <- expt@number.of.factor.levels
              factor.level.names <- expt@factor.level.names
              factor.values.kosher <- sapply(1:number.of.factors,
                                             function(x) factor.values[x] %in% factor.level.names[[x]])
              if (! all(factor.values.kosher)) {
                stop("Incorrect factor values provided")
              }
              
              treatment.names <- expt@treatment.names
              number.of.treatments <- expt@number.of.treatments
              d.func <- object@d.func
              factor.names <- expt@factor.names
              factor.level.names <- expt@factor.level.names
              treatment.names <- expt@treatment.names
              state.matrix <- object@stateTable
              tr.ratios <- object@tr.ratios
              named.factors <- paste(factor.names, factor.values, sep=":")
              f.mat <- state.matrix[, named.factors]
              ##print("our mat")
              ##print(f.mat)
              ##print(treatment.names)
              
              tab <- object@tr.assignments
              tab$Treatment <- factor(tab$Treatment, levels = treatment.names)
              tab$Counts <- factor(tab$Counts, levels = factor.level.names[[which(factor.names == "Counts")]])
              tab <- with(tab, table(Counts, Treatment)) * 
                as.numeric(factor.level.names[[which(factor.names == "Counts")]])
              counts <- colSums(tab) + as.numeric(factor.values[which(expt@factor.names == "Counts")])
              counts
          })
