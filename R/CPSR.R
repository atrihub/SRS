setMethod("randomize",
          signature(object="PocockSimonRandomizer", subject.id="character", factor.values="character"),
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
              
              number.of.treatments <- expt@number.of.treatments
              treatment.names <- expt@treatment.names
              factor.names <- expt@factor.names
              named.factors <- paste(factor.names, factor.values, sep=":")
              imbalances <- computeImbalances(object, factor.values)
              overallImbalance <- computeOverallImbalance(object, imbalances)
              ##print("Overall imbalance")
              ##print(overallImbalance)
              p.func <- object@p.func
              tr.ratios <- object@tr.ratios
              
              p.vec <- p.func(overallImbalance)
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