    calcExpectedBAF <- function(exp.intA, exp.intB, tumor.purity) {
        exp.BAF <- (1 - tumor.purity + (tumor.purity * exp.intB)) / (2 - 2 * tumor.purity + tumor.purity * (exp.intA + exp.intB))
        return(exp.BAF)
    }

    simulateBAFs <- function(nrSNPs = 100, purity = 0.1, cpnA = 2, cpnB = 1, depth = 200) {
        nBottom <- rbinom(n = 1, size = nrSNPs, prob = 0.5)
        nTop    <- nrSNPs - nBottom
        expectedBAFmin <- calcExpectedBAF(cpnA, cpnB, purity)
        expectedBAFmaj <- 1 - expectedBAFmin

        sim        <- colSums(replicate(nBottom, rbinom(n = depth, size = 1, prob = expectedBAFmin)))
        # tmp.bottom <- data.frame(varCount = sim, refCount = depth - sim, BAF = sim / depth, class = "bottom")
        tmp.bottom <- data.frame(nrSNPs = nrSNPs, purity = purity, depth = depth, cpnA = cpnA, cpnB = cpnB, BAF = sim / depth, class = "bottom")    
        
        sim        <- colSums(replicate(nTop, rbinom(n = depth, size = 1, prob = expectedBAFmaj)))
        # tmp.top    <- data.frame(varCount = sim, refCount = depth - sim, BAF = sim / depth, class = "top")
        tmp.top    <- data.frame(nrSNPs = nrSNPs, purity = purity, depth = depth, cpnA = cpnA, cpnB = cpnB, BAF = sim / depth, class = "top")    

        return(rbind(tmp.top, tmp.bottom))
    }