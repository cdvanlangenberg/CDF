list.packages <-
  c("permute", "plyr", "foreach", "doParallel", "binhf")

new.packages <-
  list.packages[!(list.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

# Load packages into session
sapply(list.packages, require, character.only = TRUE)

nsims <- 200
nperms <- 50

typecode <- 1

# Initialize vectors
ks.count.mean = 0
kp.count.mean = 0
cm.count.mean = 0
ks.count.med = 0
kp.count.med = 0
cm.count.med = 0
ks.count.95 = 0
kp.count.95 = 0
cm.count.95 = 0
ks.count.max = 0
kp.count.max = 0
cm.count.max = 0

set.seed(2311)

# start big loop
pb <-
  txtProgressBar(min = 0, max = nsims, style = 3) # this is not usefull if processing parallel

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

foreach(zz = 1:nsims) %do% {
  sapply(list.packages, require, character.only = TRUE)
  
  setTxtProgressBar(pb, zz)
  
  if (exists(c("combine", "ac.data", "as.data")))
    rm(combine, ac.data, as.data)
  
  # generate random deviates from an exponential distribution
  mean.as <- 0.2
  mean.ac <- 0.2
  nsamp <- 5
  nobs <- 30
  big <- nsamp * nobs
  
  pairs <- data.frame(
    pair = 1:nsamp ^ 2,
    grp1 = rep(seq_len(nsamp), each = nsamp),
    grp2 = rep(seq_len(nsamp) + nsamp, nsamp)
  )
  
  
  ac.data <- data.frame(x.sim = rexp(big, mean.ac))
  ac.data <- within(ac.data, {
    cdf.ac.sim <- pexp(x.sim, mean.ac)
    cdf <- cdf.ac.sim
    groupcode <- 1
  })
  
  as.data <- data.frame(x.sim = rexp(big, mean.as))
  as.data <- within(as.data, {
    cdf.as.sim <- pexp(x.sim, mean.as)
    cdf <- cdf.as.sim
    groupcode <- 2
  })
  
  
  combine <- rbind.fill(ac.data, as.data)
  combine <- within(combine, {
    typecode <- typecode
    sample <- gl(2 * nsamp, nobs)
  })
  
  combine <-
    combine[, c("typecode",
                "groupcode",
                "sample",
                "x.sim",
                "cdf",
                "cdf.as.sim",
                "cdf.ac.sim")] # this for just convenience of data visualization
  
  # pairwise summary
test.result <-
                foreach(ii = 1:nrow(pairs), .combine = rbind) %do% {
                  if (exists("sub.combine"))
                    rm(sub.combine)
                  # ii=1
                  sub.combine <-
                    combine[combine$sample %in% c(pairs[ii, 2], pairs[ii, 3]), ]
                  
                  sub.combine <- sub.combine[order(sub.combine$x.sim), ]
                  
                  sub.combine <- within(sub.combine, {
                    change <- ifelse(groupcode == shift(groupcode, 1), 0, 1)
                    change[1] = 0
                    cumsum <- cumsum(change)
                    #cdf.min <- ave(shift(cdf,1), cumsum, FUN=min )
                    cdf.ac.new <-
                      ifelse(is.na(cdf.ac.sim),
                             ave(shift(cdf, 1), cumsum, FUN = min),
                             cdf.ac.sim)
                    cdf.ac.new <-
                      ifelse(cumsum == 0 & is.na(cdf.ac.sim), 0, cdf.ac.new)
                    cdf.ac.new <-
                      ifelse(cumsum == cumsum[2 * nobs] &
                               is.na(cdf.ac.sim), 1, cdf.ac.new)
                    cdf.ac.new[2 * nobs] = 1
                    
                    cdf.as.new <-
                      ifelse(is.na(cdf.as.sim),
                             ave(shift(cdf, 1), cumsum, FUN = min),
                             cdf.as.sim)
                    cdf.as.new <-
                      ifelse(cumsum == 0 & is.na(cdf.as.sim), 0, cdf.as.new)
                    cdf.as.new <-
                      ifelse(cumsum == tail(cumsum, 1) &
                               is.na(cdf.as.sim),
                             1,
                             cdf.as.new)
                    cdf.as.new[2 * nobs] = 1
                    
                    cdf.diff <- cdf.as.new - cdf.ac.new
                    abs.cdf.diff <- abs(cdf.as.new - cdf.ac.new)
                    cdf.diffsq <- cdf.diff ^ 2
                    
                  })
                  
                  
                  test.result  <- data.frame(
                    ks = max(sub.combine$abs.cdf.diff),
                    kp = (
                      max(sub.combine$abs.cdf.diff) - min(sub.combine$abs.cdf.diff)
                    ),
                    cm = sum(sub.combine$cdf.diffsq)
                  )
                }
  # test.result  ## to print the results for a given iteration (zz)
  
  # overall summary statistics for above pairs
  
  test.stat <-
    apply(
      test.result,
      MARGIN = 2,
      FUN = function(x) {
        c(
          mean = mean(x),
          median = median(x),
          quantile(x, .95),
          max = max(x)
        )
      }
    )
  
  
  # summary.test.stat ## to print the summary of all pairs
  
  
  
  # begin permutation loop
  
  # nperms <- 10
  
  # register number of cores to speedup the process
  # note : if number of permutations are small parallel prosessing might take longer
  
  # cores <- detectCores()
  # cl <- makeCluster(cores)
  # registerDoParallel(cl)
  
  permute.results <-
    foreach(pp = 1:nperms, .combine = rbind) %dopar% {
      sapply(list.packages, require, character.only = TRUE)
      
      CTRL <-
        how(
          plots = Plots(combine$sample, type = "free"),
          within = Within(type = "none")
        )
      permuted <- shuffle(nrow(combine),  control = CTRL)
      
      p.combine <- rbind.fill(ac.data, as.data)
      p.combine <- within(p.combine, {
        typecode <- typecode
        sample <- combine$sample[permuted]
      })
      
      p.combine <-
        p.combine[order(p.combine$sample), c("typecode",
                                             "groupcode",
                                             "sample",
                                             "x.sim",
                                             "cdf",
                                             "cdf.as.sim",
                                             "cdf.ac.sim")] # this for just convenience of data visualization
      
      p.combine$groupcode <- rep(1:2, each = big)
      
      # pairwise summary for permutations
      
      p.result <-
        foreach(ii = 1:nrow(pairs), .combine = rbind) %do% {
          if (exists("p.sub.combine"))
            rm(p.sub.combine)
          #ii=2
          p.sub.combine <- p.combine[p.combine$sample %in% c(pairs[ii, 2], pairs[ii, 3]),]
          
          p.sub.combine <-
            p.sub.combine[order(p.sub.combine$x.sim), ]
          
          p.sub.combine <- within(p.sub.combine, {
            change <- ifelse(groupcode == shift(groupcode, 1), 0, 1)
            change[1] = 0
            cumsum <- cumsum(change)
            #cdf.min <- ave(shift(cdf,1), cumsum, FUN=min )
            cdf.ac.new <-
              ifelse(is.na(cdf.ac.sim),
                     ave(shift(cdf, 1), cumsum, FUN = min),
                     cdf.ac.sim)
            cdf.ac.new <-
              ifelse(cumsum == 0 & is.na(cdf.ac.sim), 0, cdf.ac.new)
            cdf.ac.new <-
              ifelse(cumsum == cumsum[2 * nobs] &
                       is.na(cdf.ac.sim), 1, cdf.ac.new)
            cdf.ac.new[2 * nobs] = 1
            
            cdf.as.new <-
              ifelse(is.na(cdf.as.sim),
                     ave(shift(cdf, 1), cumsum, FUN = min),
                     cdf.as.sim)
            cdf.as.new <-
              ifelse(cumsum == 0 & is.na(cdf.as.sim), 0, cdf.as.new)
            cdf.as.new <-
              ifelse(cumsum == cumsum[2 * nobs] &
                       is.na(cdf.as.sim), 1, cdf.as.new)
            cdf.as.new[2 * nobs] = 1
            
            cdf.diff <- cdf.as.new - cdf.ac.new
            abs.cdf.diff <- abs(cdf.diff)
            cdf.diffsq <- cdf.diff ^ 2
            
          })
          
          
          p.result  <- data.frame(
            ks.p = max(p.sub.combine$abs.cdf.diff),
            kp.p = (
              max(p.sub.combine$abs.cdf.diff) - min(p.sub.combine$abs.cdf.diff)
            ),
            cm.p = sum(p.sub.combine$cdf.diffsq)
          )
        }
      
      permute.results <-   apply(
        p.result,
        MARGIN = 2,
        FUN = function(x) {
          c(
            mean = mean(x),
            median = median(x),
            quantile(x, .95),
            max = max(x)
          )
        }
      )
      
    } # end of permutation loop
  
  # stopCluster(cl)
  
  test.stat1 <-
    apply(
      test.stat,
      2,
      FUN = function(x)
        replicate(nperms, x)
    )
  permute.results1 <- permute.results >= test.stat1
  pval.mean <-
    (apply(X = permute.results1[row.names(permute.results1) == "mean", ], MARGIN = 2, FUN = sum) +
       1) / (nperms + 1)
  pval.med <-
    (apply(X = permute.results1[row.names(permute.results1) == "median", ], MARGIN = 2, FUN = sum) +
       1) / (nperms + 1)
  pval.95 <-
    (apply(X = permute.results1[row.names(permute.results1) == "95%", ], MARGIN = 2, FUN = sum) +
       1) / (nperms + 1)
  pval.max <-
    (apply(X = permute.results1[row.names(permute.results1) == "max", ], MARGIN = 2, FUN = sum) +
       1) / (nperms + 1)
  
  # print(pval.mean) ## to check calculation
  
  if (pval.mean[1] <= 0.05)
    ks.count.mean = ks.count.mean + 1
  if (pval.med[1] <= 0.05)
    ks.count.med = ks.count.med + 1
  if (pval.95[1] <= 0.05)
    ks.count.95 = ks.count.95 + 1
  if (pval.max[3] <= 0.05)
    ks.count.max = ks.count.max + 1
  
  if (pval.mean[2] <= 0.05)
    kp.count.mean = kp.count.mean + 1
  if (pval.med[2] <= 0.05)
    kp.count.med = kp.count.med + 1
  if (pval.95[2] <= 0.05)
    kp.count.95 = kp.count.95 + 1
  if (pval.max[2] <= 0.05)
    kp.count.max = kp.count.max + 1
  
  if (pval.mean[3] <= 0.05)
    cm.count.mean = cm.count.mean + 1
  if (pval.med[3] <= 0.05)
    cm.count.med = cm.count.med + 1
  if (pval.95[3] <= 0.05)
    cm.count.95 = cm.count.95 + 1
  if (pval.max[3] <= 0.05)
    cm.count.max = cm.count.max + 1
  
  #rm(combine, ac.data, as.data, sub.combine, test.result)
  #gc()
  
}

stopCluster(cl)
close(pb)

ks.rej.mean <- ks.count.mean / nsims
ks.rej.med <- ks.count.med / nsims
ks.rej.95 <- ks.count.95 / nsims
ks.rej.max <- ks.count.max / nsims
kp.rej.mean <- kp.count.mean / nsims
kp.rej.med <- kp.count.med / nsims
kp.rej.95 <- kp.count.95 / nsims
kp.rej.max <- kp.count.max / nsims
cm.rej.mean <- cm.count.mean / nsims
cm.rej.med <- cm.count.med / nsims
cm.rej.95 <- cm.count.95 / nsims
cm.rej.max <- cm.count.max / nsims

cat(
  " nsims :",
  nsims,
  "\n",
  "nperms :",
  nperms,
  "\n",
  "nsamp :",
  nsamp,
  "\n",
  "nobs :",
  nobs,
  "\n",
  "mean.ac :",
  mean.ac,
  "\n",
  "mean.as :",
  mean.as,
  "\n",
  "ks.rej.mean :",
  ks.rej.mean,
  "\n",
  "ks.rej.med :",
  ks.rej.med,
  "\n",
  "ks.rej.95 :",
  ks.rej.95,
  "\n",
  "ks.rej.max :",
  ks.rej.max,
  "\n",
  "kp.rej.mean :",
  kp.rej.mean,
  "\n",
  "kp.rej.med :",
  kp.rej.med,
  "\n",
  "kp.rej.95 :",
  kp.rej.95,
  "\n",
  "kp.rej.max :",
  kp.rej.max,
  "\n",
  "cm.rej.mean :",
  cm.rej.mean,
  "\n",
  "cm.rej.med :",
  cm.rej.med,
  "\n",
  "cm.rej.95 :",
  cm.rej.95,
  "\n",
  "cm.rej.max :",
  cm.rej.max,
  "\n"
)
