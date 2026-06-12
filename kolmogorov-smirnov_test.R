# to test continuity/discontinuity between time period in the same region
# M is the f2-pairwise matrix, meta is a metafile with 

ks_continuity <- function(M, meta, region, t1, t2, B = 10000) {
  idx <- which(meta$region == region & meta$timeslice %in% c(t1, t2))
  ts  <- meta$timeslice[idx]; sub <- M[idx, idx]
  i1 <- which(ts == t1); i2 <- which(ts == t2)
  A <- c(get_d(sub,i1,i1), get_d(sub,i2,i2))   
  B_ <- get_d(sub, i1, i2)                      
  obsD <- as.numeric(ks.test(B_, A)$statistic)
  shift <- median(B_) - median(A)              
  perm <- replicate(B, {                        
    p <- sample(ts); a <- which(p==t1); b <- which(p==t2)
    as.numeric(ks.test(get_d(sub,a,b),
                       c(get_d(sub,a,a), get_d(sub,b,b)))$statistic)
  })
  data.frame(region, transizione = paste(t1,"->",t2),
             KS_D = obsD, shift_mediana = shift,
             p = (1 + sum(perm >= obsD)) / (B + 1))
}

ks_continuity(m, meta2, "West-North", "17-14_kaBP", "14-11_kaBP")
