# Spearman's rank correlation permuting for individuals instead of pairwise comparisons
# M is the matrix of pairwise f2 values, D a matrix of distances between samples (can be geographical or spatial distance or other), in the same format as the f2 matrix.
# n_perm is the number of permutation desired


node_perm_cor <- function(M, D, n_perm = 10000, method = "spearman") {
  iu <- lower.tri(M)
  y  <- M[iu]
  x  <- D[iu]
  obs <- cor(x, y, method = method)
  n <- nrow(M)
  perm <- replicate(n_perm, {
    p <- sample.int(n)            
    Dp <- D[p, p]                  
    cor(Dp[iu], y, method = method)
  })
  list(rho = obs,
       p_two_sided = (1 + sum(abs(perm) >= abs(obs))) / (n_perm + 1),
       n_pairs = sum(iu))
}

# use 
node_perm_cor(M, D, n_perm = 10000)
