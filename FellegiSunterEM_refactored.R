#start <- proc.time()
NOFIELD <- 3

# Counts
COUNTSTRUTH <- c(96995, 967, 954, 11, 946, 7, 15, 0, 0, 0, 0, 1, 0, 4, 2, 98)
OBSCOUNTS <- COUNTSTRUTH[1:2^NOFIELD] + COUNTSTRUTH[(2^NOFIELD + 1):(2^(NOFIELD + 1))]
COUNTSEXP <- c(0.999 * OBSCOUNTS, 0.001 * rev(OBSCOUNTS))

# Design matrix generators
gen_multinomial_design <- function(nf) {
  outer(0:(2^nf - 1), 0:(nf - 1), function(i, k) (i %/% 2^k) %% 2)
}

gen_double_poisson_design <- function(multi_mat) {
  nf <- ncol(multi_mat)
  n <- nrow(multi_mat)
  mat <- matrix(0, nrow = 2 * n, ncol = 2 * nf + 2)
  mat[1:n, 1] <- 1
  mat[1:n, 2:(nf + 1)] <- multi_mat
  mat[(n + 1):(2 * n), nf + 2] <- 1
  mat[(n + 1):(2 * n), (nf + 3):(2 * nf + 2)] <- multi_mat
  mat
}

MULTINOMIALDESIGNMAT <- gen_multinomial_design(NOFIELD)
DOUBLEPOISDESMAT <- gen_double_poisson_design(MULTINOMIALDESIGNMAT)

# Data and formula
POISREG <- data.frame(DOUBLEPOISDESMAT, COUNTSEXP)
predictors <- paste0("POISREG[,", 1:(2 * NOFIELD + 2), "]")
POISFORMULA <- as.formula(paste("POISREG[,", 2 * NOFIELD + 3, "] ~",
                                paste(predictors, collapse = "+"), "-1"))

# Initial regression
POISREGOUT <- glm(formula = POISFORMULA, family = poisson(),
                  maxit = 1000, epsilon = 1e-8)

# EM algorithm
for (i in 1:800) {
  for (h in 1:2^NOFIELD) {
    den <- POISREGOUT$fitted.values[h] + POISREGOUT$fitted.values[h + 2^NOFIELD]
    POISREG[h, 2 * NOFIELD + 3] <- (POISREGOUT$fitted.values[h] / den) * OBSCOUNTS[h]
    POISREG[h + 2^NOFIELD, 2 * NOFIELD + 3] <- (POISREGOUT$fitted.values[h + 2^NOFIELD] / den) * OBSCOUNTS[h]
  }
  POISREGOUT <- glm(formula = POISFORMULA, family = poisson(),
                    maxit = 1000, epsilon = 1e-8)
}

# Goodness-of-fit (Chi-squared)
fitvals <- POISREGOUT$fitted.values
chisq <- sum((OBSCOUNTS - (fitvals[1:2^NOFIELD] + fitvals[(2^NOFIELD + 1):(2^(NOFIELD + 1))]))^2 /
               (fitvals[1:2^NOFIELD] + fitvals[(2^NOFIELD + 1):(2^(NOFIELD + 1))]))

# Marginal and conditional probabilities
compute_conditionals <- function(cells, nf) {
  indices_u <- 1:(2^nf)
  indices_m <- (2^nf + 1):(2^(nf + 1))
  
  EUTOTAL <- sapply(1:nf, function(k) {
    sum(cells[indices_u][(indices_u %% 2^k == 0) | (indices_u %% 2^k > 2^(k - 1))])
  })
  
  EMTOTAL <- sapply(1:nf, function(k) {
    sum(cells[indices_m][(indices_m %% 2^k == 0) | (indices_m %% 2^k > 2^(k - 1))])
  })
  
  total_match <- sum(cells[indices_m])
  total_unmatch <- sum(cells[indices_u])
  PMATCH <- total_match / (total_match + total_unmatch)
  
  list(
    PMATCHNOM = PMATCH,
    PUCONDNOM = EUTOTAL / total_unmatch,
    PMCONDNOM = EMTOTAL / total_match
  )
}


probabilities <- compute_conditionals(cells = POISREGOUT$fitted.values, nf = NOFIELD)

# Output summaries
list(
  coefficients = POISREGOUT$coefficients,
  chisquare = chisq,
  probabilities = probabilities
)
