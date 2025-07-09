generate_counts <- function(N, p, m, u) {
  n <- length(m)  # Number of comparison fields
  
  # Generate match status
  match_status <- rbinom(N, 1, p)
  
  # Build a matrix of probabilities row-by-row
  probs <- t(sapply(match_status, function(ms) if (ms == 1) m else u))
  
  # Generate binary values
  match_vectors <- as.data.frame(matrix(rbinom(N * n, size = 1, prob = as.vector(probs)), nrow = N, byrow = FALSE))
  
  match_vectors <- cbind(match_vectors, match_status)  # Append match/non-match column
  
  # Compute indices in standard binary order
  vector = 2^seq(0,n)
  indices <- rowSums(mapply(`*`, match_vectors, vector))
  # Count occurrences
  counts <- tabulate(indices + 1, nbins = 2^(n + 1))
  
  return(counts)
}

# Example usage
set.seed(42)
N <- 100000  # Number of observations
p <- 0.001    # Probability of a match
m <- c(0.99, 0.95, 0.99)  # Match probabilities for each field
u <- c(0.01, 0.01, 0.01)  # Non-match probabilities for each field

#Generate random true counts
#COUNTSTRUTH <- generate_counts(N, p, m, u)

#Choose counts
COUNTSTRUTH = c(96995, 967, 954, 11, 946, 7, 15, 0, 0, 0, 0, 1, 0, 4, 2, 98)
print(COUNTSTRUTH)

# Conditional-independence model from Fellegi-Sunter paper
# Supports both supervised & unsupervised learning approaches

NOFIELD <- 3  # Number of comparison fields (modifiable)

# Aggregated counts over matching status
OBSCOUNTS <- COUNTSTRUTH[1:2^NOFIELD] + COUNTSTRUTH[(2^NOFIELD + 1):(2^(NOFIELD + 1))]

# Generate Multinomial Design Matrix
MULTINOMIALDESIGNMAT <- matrix(0, nrow = 2^NOFIELD, ncol = NOFIELD)
for (i in 1:2^NOFIELD) {
  for (k in 1:NOFIELD) {
    MULTINOMIALDESIGNMAT[i, k] <- ifelse(((i - 1) %/% 2^(k-1)) %% 2 > 0, 1, 0)
  }
}

# Create Poisson Design Matrix
CONTRASTMAT <- 2 * MULTINOMIALDESIGNMAT - 1
DOUBLEPOISDESMAT <- matrix(0, nrow = 2 * 2^NOFIELD, ncol = 2 * NOFIELD + 2)
DOUBLEPOISDESMAT[1:2^NOFIELD, 1] <- 1
DOUBLEPOISDESMAT[1:2^NOFIELD, 2:(NOFIELD + 1)] <- CONTRASTMAT
DOUBLEPOISDESMAT[(2^NOFIELD + 1):(2*2^NOFIELD), NOFIELD + 2] <- 1
DOUBLEPOISDESMAT[(2^NOFIELD + 1):(2*2^NOFIELD), (NOFIELD + 3):(2*NOFIELD + 2)] <- CONTRASTMAT

# Poisson regression formula
ARGUMENTS <- paste0("POISREG[, ", 1:(2*NOFIELD + 2), "]")
POISFORMULA <- as.formula(paste("POISREG[, 2*NOFIELD + 3] ~", paste(ARGUMENTS, collapse = "+"), "-1"))

# Perform Poisson regression on truth
POISREG <- as.data.frame(cbind(DOUBLEPOISDESMAT, COUNTSTRUTH))
POISREGOUT <- glm(formula = POISFORMULA, family = poisson, maxit = 1000, epsilon = 1e-6)

# Compute marginal & conditional probabilities
compute_probs <- function(coeffs, range) {
  exp(2 * coeffs[range]) / (1 + exp(2 * coeffs[range]))
}
mprobtruth <- compute_probs(POISREGOUT$coefficients, (NOFIELD + 3):(2*NOFIELD + 2))
uprobtruth <- compute_probs(POISREGOUT$coefficients, 2:(NOFIELD + 1))

# Neyman-Pearson goodness-of-fit chi-square
neymantruth <- sum((COUNTSTRUTH - POISREGOUT$fitted.values)^2 / POISREGOUT$fitted.values)

# EM Algorithm for mixture likelihood estimation
COUNTSEXP <- c(.99 * OBSCOUNTS, .01 * rev(OBSCOUNTS))
POISREG <- as.data.frame(cbind(DOUBLEPOISDESMAT, COUNTSEXP))
POISREGOUT <- glm(formula = POISFORMULA, family = poisson, maxit = 1000, epsilon = 1e-6)

for (i in 1:2000) {
  for (h in 1:2^NOFIELD) {
    total_fitted <- POISREGOUT$fitted.values[h] + POISREGOUT$fitted.values[h + 2^NOFIELD]
    POISREG[h, 2*NOFIELD + 3] <- (POISREGOUT$fitted.values[h] / total_fitted) * OBSCOUNTS[h]
    POISREG[h + 2^NOFIELD, 2*NOFIELD + 3] <- (POISREGOUT$fitted.values[h + 2^NOFIELD] / total_fitted) * OBSCOUNTS[h]
  }
  POISREGOUT <- glm(formula = POISFORMULA, family = poisson, maxit = 1000, epsilon = 1e-6)
}

# Compute probabilities from EM algorithm
mprobem <- compute_probs(POISREGOUT$coefficients, (NOFIELD + 3):(2*NOFIELD + 2))
uprobem <- compute_probs(POISREGOUT$coefficients, 2:(NOFIELD + 1))

# Neyman-Pearson goodness-of-fit for EM estimates
neymanem <- sum((OBSCOUNTS - (POISREGOUT$fitted.values[1:2^NOFIELD] + 
                                POISREGOUT$fitted.values[(2^NOFIELD + 1):(2^(NOFIELD + 1))]))^2 / 
                  (POISREGOUT$fitted.values[1:2^NOFIELD] + POISREGOUT$fitted.values[(2^NOFIELD + 1):(2^(NOFIELD + 1))]))

# Compute marginal & conditional probabilities
ETOTAL <- c(sum(POISREGOUT$fitted.values[1:2^NOFIELD]), sum(POISREGOUT$fitted.values[(2^NOFIELD + 1):(2^(NOFIELD + 1))]))
EMTOTAL <- sapply(1:NOFIELD, function(k) sum(POISREGOUT$fitted.values[((2^NOFIELD + 1):2^(NOFIELD + 1))[(i %% 2^k == 0) | (i %% 2^k > 2^(k-1) )]]))
EUTOTAL <- sapply(1:NOFIELD, function(k) sum(POISREGOUT$fitted.values[(1:2^NOFIELD)[(i %% 2^k == 0) | (i %% 2^k > 2^(k-1) )]]))
PUCOND <- EUTOTAL / ETOTAL[1]
PMCOND <- EMTOTAL / ETOTAL[2]

PMATCHNOM <- ETOTAL[2] / sum(ETOTAL)
PUCONDNOM <- PUCOND
PMCONDNOM <- PMCOND
