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

NOFIELD <- length(m)  # Number of comparison fields 

# Aggregated counts over matching status
OBSCOUNTS <- COUNTSTRUTH[1:2^NOFIELD] + COUNTSTRUTH[(2^NOFIELD + 1):(2^(NOFIELD + 1))]
