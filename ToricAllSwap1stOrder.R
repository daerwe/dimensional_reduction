NOFIELD <- 3 #This program only works for 3 comparison fields
NROWS <- 2^(NOFIELD - 1) * NOFIELD *2^NOFIELD
NCOLS <- NOFIELD - 1
TORICDESIGNMATALLSWAP1ST <- matrix(nrow = NROWS, ncol = NCOLS)

# Generate binary design matrix
generate_design_matrix <- function(n) {
  m <- 2^n
  binary_matrix <- matrix(0, nrow = m - 1, ncol = n)
  for (i in 2:m) {
    binary_matrix[i - 1, ] <- as.integer(intToBits(i - 1))[1:n]
  }
  rbind(rep(0, n), binary_matrix)
}

# Swapping function
apply_swap <- function(mat, design, cols) {
  swapped <- mat
  for (col in cols) {
    idx1 <- which(design[, col] == 1)
    idx0 <- which(design[, col] == 0)
    swapped[idx1, ] <- mat[idx0, ]
    swapped[idx0, ] <- mat[idx1, ]
  }
  swapped
}

# Loop over each survivor field
for (ISURV in 1:NOFIELD) {
  survivors <- setdiff(1:NOFIELD, ISURV)
  nominal <- generate_design_matrix(NOFIELD)
  surv_mat <- nominal[, survivors, drop = FALSE] * (NOFIELD - length(survivors) + 1)
  phantom_mat <- nominal[2^NOFIELD:1, ISURV, drop = FALSE] - 1
  
  # Add phantom columns
  for (i in seq_along(survivors)) {
    surv_mat[, i] <- surv_mat[, i] + phantom_mat
  }
  
  start_idx <- function(swap_idx) {
    ((ISURV - 1) * 2^(NOFIELD - 1) * 2^NOFIELD) + ((swap_idx - 1) * 2^NOFIELD) + 1
  }
  
  TORICDESIGNMATALLSWAP1ST[start_idx(1):(start_idx(1) + 2^NOFIELD - 1), ] <- surv_mat
  TORICDESIGNMATALLSWAP1ST[start_idx(2):(start_idx(2) + 2^NOFIELD - 1), ] <- apply_swap(surv_mat, nominal, 1)
  TORICDESIGNMATALLSWAP1ST[start_idx(3):(start_idx(3) + 2^NOFIELD - 1), ] <- apply_swap(surv_mat, nominal, 2)
  TORICDESIGNMATALLSWAP1ST[start_idx(4):(start_idx(4) + 2^NOFIELD - 1), ] <- apply_swap(apply_swap(surv_mat, nominal, 1), nominal, 2)
}

