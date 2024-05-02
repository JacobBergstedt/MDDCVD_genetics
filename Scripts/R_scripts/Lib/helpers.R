
cor_mat_from_long_format <- function(dat, values) {
  
  traits <- unique(c(dat$Trait1, dat$Trait2))
  n_traits <- length(traits)
  mat <- matrix(NA, ncol = n_traits, nrow = n_traits)
  colnames(mat) <- rownames(mat) <- traits
  
  
  
  for (col in colnames(mat)) {
    for (row in rownames(mat)) {
      
      print(paste0(col, "_", row))
      if (col != row) {
        
        is_col <- dat$Trait1 == col | dat$Trait2 == col
        is_row <- dat$Trait1 == row | dat$Trait2 == row
        this <- is_col & is_row
        gcov <- dat[[values]][this]
        
      } else {
        
        gcov <- 1
        
      }
      
      mat[col, row] <- gcov
      
    }
  }
  
  mat
  
}

