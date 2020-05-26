orientate <- function(bx, by) {
  tmp <- (bx < 0 & by < 0) | (bx < 0 & by > 0)
  bx <- abs(bx)
  by[tmp] <- -by[tmp]
  return(list(bx = bx, by = by))
}
