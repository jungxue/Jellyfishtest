

x <-c(0,3,6)
y <-c(0,3,0)

getTriangleArea <- function(x, y) {
  # Write your code here
  base <- range(y)[2]-range(y)[1]
  height <- range(x)[2]-range(x)[1]
  Area <- (base*height)/2
  return(Area)
}
getTriangleArea(x,y)

getTriangleArea <- function(x, y) {
  matrixX <- matrix(c(x,y,rep(1,3)),ncol=3)
  Area <- abs(0.5*det(matrixX))
  return(Area)
}
getTriangleArea(x,y)

crossdist(x,y)
