######### SECTION Functions
## ANCHOR Polynomial
polynomialSum <- function(...) {
  x = list(...)
  n = length(x)
  for (p in x) {
     n = max(n, length(p))
  }
  sum = rep(0, n)
  for (p in x) {
     sum <- sum + c(p, rep(0, n-length(p)))
  }
  return(sum)
}

polynomialProduct <- function(...) {
  l = list(...)
  if (length(l)==2) {
    x = l[[1]]
    y = l[[2]]
    product <- rep(0, length(x) + length(y) - 1)
    # length(product) <- length(x) + length(y)
    for (i in 1:length(x)) {
      for (j in 1:length(y)) {
          product[i+j-1] <- product[i+j-1] + x[i] * y[j]
      }
    }
    return(product)
  } 
  else if (length(l)==1) {
    return(l[[1]])
  }
  else {
    return(polynomialProduct(polynomialProduct(l[[1]],l[[2]]),l[[3:length(l)]]))
  }
}

polynomialEquation <- function (coefficients,x="x") {
  eq <- ""
  for (i in 1:length(coefficients)) {
    if(coefficients[i] != 0) {
      if (i==1) {
        eq <- paste(eq, coefficients[1], sep = "")
      }
      else if (i==2) {
        eq <- paste(eq,ifelse(eq!="","+",""),coefficients[i],"*",x, sep = "")
      } else {
        eq <- paste(eq,ifelse(eq!="","+",""),coefficients[i],"*",x,"^",i-1, sep = "")
      }
    }
  }
  return(ifelse(eq=="","0",eq))
}

## ANCHOR Quadratic
quadraticInterpolate <- function(z,x,y) {
  return(
  y[1] + (y[2]-y[1])/(x[2]-x[1])*(z-x[1])+ ((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3])*(z-x[1])*(z-x[2])
  )
}

quadraticCoefficients <- function(x,y) {
  return(
    polynomialSum(c(y[1]), (y[2]-y[1])/(x[2]-x[1])*c(-x[1],1), ((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3])*polynomialProduct(c(-x[1],1),c(-x[2],1)))
  )
}

quadraticEquation <- function(x,y) {
  return(
    paste(as.character(y[1]),"+",(y[2]-y[1])/(x[2]-x[1]),"(x-",x[1],")","+",((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3]),"(x-",x[1],")","(x-",x[2],")", sep="")
  )
}

## ANCHOR Patching
patch <- function(x,y,p) {
  return(
    (1 - 3*p^2 + 2*p^3) * (x - y) + y
  )
}

patchPolynomial <- function(x,y,p) {
  return(
    polynomialSum(polynomialProduct(polynomialSum(c(1),-3*polynomialProduct(p,p),2*polynomialProduct(p,p,p)),polynomialSum(x,-y)), y)
  )
}

patchEquation <- function(x,y,p) {
  return(
    paste("(1-3(",p,")^2+2(",p,")^3)*(",x,")+(3(",p,")^2-2(",p,")^3)*(",y,"))", sep = "")
  )
}

## ANCHOR Interpolate
interpolate <- function(z,x,y) {
  f <- rep(NA, length(z))
  j <- 0
  for (i in 1:length(x)) {
    k <- j + 1
    while(z[j+1] <= x[i]) {
      j <- j + 1
      if(j == length(z)) {
        break
      }
    }

    if(j >= k) {
      if (i == 1 & z[j] == x[1]) {
       f[j] <- x[1]
      }
      else if (i == 2) {
        f[k:j] <- quadraticInterpolate(z[k:j],x[1:3],y[1:3])
      }
      else if (i==length(x)) {
        f[k:j] <- quadraticInterpolate(z[k:j],x[(i-2):i],y[(i-2):i])
      }
      else {
        f[k:j] <- patch(
            quadraticInterpolate(z[k:j],x[(i-2):i],y[(i-2):i]), 
            quadraticInterpolate(z[k:j],x[(i-1):(i+1)],y[(i-1):(i+1)]), 
            (z[k:j]-x[i-1])/(x[i]-x[i-1])
          )
      }
    }

    if(j == length(z)) {
      break
    }
  }

  return(f)
}

piecewiseRange <- function(min, max, x="x") {
  return(paste("(",x,">",min," & ",x,"<",max,")", sep = ""))
}

interpolationEquation <- function(x,y) {
  eq <- paste(piecewiseRange(x[1],x[2]),"*(",quadraticEquation(x[1:3],y[1:3]),")", sep = "")
  for (i in 2:(length(x)-1)) {
    if (i==length(x)-1) {
      eq <- paste(eq," + ",piecewiseRange(x[i],x[i+1]),"*(",quadraticEquation(x[(i-1):(i+1)],y[(i-1):(i+1)]),")", sep = "")
    }
    else {
      eq <- paste(eq," + ",piecewiseRange(x[i],x[i+1]),"*(",patchEquation(quadraticEquation(x[(i-1):(i+1)],y[(i-1):(i+1)]),quadraticEquation(x[i:(i+2)],y[i:(i+2)]),paste("(x-",x[i],")/",x[i+1]-x[i], sep = "")),")", sep = "")
    }
  }
  return(eq)
}

interpolationEquationExpanded <- function(x,y) {
  eq <- paste(piecewiseRange(x[1],x[2]),"*(",polynomialEquation(quadraticCoefficients(x[1:3],y[1:3])),")", sep = "")
  for (i in 2:(length(x)-1)) {
    if (i==length(x)-1) {
      eq <- paste(eq," + ",piecewiseRange(x[i],x[i+1]),"*(",polynomialEquation(quadraticCoefficients(x[(i-1):(i+1)],y[(i-1):(i+1)])),")", sep = "")
    }
    else {
      eq <- paste(eq," + ",piecewiseRange(x[i],x[i+1]),"*(",polynomialEquation(patchPolynomial(quadraticCoefficients(x[(i-1):(i+1)],y[(i-1):(i+1)]),quadraticCoefficients(x[i:(i+2)],y[i:(i+2)]),c(-x[i],1)/(x[i+1]-x[i]))),")", sep = "")
    }
  }
  return(eq)
}

graphInterpolation <- function(x,y,xx=seq(min(x),max(x),0.05)) {
  plot(xx, interpolate(xx,x,y),type="l",xlab="x", ylab="y")
  points(x,y,col="red")
}
######### !SECTION

######### SECTION Main
x <- c(1,2,3,4,5)
y <- c(1,2^2,3^2,2^2,1)

print(interpolationEquation(x,y))
print(interpolationEquationExpanded(x,y))
graphInterpolation(x,y)