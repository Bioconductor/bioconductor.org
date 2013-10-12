Each section provides a function that supposedly works as expected, but quickly proves to misbehave. 
The exercise aims at first writing some dedicated testing functions that will identify the problems 
and then update the function so that it passes the specific tests. This practice is called unit testing 
and we use the `RUnit` package for this. For details on unit testing using `RUnit` 
see <http://bioconductor.org/developers/how-to/unitTesting-guidelines/>.

# Subsetting

## The buggy function

    ## Example
    isIn <- function(x, y) {
        sel <- match(x, y)
        y[sel]
    }
    
    ## Expected
    x <- sample(LETTERS, 5)
    isIn(x, LETTERS)
    
    ## Bug!
    isIn(c(x, "a"), LETTERS)

## A unit test and a solution

    ## Unit test:
    library("RUnit")
    test_isIn <- function() {
        x <- c("A", "B", "Z")
        checkIdentical(x, isIn(x, LETTERS))
        checkIdentical(x, isIn(c(x, "a"), LETTERS))
    
    }
    
    test_isIn()
    
    ## updated function
    isIn <- function(x, y) {
        sel <- x %in% y
        x[sel]
    }
    
    test_isIn()

# Character matching

## The buggy function

    ## Example
    isExactIn <- function(x, y)
        y[grep(x, y)]
    
    ## Expected
    isExactIn("a", letters)
    
    ## Bugs
    isExactIn("a", c("abc", letters))
    isExactIn(c("a", "z"), c("abc", letters))

## A unit test and a solution

    ## Unit test:
    library("RUnit")
    test_isExactIn <- function() {
        checkIdentical("a", isExactIn("a", letters))
        checkIdentical("a", isExactIn("a", c("abc", letters)))
        checkIdentical(c("a", "z"), isExactIn(c("a", "z"), c("abc", letters)))
    }
    
    test_isExactIn()
    
    ## updated function:
    isExactIn <- function(x, y)
        x[x %in% y]
    
    test_isExactIn()

# If conditions with length > 1

## The buggy function

    ## Example
    ifcond <- function(x, y) {
        if (x > y) {
            ans <- x*x - y*y
        } else {
            ans <- x*x + y*y
        } 
        ans
    }
    
    ## Expected
    do(3, 2)
    do(2, 2)
    do(1, 2)
    
    ## Bug!
    do(3:1, c(2, 2, 2))

## A unit test and a solution

    ## Unit test:
    library("RUnit")
    test_ifcond <- function() {
        checkIdentical(5, ifcond(3, 2))
        checkIdentical(8, ifcond(2, 2))
        checkIdentical(5, ifcond(1, 2))
        checkIdentical(c(5, 8, 5), ifcond(3:1, c(2, 2, 2)))
    }
    
    test_ifcond()
    
    ## updated function:
    ifcond <- function(x, y)
        ifelse(x > y, x*x - y*y, x*x + y*y)
    
    test_ifcond()

# Know your inputs

## The function

    ## Example
    distances <- function(point, pointVec) {
        x <- point[1]
        y <- point[2]
        xVec <- pointVec[,1]
        yVec <- pointVec[,2]
        dist <- sqrt((xVec - x)^2 + (yVec - y)^2)
        return(dist)
    }
    
    ## Expected
    x <- rnorm(5)
    y <- rnorm(5)
    
    m <- cbind(x, y)
    p <- m[1, ]
    
    distances(p, m)
    
    ## Bug!
    dd <- data.frame(x, y)
    q <- dd[1, ]
    
    distances(q, dd)

## A unit test and a solution

    ## Unit test:
    library("RUnit")
    test_distances <- function() {
        x <- y <- c(0, 1, 2)
        m <- cbind(x, y)
        p <- m[1, ]
        dd <- data.frame(x, y)
        q <- dd[1, ]
        expct <- c(0, sqrt(c(2, 8)))
        checkIdentical(expct, distances(p, m))
        checkIdentical(expct, distances(q, dd))
    }
    
    test_distances()
    
    ## updated function
    distances <- function(point, pointVec) {
        point <- as.numeric(point)
        x <- point[1]
        y <- point[2]
        xVec <- pointVec[,1]
        yVec <- pointVec[,2]
        dist <- sqrt((xVec - x)^2 + (yVec - y)^2)
        return(dist)
    }
    
    test_distances()

# Iterate on 0 length

## The buggy function

    ## Example
    sqrtabs <- function(x) {
        v <- abs(x)
        sapply(1:length(v), function(i) sqrt(v[i]))
    }
    
    ## Expected
    all(sqrtabs(c(-4, 0, 4)) == c(2, 0, 2))
    
    ## Bug!
    sqrtabs(numeric())

## A unit tests and a solution

    ## Unit test:
    library(RUnit)
    test_sqrtabs <- function() {
        checkIdentical(c(2, 0, 2), sqrtabs(c(-4, 0, 4)))
        checkIdentical(numeric(), sqrtabs(numeric()))
    }
    test_sqrtabs()
    
    ## updated function:
    sqrtabs <- function(x) {
      v <- abs(x)
      sapply(seq_along(v), function(i) sqrt(v[i]))
    }
    test_sqrtabs()                          # nope!
    
    sqrtabs <- function(x) {
      v <- abs(x)
      vapply(seq_along(v), function(i) sqrt(v[i]), 0)
    }
    test_sqrtabs()                          # yes!
