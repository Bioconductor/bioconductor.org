Each section provides a function that supposedly works as expected,
but quickly proves to misbehave.  The exercise aims at first writing
some dedicated testing functions that will identify the problems and
then update the function so that it passes the specific tests. This
practice is called unit testing and we use the `RUnit` package for
this. See the
[Unit Testing How-To](/developers/how-to/unitTesting-guidelines/)
guide for details on unit testing using `RUnit`.

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

# Know your inputs

## The function

    ## Example
    distances <- function(point, pointVec) {
        x <- point[1]
        y <- point[2]
        xVec <- pointVec[,1]
        yVec <- pointVec[,2]
        sqrt((xVec - x)^2 + (yVec - y)^2)
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

