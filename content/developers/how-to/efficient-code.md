# Robust and Efficient Code

_R_ can be a robust, fast and efficient programming language, but some
coding practices can be very unfortunate. Here are some suggestions.

## Guiding Principles

1. The primary principle is to make sure your code, and attempts to
   make it efficient, are correct first. Use `identical()` or
   `all.equal()` to ensure correctness.

2. Write robust code. Avoid efficiencies that do not easily handle
   edge cases such as 0 length or `NA` values.

3. Know when to stop trying to be efficient. If the code takes a
   fraction of a second to evaluate, there is no sense in trying for
   further improvement. Use `system.time()` or a package like
   [microbenchmark][] to quantify performance gains.

## Common Solutions

1. Write `seq_len(n)` or `seq_along(x)` rather than `1:n` or
   `1:length(x)`. This protects against the case when `n` or
   `length(x)` is 0 (which often occurs as an unexpected 'edge case'
   in real code) or negative.

2. Preallocate-and-fill (usually via `lapply()` or `vapply()`) rather
   than copy-and-append. If creating a vector or list of results, use
   `lapply()` (to create a list) or `vapply()` (to create a vector)
   rather than a `for` loop. For instance,

        n <- 10000
        x <- vapply(seq_len(n), function(i) {
            ## ...
        }, integer(1))

    manages the memory allocation of `x` and compactly represents the
    transformation being performed. A `for` loop might be appropriate
    if the iteration has side effects (e.g., displaying a plot) or
    where calculation of one value depends on a previous value. When
    creating a vector in a `for` loop, always pre-allocate the result

        x <- integer(n)
        if (n > 0) x[1] <- 0
        for (i in seq_len(n - 1)) {
            ## x[i + 1] <- ...
        }

    **Never** adopt a strategy of 'copy-and-append'

        not_this <- function(n) {
            x <- integer()
            for (i in seq_len(n))
                x[i] = i
            x
        }

    This pattern copies the current value of `x` each time through the
    loop, making `n^2 / 2` total copies, and scale very poorly even
    for trivial computations:

        > system.time(not_this(1000)
           user  system elapsed 
          0.004   0.000   0.004 
        > system.time(not_this(10000))
           user  system elapsed 
          0.169   0.000   0.168 
        > system.time(not_this(100000))
           user  system elapsed 
         22.827   1.120  23.936 

3. Vectorize, rather than iterate. A single call `y <- sqrt(x)` with a
   vector `x` of length `n` is an example of a vectorized function. A
   call such as `y <- sapply(x, sqrt)` or a `for` loop
   `for (i in seq_along(x)) y[i] <- sqrt(x[i])` is an iterative version
   of the same call, and should be avoided. Often, iterative
   calculations that can be vectorized are "hidden" in the body of a
   `for` loop

        for (i in seq_len(n)) {
            ## ...
            tmp <- foo(x[i])
            y <- tmp + ## ...
        }

    and can be 'hoisted' out of the loop

        tmp <- foo(x)
        for (i in seq_len(n)) {
            ## ...
            y <- tmp[i] + ##
        }

    Often this principle can be applied repeatedly, and an iterative
    `for` loop becomes a few lines of vectorized function calls.

[microbenchmark]: https://cran.r-project.org/web/packages/microbenchmark
