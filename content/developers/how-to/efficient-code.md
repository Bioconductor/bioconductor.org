# Robust and Efficient Code

_R_ can be a robust, fast and efficient programming language, but some
coding practices can be very unfortunate. Here are some suggestions.

## Guiding Principles

1. The primary principle is to make sure your code is **correct**. Use
   `identical()` or `all.equal()` to ensure correctness, and
   [unit tests][] to ensure consistent results across code revisions.

2. Write robust code. Avoid efficiencies that do not easily handle
   edge cases such as 0 length or `NA` values.

3. Know when to stop trying to be efficient. If the code takes a
   fraction of a second to evaluate, there is no sense in trying for
   further improvement. Use `system.time()` or a package like
   [microbenchmark][] to quantify performance gains.

## Common Advice

### Vectorize

Vectorize, rather than iterate (`for` loops, `lapply()`, `apply()` are
common iteration idioms in _R_). A single call `y <- sqrt(x)` with a
vector `x` of length `n` is an example of a vectorized function. A
call such as `y <- sapply(x, sqrt)` or a `for` loop `for (i in
seq_along(x)) y[i] <- sqrt(x[i])` is an iterative version of the same
call, and should be avoided. Often, iterative calculations that can be
vectorized are "hidden" in the body of a `for` loop

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

### 'Pre-allocate and fill' if iterations are necessary

Preallocate-and-fill (usually via `lapply()` or `vapply()`) rather
than copy-and-append. If creating a vector or list of results, use
`lapply()` (to create a list) or `vapply()` (to create a vector)
rather than a `for` loop. For instance,

    n <- 10000
    x <- vapply(seq_len(n), function(i) {
        ## ...
    }, integer(1))

manages the memory allocation of `x` and compactly represents the
transformation being performed. A `for` loop might be appropriate if
the iteration has side effects (e.g., displaying a plot) or where
calculation of one value depends on a previous value. When creating a
vector in a `for` loop, always pre-allocate the result

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
loop, making `n^2 / 2` total copies, and scale very poorly even for
trivial computations:

    > system.time(not_this(1000)
       user  system elapsed 
      0.004   0.000   0.004 
    > system.time(not_this(10000))
       user  system elapsed 
      0.169   0.000   0.168 
    > system.time(not_this(100000))
       user  system elapsed 
     22.827   1.120  23.936 

### Avoid `1:n` style iterations

Write `seq_len(n)` or `seq_along(x)` rather than `1:n` or
`1:length(x)`. This protects against the case when `n` or `length(x)`
is 0 (which often occurs as an unexpected 'edge case' in real code) or
negative.

### Re-use existing functionality

For common input formats see [common _Bioconductor_ import and classes][]

If there are problems, e.g., in performance or parsing your particular
file type, ask for input from other developers on the bioc-devel
mailing list. Common disadvantages to 'implementing your own' are the
introduction of non-standard data representations (e.g., neglecting to
translate coordinate systems of file formats to Bioconductor objects)
and user bewilderment.

### Re-use existing classes

Re-use enhances interoperability between _Bioconductor_ packages while
providing robust code for data manipulation.

Use `GenomicRanges::GRanges` (and `GRangesList`) to represent 1-based,
closed-interval genomic coordinates.

Use `SummarizedExperiment::SummarizedExperiment` (with or without
ranges as row data) to coordinate rectangular feature x sample data
(e.g., RNAseq count matrix) with feature and sample description. Use
`SummarizedExperiment` rather than the older `ExpressionSet`,
especially for sequence data.

For more existing classes see [common _Bioconductor_ import and classes][]

### Essential S4 interface

For any class you define, implement and use a 'constructor' for object
creation. A constructor is usually plain-old-function (rather than,
e.g., a generic with methods). It provides documented and
user-friendly arguments, while allowing for developer-friendly
implementation. Use the constructor throughout your own code,
examples, and vignette.

Implement a `show()` method to effectively convey information to your
users without overwhelming them with detail.

Accessors (simple functions that return components of your object)
rather than direct slot access (using `@`) help to isolate the
implementation of your class from its interface. Generally `@` should
only be used in an accessor, all other code should use the
accessor. The accessor does not need to be exported from the class if
the user has no need or business accessing parts of your data
object. Plain-old-functions (rather than generic + method) are often
sufficient for accessors; it's often useful to employ (consistently) a
lightweight name mangling scheme (e.g., starting the accessor method
name with a 2 or 3 letter acronym for your package) to avoid name
collisions between similarly named functions in other packages.

[microbenchmark]: https://cran.r-project.org/web/packages/microbenchmark
[unit tests]: /developers/how-to/unitTesting-guidelines/
[common _Bioconductor_ import and classes]: /developers/how-to/commonImportsAndClasses