# Debugging C code

Creation Date: 8 Feb 2012;
Last Edit Date: 10 Feb 2012

The following applies to non-Windows operating systems.  This is not
for the faint of heart, and requires some C-level familiarity.

## Set-up

For debugging package code at the C level, one usually starts (on
non-Windows) by installing the package without any compiler
optimizations, e.g., by following

    RShowDoc("R-admin")

section 6.3.3. setting for instance

    CFLAGS="-g -O0"

in .R/Makevars.

The next, essential, step is to write a short script that reliably and
quickly reproduces the error. Call this script `buggy.R`.

## gdb for interactive debugging

When the bug is isolated and easily produced in `buggy.R` start R with
a C-level debugger such as gdb (this is an old-school command-line
style debugger, not to everyone's taste)

    R -d gdb -f buggy.R

You'll end up at the gdb prompt

    (gdb)

and a typical operation is to (r)un or (c)ontinue execution

    (gdb) r

to run buggy.R. You'll end up back in C when there is a segfault, or
you press cntrl-C (`^C`, below), or when you've inserted a
(b)reakpoint at some C-level function that you suspect is buggy, e.g.,

    > ^C
    (gdb) b some_buggy_fun
    (gdb) c

When you do end up back in the debugger, you can print C variables or
the C representation of R variables (provided R isn't too confused by
this point)

    (gdb) p c_var
    (gdb) call Rf_PrintValue(some_R_variable)

You can also view a (b)ack(t)race of the call stack, navigate (u)p and
(d)own the call stack, etc. See

    (gdb) help

and our mutual friend Google for additional information.

## valgrind

The valgrind debugger can be used to spot memory access problems,
which are a common source of faults in C code. Start `R` with

    R -d valgrind -f buggy.R

This runs quite slowly, and will flag invalid memory read and write
locations. The former typically contribute to bad data, the later to
memory corruption and spectacular failures. The output will require C
familiarity to interpret. Again it is helpful to run the buggy code
with a package that has been installed without compiler
optimizations. See section 4.3 of `RShowDoc("R-exts")`

# Case study

As a case study, a colleague reported that their complicated program
would, on one particular computer, produce a segmentation fault or
just stop responding. The same series of actions wouldn't cause
problems on other computers. This sounds like a classic memory
problem, with the segfault and difficulty of reproduction. 

The first advice was to develop a simple script that reproduced the
problem: the original report had too many moving parts.  A big insight
was that the bug could be produced by running part of the code that
uses RCurl, followed by a call to the garbage collector, `gc()`. The
role of the garbage collector suggests again memory corruption of some
sort, and in particular that perhaps RCurl is allocating (at the C
level) an R object but not properly PROTECT'ing it from garbage
collection. We suspect RCurl rather than R or libcurl (other possible
players) because it is the least tested of the code. We could be
wrong, of course... After many iterations, my colleague arrived at
buggy24.R:

    library(RCurl)
    
    foo <- function() {
        url <- "https://google.com"
        curl <- getCurlHandle()
        opts <- list(followlocation=NULL, ssl.verifypeer=TRUE)
        d <- debugGatherer()
        getURL(url,customrequest="GET",curl=curl,debugfunction=d$update,.opts=opts)
    }
    
    execute <- function() {
        foo()
        gc()
    }
    
    execute()

This is pretty simple, and doesn't require access to any special
resources (like the server that was originally being queried). This
script doesn't cause a segfault when run on all systems, but running
valgrind (having installed RCurl without any optimizations) shows...

    > R -d valgrind -f buggy24.R
	...
    ==10859== Conditional jump or move depends on uninitialised value(s)
    ==10859==    at 0x11BF00F6: getCurlPointerForData (curl.c:798)
    ==10859==    by 0x11BF0E80: R_curl_easy_setopt (curl.c:164)
    ==10859==    by 0x11BF17AD: R_curl_easy_perform (curl.c:89)
    ==10859==    by 0x4ED5499: do_dotcall (dotcode.c:588)
    ==10859==    by 0x4F1CAA4: Rf_eval (eval.c:593)
    ==10859==    by 0x4F2BD5C: do_set (eval.c:1828)
    ==10859==    by 0x4F1C8B7: Rf_eval (eval.c:567)
    ==10859==    by 0x4F2B957: do_begin (eval.c:1514)
    ==10859==    by 0x4F1C8B7: Rf_eval (eval.c:567)
    ==10859==    by 0x4F297E9: Rf_applyClosure (eval.c:960)
    ==10859==    by 0x4F1CBA5: Rf_eval (eval.c:611)
    ==10859==    by 0x4F2BD5C: do_set (eval.c:1828)
    
Look around the C source code in RCurl's curl.c, as suggested by the
backtrace, just to get oriented. Then do

    R -d gdb -f buggy24.R

to run the script under gdb. Run our test script

    (gdb) r

No error. Don't give up, set a break point

    (gdb) b curl.c:798

and run again

    (gdb) r
    Breakpoint 1, getCurlPointerForData (el=0x79e038,
        option=CURLOPT_WRITEFUNCTION, isProtected=FALSE, curl=0x1d9bdc0)
        at curl.c:798
    798    curl.c: No such file or directory.
    (gdb)

That 'no such file' means that gdb doesn't know where to find the
RCurl package src/ directory, so tell it and (l)ist the context, and
(p)rint the value of the C variable `isProtected`, which seems to be
the source of the valgrind warning

    (gdb) dir ~/tmp/RCurl/src
    (gdb) l
    793                        }
    794                    }
    795                }
    796                break;
    797              case CLOSXP:
    798              (gdb) l
    793                        }
    794                    }
    795                }
    796                break;
    797              case CLOSXP:
    798                  if(!isProtected) {
    799                R_PreserveObject(el);
    800                }
    801                ptr = (void *) el;
    802                break;
    (gdb) p isProtected
    $5 = FALSE

`isProtected` has a value (it has to!), and furthermore the value of
FALSE results in PROTECT'ing the object `el` across C calls (this is
what `R_PreserveObject` does). This is pretty interesting, because
we're aware that garbage collection triggers the segfault. valgrind is
telling us that the value of `isProtected` isn't actually the result
of an assignment, it could be the result of accessing an array out of
bounds. Let's head up the call stack and see where this value is
coming from

    (gdb) up
    #1  0x00007ffff426e273 in R_curl_easy_setopt (handle=0x15d9600,
        values=0x1445788, opts=0xf3d418, isProtected=0xb7d308, encoding=0x776db0)
        at curl.c:164
    164            val = getCurlPointerForData(el, opt, LOGICAL(isProtected)[ i % n ], obj);
    (gdb) l
    159        /* Loop over all the options we are setting. */
    160        for(i = 0; i < n; i++) {
    161            opt = INTEGER(opts)[i];
    162            el = VECTOR_ELT(values, i);
    163                 /* Turn the R value into something we can use in libcurl. */
    164            val = getCurlPointerForData(el, opt, LOGICAL(isProtected)[ i % n ], obj);
    165   
    166                    if(opt == CURLOPT_WRITEFUNCTION && TYPEOF(el) == CLOSXP) {
    167                data->fun = val; useData++;
    168                status =  curl_easy_setopt(obj, CURLOPT_WRITEFUNCTION, &R_curl_write_data);
    (gdb)


We're entering the function `getCurlPointerForData` with the value
`LOGICAL(isProtected)[ i % n ]`. Here, `isProtected` is now an R
object, not a C variable. Looking at the surrounding code, that `i %
n` doesn't look right -- it's probably meant to recycle `isProtected`
in the case where a shorter logical variable is provided than the
vector of elements requiring protection, but the value of `n` is not
necessarily the length of `isProtected`. Let's have a look at what
we've got, using a C-level R function `Rf_PrintValue` to print R
values (SEXP's) in an R fashion

    (gdb) p isProtected
    $1 = (SEXP) 0xaad8a0
    (gdb) call Rf_PrintValue(isProtected )
    [1] FALSE

`isProtected` is a logical vector of length 1.

    (gdb) p i
    $7 = 1
    (gdb) p n
    $8 = 6
    (gdb) p i % n
    $9 = 1

...and we're trying to access element 1 of it. But the C
representation of R vectors is zero-based, so the only valid value of
the index is 0 -- we're out of bounds! This could well be our bug, and
it's time to try fixing it (naively,
`LOGICAL(isProtected)[ i % LENGTH(isProtected) ]`) to confirm our
diagnosis, or report to the `packageDescription("RCurl")$Maintainer`
who might have a better sense of the overall structure and intention
of the code.
