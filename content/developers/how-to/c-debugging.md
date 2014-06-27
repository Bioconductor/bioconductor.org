# Debugging C/C++ code

Creation Date: 8 Feb 2012;
Last Edit Date: 27 Jun 2014

The following applies to non-Windows operating systems.  This is not
for the faint of heart, and requires some C-level familiarity.

For those who learn better visually, see the video by a Bioconductor
alum about [using gdb to debug R packages with native
code](http://vimeo.com/11937905).

The benefit of the [diagnose-a-crash](#diagnose-a-crash) and [Case
study](#case-study) examples is all the steps and logic are written
out; one does not need to rewind a video to review the steps.

## Set-up

The first, essential, step is to write a short script that reliably and
quickly reproduces the error. Call this script `buggy.R`.

For debugging package code at the C/C++ level, one usually starts by
installing the package without any compiler optimizations, e.g., by
following

    RShowDoc("R-admin")

section 6.3.3. setting for instance

    CFLAGS=-ggdb -O0

in .R/Makevars. See the relevant [Package Guidelines
section](/developers/package-guidelines/#development-compiler-settings)
for more examples and information.

## Detecting memory errors (Valgrind)

[Valgrind](http://valgrind.org/) is a mature suite of tools for
low-level program analysis. Valgrind's [memory error checker
(Memcheck)](http://valgrind.org/docs/manual/mc-manual.html) is the
premier tool for diagnosing C/C++ memory errors.

Valgrind can be used to spot memory access problems, which are a
common source of segfaults in C/C++ code. When the bug is isolated and
easily produced in `buggy.R`, start `R` with:

    R -d valgrind -f buggy.R

This runs quite slowly, and will flag invalid memory read and write
locations. The former typically contribute to bad data, the latter to
memory corruption and spectacular failures. The output will require C
familiarity to interpret. It is helpful to run the buggy code with a
package that has been installed without compiler optimizations. See
section 4.3 of `RShowDoc("R-exts")` and the relevant [Package
Guidelines
section](/developers/package-guidelines/#development-compiler-settings).

## Interactive debugging (gdb or lldb)

If you have never used a command-line debugger, there are many fine
quick-start guides on the web; it is not as daunting as it seems.

On Linux the preferred debugger is
[gdb](http://www.gnu.org/software/gdb/), but
[lldb](http://lldb.llvm.org/) is default for the Mavericks
platform. The interfaces are similar, but if you are accustomed to
gdb, see the [gdb to lldb command
map](http://lldb.llvm.org/lldb-gdb.html).

Start `R` with a C-level debugger such as gdb.

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

<a name="diagnose-a-crash"></a>
### Finding clues for causes of a program crash

Perhaps the most useful function of debuggers is to provide a
breadcrumb trail ("backtrace") of the routines that led to the crash
of a buggy program. With this knowledge we can narrow our inquiry to
code that affects parts of the program state relevant at the time of
the crash.

It is worth reiterating that it is *essential* optimizations be turned
off and the compiler be instructed to include debugging symbols if one
hopes to have a fruitful debugging session. See the relevant [Package
Guidelines
section](/developers/package-guidelines/#development-compiler-settings).

Though the debugger output in the example may differ slightly compared
to output from other computing environments, the underlying techniques
apply for diagnosing program crashes on any platform. See the [Case
study](#case-study) for a real-world example that uses Valgrind and
gdb in conjunction.

We will use a contrived example to demonstrate how to identify
potential places in our code that cause a crash. You should be able to
use the example files exactly as they appear. Some extraneous output
has been omitted for brevity.

* This particular example uses lldb on code compiled with clang on
  Mavericks. The example does not lead to a program crash on Linux
  with GCC or clang.

C++ file `buggy.cpp`:

    #include <map>
    #include <utility>
    #ifdef __cplusplus
    #define R_NO_REMAP
    #endif
    #include <Rdefines.h>

    extern "C" SEXP buggy_function();

    SEXP buggy_function() {
        std::map<int, int> m;
        m.insert(std::make_pair(5, 7));
        m.insert(std::make_pair(9, 42));
        std::map<int,int>::const_iterator it = m.begin();
        ++it;
        ++it;
        ++it;

        return R_NilValue;
    }

Compile with `R CMD SHLIB buggy.cpp -o buggy.so`.

`source()`ing this file (`buggy.R`) within an `R` session (or entering
the commands in an `R` session) will lead to a program crash:

    dyn.load("buggy.so")
    .Call("buggy_function")

Unfortunately `R`'s diagnostics are not very illuminating:

    > source("buggy.R")

     *** caught segfault ***
    address 0x2, cause 'memory not mapped'

    Traceback:
     1: .Call("buggy_function")
     2: eval(expr, envir, enclos)
     3: eval(ei, envir)
     4: withVisible(eval(ei, envir))
     5: source("buggy.R")

Now we turn to the debugger. Start `R` with the `lldb` debugger (or
equivalent for your platform):

    R -d lldb
    (lldb) run
    ## R startup messages elided
    ## now in R session
    > source("buggy.R")

At this point `R` crashes, lldb produces some output, and we are back
at the lldb prompt. The lldb output looks like this (showing us the
frame (#0) in the call stack where the crash occurred):

    Process 21657 stopped
    * thread #1: tid = 0xbcb4ab, 0x00000001028fcbb0 buggy.so`buggy_function [inlined] std::__1::__tree_node_base<void*>* std::__1::__tree_min<std::__1::__tree_node_base<void*>*>(std::__1::__tree_node_base<void*>*) at __tree:134, queue = 'com.apple.main-thread', stop reason = EXC_BAD_ACCESS (code=1, address=0x2)
        frame #0: 0x00000001028fcbb0 buggy.so`buggy_function [inlined] std::__1::__tree_node_base<void*>* std::__1::__tree_min<std::__1::__tree_node_base<void*>*>(std::__1::__tree_node_base<void*>*) at __tree:134
       131 	_NodePtr
       132 	__tree_min(_NodePtr __x) _NOEXCEPT
       133 	{
    -> 134 	    while (__x->__left_ != nullptr)
       135 	        __x = __x->__left_;
       136 	    return __x;
       137 	}

It looks like the debugger is telling us there was a memory access
error when fetching a tree node. (Trees are a common underlying data
structure for the standard library `map`). The output is voluminous
and looks confusing, but only the gist is important right now.

Still in the same lldb session, enter the `bt` command (for
"backtrace") at the lldb prompt and we see all the stack frames (and
function calls) that preceded the crash. Frames are listed in
ascending order, starting with the frame in which crash
happened. (Note frame #0 here is the same as frame #0 given above.)
This means when diagnosing a crash, it usually makes sense to start
with lower-numbered frames and proceed upwards.

    (lldb) bt
    * thread #1: tid = 0xbcb4ab, 0x00000001028fcbb0 buggy.so`buggy_function [inlined] std::__1::__tree_node_base<void*>* std::__1::__tree_min<std::__1::__tree_node_base<void*>*>(std::__1::__tree_node_base<void*>*) at __tree:134, queue = 'com.apple.main-thread', stop reason = EXC_BAD_ACCESS (code=1, address=0x2)
      * frame #0: 0x00000001028fcbb0 buggy.so`buggy_function [inlined] std::__1::__tree_node_base<void*>* std::__1::__tree_min<std::__1::__tree_node_base<void*>*>(std::__1::__tree_node_base<void*>*) at __tree:134
        frame #1: 0x00000001028fcbb0 buggy.so`buggy_function [inlined] std::__1::__tree_node_base<void*>* std::__1::__tree_next<std::__1::__tree_node_base<void*>*>(std::__1::__tree_node_base<void*>*) + 20 at __tree:158
        frame #2: 0x00000001028fcb9c buggy.so`buggy_function [inlined] std::__1::__tree_const_iterator<std::__1::__value_type<int, int>, std::__1::__tree_node<std::__1::__value_type<int, int>, void*>*, long>::operator++() at __tree:747
        frame #3: 0x00000001028fcb9c buggy.so`buggy_function [inlined] std::__1::__map_const_iterator<std::__1::__tree_const_iterator<std::__1::__value_type<int, int>, std::__1::__tree_node<std::__1::__value_type<int, int>, void*>*, long> >::operator++() at map:750
        frame #4: 0x00000001028fcb9c buggy.so`buggy_function + 188 at buggy.cpp:17
        frame #5: 0x0000000100073a13 libR.dylib`do_dotcall(call=<unavailable>, op=<unavailable>, args=<unavailable>, env=<unavailable>) + 323 at dotcode.c:578

frame #5 mentions `do_dotcall`, which is the native function (in the
`R` library) that corresponds to the `.Call("buggy_function")` line in
`buggy.R` where we call our C entry point. We can reasonably conclude
the useful information for our bug is probably in frames #0-4.

Here is a possible chain of thought that leads to the right conclusion:

1. frames #0-2 look like they are dealing with tree/map internals;
   ignore for the moment.

3. frame #3 indicates we are probably talking about our map
   const_iterator variable declared at line #14 in buggy.cpp
   (`std::map<int,int>::const_iterator it = m.begin();`).

4. frame #4 is the key: it tells us line (#17) in the `buggy.cpp` file
   (`++it;`) is where execution went from C++ code *we* wrote into the
   map iterator internals that produced the error.

5. Eureka! By carefully reading the code in `buggy.cpp` we realize
   that after the inserts the size of map `m` is 2. That means after
   incrementing iterator `it` at line #16 (`++it;`), the value of `it`
   is the [special
   *past-the-end*](http://www.cplusplus.com/reference/map/map/end/)
   value.  Incrementing an iterator beyond *past-the-end* (the third
   `++it;` at line #17) is [undefined
   behavior](/developers/how-to/mavericks-howto/#undefined-behavior)!

If we revise `buggy.cpp` to not increment `it` beyond *past-the-end*
by removing the third `++it;` the program runs without
complaint. Problem solved!

As you can see, the debugger was not able to immediately tell us *why*
the program crashed, just *where* the program crashed. We used the
information about where the crash happened to home in on the parts of
our code that affected program state relevant at the time of the
crash. Obviously this example is contrived; in a real-world scenario
the extra help afforded by the insight about relevant program state is
invaluable.

<a name="case-study"></a>
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
