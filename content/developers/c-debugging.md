# Debugging C code

Creation Date: 8 Feb 2012

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

    gdb>

and a typical operation is to (r)un or (c)ontinue execution

    gdb> r

to run buggy.R. You'll end up back in C when there is a segfault, or
you press cntrl-C (`^C`, below), or when you've inserted a
(b)reakpoint at some C-level function that you suspect is buggy, e.g.,

    > ^C
    gdb> b some_buggy_fun
    gdb> c

When you do end up back in the debugger, you can print C variables or
the C representation of R variables (provided R isn't too confused by
this point)

    gdb> p c_var
    gdb> call Rf_PrintValue(some_R_variable)

You can also view a (b)ack(t)race of the call stack, navigate (u)p and
(d)own the call stack, etc. See

    gdb> help

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

Good luck! 
