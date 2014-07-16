<!--
# build me as follows:
library(knitr)
knit2html("mavericks-howto.md")
# generates mavericks-howto.html

-->


# C++/Mavericks Best Practices

This guide mostly addresses how to adapt C/C++ code in Bioconductor
packages to build on Mac OS X 10.9 (Mavericks). If your package does
not use C or C++ code, you can assume this document is not relevant to
you.

## Table of Contents

This document is written to be read from beginning to end, except for
the last section.

* [Orientation](#orientation): how to approach developing for
Mavericks.
* [C++11](#CXX11): important info about C++11 and clang's standard
library implementation.
* [Issues specific to the Mavericks
environment](#issues-with-Mavericks-env): special considerations for
tools available for Mavericks.
* [C++ Best Practices](#best-practices): generally applicable best
practices for C++ in Bioconductor packages.
* [Issues with code from external
sources](#issues-external-code-sources): tips for using code not
written by the package contributor.
* [Lessons learned from specific issues](#lessons-learned): a
repository of knowledge about specific problems.

You might skip to the [lessons-learned](#lessons-learned) section to
check if your issue has already been explored.

<a name="orientation"></a>
## Orientation

Note: For simplicity, this guide uses 'GCC' (GNU Compiler Collection)
and 'clang' to refer to each respective collection of tools (including
C++ compilers), rather than simply the 'GCC C compiler' or 'clang C
compiler'. 'Mavericks environment' refers to the combination of clang
and Xcode versions available by default for Mavericks.

With the release of the Mavericks
build of R, CRAN and Bioconductor have adopted Apple's preferred
toolchain for building packages for the Mavericks
platform. Bioconductor packages are built on the Mavericks platform
using the OS X default combination of [clang](http://clang.llvm.org/)
and [Xcode](https://developer.apple.com/xcode/).

The introduction of the Mavericks environment has revealed a number of
issues with building packages, most of which are due to C/C++ code
that relies too heavily on the GCC way of doing things. Most problems
revealed by the transition to Mavericks are caused by C++ coding
practices that are universally recognized as problematic, and are
addressed by adhering to established [best practices](#best-practices).

Here are some common sources of problems, posed as how the Mavericks
environment contrasts with GCC:

* [Undefined behavior](#undefined-behavior). Under GCC the undefined
  behavior might fail silently without affecting program execution,
  whereas the same code in the Mavericks environment leads to a
  segfault ("segmentation fault").

* [Issues related to naming](#name-resolution-errors), particularly in
  circumventing the protections offered by C++ namespaces. GCC header
  organization seems to be more forgiving of loose naming, whereas the
  Mavericks environment is more exacting.

* [C++11]((#CXX11) as the default language specification, and clang's
  [libc++](http://libcxx.llvm.org/). libc++ is an implementation of the
  C++11 standard library written from scratch. At the time of this
  writing the libc++ implementation is only available for Mac OS X.

* [Code not written directly by the package
  contributor(s)](#issues-external-code-sources). Many Bioconductor
  packages that do not transition well to the Mavericks environment
  include third-party (e.g., Boost) or generated (e.g., SWIG)
  code. Many external code sources make assumptions that are not valid
  for the Mavericks environment. The compounded difficulty of
  diagnosing problems in code not written by the package author limits
  the Bioconductor team's ability to help with your package.

### What is different about developing for the Mavericks environment?

The biggest change is the introduction of clang's
[libc++](http://libcxx.llvm.org/) implementation of the C++11 standard
library and the library headers associated with Xcode. clang is
gaining market share for several reasons, helped by the fact that
clang is intended only as a compiler for C-based languages. Advocates
of clang hold that clang:

* offers better diagnostic information for errors and warnings
* has quicker compilation times
* sometimes yields smaller binaries
* in some cases yields faster execution speeds (disputed)

For guidance on compiler flags to use while developing, see the
[relevant section of the Package Guidelines
page](/developers/package-guidelines/#c-code).

#### Differences in unspecified behavior, memory addressing policies

Some errors encountered during the transition seem to be attributable
to reliance on non-portable unspecified behaviors. See the section of
this guide about [unspecified behavior](#undefined-behavior) for more
information.

For example, many aspects of C/C++ memory addressing are
implementation-dependent, which means expected behavior is not
prescribed by the C/C++ standard ("unspecified") and is therefore up
to compiler writers to decide.

The foremost difference regarding memory is clang seems to be more
restrictive about out-of-bounds memory addressing.

<a name="find-bugs"></a>
#### How to find bugs

With GCC the preferred debugger is
[gdb](http://www.gnu.org/software/gdb/), but many prefer
[lldb](http://lldb.llvm.org/) for debugging in the Mavericks
environment. Support for [lldb on other
platforms](http://lldb.llvm.org/status.html) is limited at this time.

Because a number of bugs in packages are related to memory addressing
or layout errors, relying on a debugger alone might not be sufficient
to track down memory errors. [Valgrind](http://valgrind.org/) is the
premier tool for detecting memory errors.

See the Bioconductor guide on [debugging C/C++
code](/developers/how-to/c-debugging/) for examples of using a
debugger and Valgrind.

#### What if I cannot access a Mavericks machine?

There is no substitute for using a Mavericks machine to troubleshoot
packages that fail on Mavericks. Many of the errors seen are only
reproducible with the combination of clang, Xcode, and OS X 10.9.

But there are several options short of procuring a Mavericks machine:

1. As an *exploratory measure*, install a more recent version of GCC
and compile your package with `-std=c++11` or `-std=gnu++11` compiler
arguments (see [Package
Guidelines](/developers/package-guidelines#development-compiler-settings)
for info about development compiler flags); as of [version
4.8.1](https://gcc.gnu.org/gcc-4.8/cxx0x_status.html), GCC implements
all major features of the 2011 ISO C++ standard. Diagnostics for
errors and warnings have also greatly improved with recent GCC
versions. Using a [C++11](#CXX11) implementation might reveal warnings
or errors that point to the same issues encountered on Mavericks.

2. Install clang; this is of limited value because many errors are
unique to the Mavericks environment.

3. Use Valgrind for memory addressing problems; because so many errors
on the Mavericks platform are related to memory addressing problems,
many errors should be equally discoverable using Valgrind on Linux.

4. If you are unable to diagnose your problem using the combination of
[build system output](http://www.bioconductor.org/checkResults/) and
Valgrind, feel free to contact the [bioc-devel mailing
list](/help/mailing-list/).

<a name="CXX11"></a>
## C++11

Although the default version of clang on Mavericks includes support
for all C++11 features, Bioconductor support for C++11 is dependent on
the platform with the oldest toolchain. Because the current Snow
Leopard (Mac OS X 10.6.8) toolchain does not support any C++11
features, Bioconductor packages generally should not use C++11
features. Eventually, when Mavericks is more widely adopted, support
for Snow Leopard will be dropped.

C++11 is *not* completely backward-compatible with older standards.

It is possible to tell clang to use older versions of the standard
library (the default is [libc++](http://libcxx.llvm.org/)), but
relying on OS version-specific compilation settings is not a workable
long-term solution. This approach greatly increases the maintenance
burden for package authors and limits the Bioconductor team's ability
to offer support.

* An insidious and more catastrophic consequence of using non-default
  standard libraries is the issue of binary incompatibility. Packages
  linked against one standard library are liable to crash
  (mysteriously) when interfacing with packages linked against a
  different standard library. This is also true for programs at the OS
  level: any program compiled and linked against libstdc++ on the
  Mavericks platform is, by default, assumed to be incompatible with
  programs compiled and linked against libc++.

Code should be adapted to avoid constructs that are backward- or
forward-incompatible. See the [forward-incompatibility
problems](#forward-incompatible) section for examples.

<a name="issues-with-Mavericks-env"></a>
## Issues specific to Mavericks environment

### C Linkage

C++ uses `extern "C"` to give declarations C linkage, and hence make
the declarations accessible to C code. Some `R` headers when
`#include`d in C++ will include C++ system headers that should **not**
have C linkage. According to the relevant [Writing R Extensions Manual
section](http://cran.r-project.org/doc/manuals/R-exts.html#Interfacing-C_002b_002b-code),
`R` header files should **not** be included within `extern "C"`
blocks.

A typical symptom of bad linkage is at package *load time* (not
compilation or link time) an error says a particular symbol cannot be
found.

* C++ mangles names, so the symbol name `R` says it cannot find will
  often be unrecognizable. Use the `c++filt` program installed on your
  system or an online name demangler to produce a human-readable
  version of the symbol name. Note that mangled names are
  environment-specific so a demangler meant for GCC symbol names on
  Linux will not demangle clang symbol names from a Mac.

Solution: All `R` headers should be `#include`d **outside** `extern
"C"` blocks.

Example of correct `#include` of R headers:

    #include <R.h>
    extern "C" {
      void foo(); // function 'foo' and other code in this block has C linkage
      ...
    }
    extern "C" void bar(); // function 'bar' has C linkage

<a name="openmp"></a>
### OpenMP

As of this writing the Mavericks environment does not support OpenMP,
and it is unknown if the tools released by Apple ever will.

Code should not rely on the availability of OpenMP. Independent of
concerns over OpenMP support, code should be written from the start to
degrade nicely in a single-threaded environment.

See the [Writing R Extensions Manual
section](http://cran.r-project.org/doc/manuals/R-exts.html#OpenMP-support)
for information about OpenMP code in R packages, and detecting
support.

Solution: use preprocessor if-else directives so code degrades
gracefully if OpenMP support is not available:

```
#ifdef SUPPORT_OPENMP
    // multithreaded OpenMP version of code
#else
    // single-threaded version of code
#endif
```

See the
[ShortRead](http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html)
package for an example of good practices around support for OpenMP.

<a name="best-practices"></a>
## C++ Best Practices

These C++ practices are applicable for most C++ projects, but have
been identified by the Bioconductor community as particularly helpful
in avoiding issues in the Mavericks environment.

### Use Rcpp

The [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html)
CRAN package allows seamless integration of C++ with `R`, and is
cross-platform. The package affords many of the same benefits for the
`R` C interface that make C++ so appealing as a language, while
eliminating many of the pitfalls of programming to the `R` interface.

The package is well documented, and has an extensive repository of
working examples for many tasks: the [Rcpp
Gallery](http://gallery.rcpp.org/).

<a name="name-resolution-errors"></a>
### Avoid name resolution errors

A name resolution error occurs when the compiler encounters an
identifier (e.g., variable or function name) that is ambiguous (that
is, there is a "collision" between two or more identifiers), or name
lookup rules lead the compiler to resolve a name incorrectly. clang is
more exacting about identifiers.

A typical symptom of a name resolution error is the compiler complains
that the type or number of arguments a function got is different from
what it expects, and the compiler points to a C++ header file in the
standard library.

There are two primary issues with name resolution when writing R
packages:

#### Re-mapping of identifiers from R headers

For convenience, `R` aliases common identifiers from `R`
headers. E.g., `Rf_length(SEXP)` becomes `length(SEXP)`. While this
may be convenient for `C`, the organization of headers in the
Mavericks environment seems to lead to more collisions than with
GCC. See relevant section of the [Writing R Extensions
Manual](http://cran.r-project.org/doc/manuals/R-exts.html#The-R-API)

Solution: prevent re-mapping of `R` identifiers for C++ code by
defining the R_NO_REMAP symbol. This can be done at the package level
with a -DR_NO_REMAP preprocessor flag, or on a file-by-file basis with
`#define R_NO_REMAP`. Use fully-qualified versions of `R` identifiers,
usually by prepending `Rf_`.

Example excerpt from header file that prevents re-mapping:

    file CxxCode.h
    ------------------

    #ifndef CXX_CODE_H
    #define CXX_CODE_H
    #ifdef __cplusplus
    #define R_NO_REMAP
    #endif

    ...
    
    void foo(SEXP s) {
        if(Rf_length(s) > 1) // fully qualified: 'Rf_length'
            ...
    }

    #endif

Note `Rf_length` is just one example of the many `R` identifiers that
might conflict with names in C++ standard library headers.

#### Namespace hygiene

Namespaces were introduced in C++ to limit the incidence of name
collisions. Many authors new to C++, however, use
[using-directives](http://en.cppreference.com/w/cpp/language/namespace#Using-directives)
(particularly the '`using namespace std`' directive) unnecessarily,
thereby reintroducing the very problems namespaces are meant to solve.

As pointed out in the [cppreference Notes
section](http://en.cppreference.com/w/cpp/language/namespace#Notes)
the use of the `using namespace std` directive introduces the entire
std namespace for name resolution. There is a high likelihood that
among all the headers in the standard library there is an identifier
that conflicts with the identifiers in your package.

Solution: avoid the '`using namespace std`' directive completely if
possible, *especially* in header files. Prefer
[using-*declarations*](http://en.cppreference.com/w/cpp/language/namespace#Using-declarations)
over using-directives or simply use fully-qualified versions of
standard library identifiers. New C++ authors overestimate how much
including the scope resolution operator (namely '`std::`') affects
readability.

Example of introducing std namespace identifiers:

```
#include <map>
#include <utility>
// Suppose we want to access the std::map and std::make_pair identifiers
```

Many new C++ authors will use a using-directive to introduce the
identifiers they need. **Avoid** if possible:

```
using namespace std; // introduces entire std namespace for resolution
```

One alternative is using-declarations (e.g., '`using std::map;`',
which allow hand-picking of identifiers to introduce (as opposed to
the entire std namespace); here we just want `std::map` and
`std::make_pair`. Even if the list of identifiers we want is quite
long we just need a single using-declaration for each one:

```
using std::map; // 'map' and 'make_pair' introduced at declaration scope
using std::make_pair;
```

Using-declarations can also be block-scoped. This is preferred over
using-declarations at global scope, as it prevents the unnecessary
introduction of names at the global scope, a tenet of good namespace
hygiene:

```
void foo() {
     using std::map;
     using std::make_pair;
     map<int, int> m;
     m.insert(make_pair(5, 7));
     ...
}
```

A perfectly good alternative is to simply precede standard library
identifiers with '`std::`', which most C++ programmers are accustomed
to reading:

```
void foo() {
     std::map<int, int> m;
     m.insert(std::make_pair(5, 7));
}
```

<a name="undefined-behavior"></a>
### Avoid undefined behavior and non-portable unspecified behavior

There are two major categories of behavior that is not prescribed by
the C or C++ standards:

* *undefined behavior* is specified to be arbitrary; code that
   produces the behavior might cause the program to crash, or it
   might execute without complaint ("silently"). The effect can also
   differ from one program execution to the next. Well-known examples
   include division by zero, indexing outside of array bounds, and
   dereferencing a null pointer.

* *unspecified behavior* is consistent and documented, but decided by
   an implementation. These are behaviors that are either not
   mentioned at all by the respective standard, or are mentioned to
   say that they are implementation-dependent. Well-known examples
   include the size of the `int` type and the size of pointers.

A typical symptom of problematic undefined or unspecified behaviors is
a segfault that only appears in the Mavericks environment. The reason
the problem was not discovered before might be that GCC silently
allows the code to execute instead of crashing the program.

Solution: code defensively to avoid problematic constructs and [use
debuggers](#find-bugs) to find the code that leads to errors.

<a name="issues-external-code-sources"></a>
## Issues with code from external sources

Some packages need to use code not written directly by the
contributor(s). The most common scenario is to include the source code
of a library written by a third party. Some packages also use code
produced by code generation tools, e.g.,
[SWIG](http://www.swig.org/). First, see the relevant [Package
Guidelines section](../../package-guidelines#third-party-code) for
guidance on code from external sources.

See the [lessons learned](#lessons-learned) section of this guide for
suggestions about specific code sources.

### Generated Code

Some packages use code that is generated by third-party tools, i.e.,
code written by machines. SWIG is a common example.

A problem with code written by machines is that the code is meant to
be *read* by machines. The top of each code file produced by SWIG, for
example, states that the code is not meant to be read or edited by
hand.

Because many code generation tools make assumptions that are invalid
for the Mavericks environment, the code needs human attention to fix
errors; but because of the inscrutable nature of machine-written code
it is very difficult to isolate the errors.

Solution: re-generate problem code if possible, otherwise fix by
hand. Fixing by hand is **strongly discouraged**.

### Third-party code

Some packages include third-party libraries that were not written in a
compiler-independent way, and so do not build out-of-the-box in the
Mavericks environment.

Solutions in approximate order of preference:

1. Check if an existing [CRAN](http://cran.r-project.org/) or
   Bioconductor package provides the same functionality, while meeting
   the performance needs of your use case. Eliminating third-party
   code from your package greatly reduces the maintenance burden.

2. Check if the library has been updated. Some libraries with an
   active user community undergo updates that add support for more
   compilers/environments.

3. Check if the maintainers are aware the library does not work for
   the Mavericks environment, and find out if support is
   forthcoming. It is usually easy to directly contact authors of
   libraries maintained by an individual or a small group.

4. Use an actively maintained alternate library that provides
   equivalent functionality. Sometimes if a library is no longer
   maintained, it is because the library has been abandoned for an
   alternative project that provides the same functionality.

5. Update the library code included in your package by
   hand. **Strongly discouraged**. Maintainer assumes responsibility
   for keeping code up-to-date with mainline source project. If
   undertaken, record descriptions of changes that are needed so when
   the codebase is updated the changes can be easily reproduced.

<a name="lessons-learned"></a>
## Lessons learned from specific issues

This section serves as a loosely organized repository for knowledge
gained about specific problems and their solutions. It is not expected
to be comprehensive. Items will be added as the knowledge base
grows. Bioconductor is eager for suggestions; please write the
[bioc-devel mailing list](/help/mailing-list/) if you have any!

If you have no idea where to start for diagnosing your broken package
it might be worthwhile to skim over all of this section.

Where relevant, the issue is marked as being discoverable at compile
time or runtime.

Where applicable, a link to a live code demo is provided.

<a name="forward-incompatible"></a>
### Forward-incompatibility problems with C++11

C++11 is not completely backward-compatible. In particular, the API
for some parts of the standard library has changed slightly, in
perhaps subtle ways. Most issues require minimal tweaking to fix.

#### Container iterator const-ness

Type: *compile time*

A number of operations on standard library containers now require
iterators to be `const`. Two ready examples are the
[`insert`](http://www.cplusplus.com/reference/vector/vector/insert/)
and [`erase`](http://www.cplusplus.com/reference/vector/vector/erase/)
methods that take iterator parameters.

### Iterating standard library containers

Type: *runtime*

Generally, using the special
[*past-the-end*](https://gcc.gnu.org/onlinedocs/libstdc++/manual/iterators.html)
iterator value other than for equality checks (i.e., `==` or `!=`)
results in undefined behavior. Particularly, in the Mavericks
environment:

* Dereferencing an iterator with the *past-the-end* value results in a
  segfault
* Incrementing an iterator *beyond* the *past-the-end* value results
  in a segfault

For a walkthrough example of incrementing beyond *past-the-end*, see
the [diagnose-a-crash
example](/developers/how-to/c-debugging/#diagnose-a-crash) on the
debug C/C++ page.

<a name="external-code-sources"></a>
### External code sources

Common examples of external code sources are
[SWIG](http://www.swig.org/), [Boost](http://www.boost.org/), and
numerous file format libraries. Code from external sources is
sometimes written in non-compiler-independent way. Check the
documentation to see if the Mavericks environment is supported.

#### Boost

[Boost](http://www.boost.org/) is a source of free, peer-reviewed C++
libraries that enhance the language. Many parts of Boost are
"header-only", which means they do not need to be separately compiled
and the headers merely need to appear in the search path in order for
client code to use them.

Many Boost libraries are platform-independent, but not all. Some Boost
libraries are either in the process of adding Mavericks environment
support, or the library authors have announced Mavericks environment
support will *not* be added.

Solutions in approximate order of preference:

1. Use the [BH](http://cran.r-project.org/web/packages/BH/index.html)
   package on CRAN, if possible. The BH package provides several Boost
   header-only libraries. Using the BH package means the maintenance
   cost of using Boost in your package is virtually nothing.

2. Update the Boost libraries you include with your package. Boost
   libraries sometimes contain bugs, or are later updated to add
   support for other platforms. It is the Bioconductor package
   maintainer's responsibility to keep all code in the package
   updated.

3. Contact the authors of the specific Boost library. If you cannot
   find an announcement regarding support for the Mavericks
   environment, it might be worth contacting the library authors to
   inquire.

#### SWIG

SWIG generates code to interface between code written in C/C++ and
other languages. At the time of this writing, SWIG support for clang
is limited, and SWIG particularly has problems with clang's libc++
version of the (C++11) standard library. Some of the problems are
limited to issues that can be addressed by tweaking function
signatures. Other problems are deeply embedded in the way SWIG
produces code.

At the time of this writing, [this thread on the SWIG-devel mailing
list](http://comments.gmane.org/gmane.comp.programming.swig.devel/22966)
seems to have the most in-depth discussion about working with SWIG on
Mavericks.

Solutions in approximate order of preference:

1. Eliminate SWIG code, if possible. This will probably do the most
   for reducing maintenance burden.

2. Re-generate SWIG code with the newest version of SWIG. At the time
   of this writing SWIG was recently updated to include partial
   support for C++11, which might alleviate problems with clang's
   [libc++](http://libcxx.llvm.org/). See the [SWIG document about
   C++11 support](http://www.swig.org/Doc3.0/CPlusPlus11.html). It is
   possible the new version will not produce the problematic
   code. Note code must be valid for all supported compilers.

3. Troubleshoot and fix errors by hand. **Strongly discouraged**. If
   undertaken, record descriptions of the changes that were needed so
   if the code is regenerated the changes can be easily reproduced.

   Read [SWIG documentation](http://www.swig.org/doc.html) to find
   guidance about troubleshooting. (For example, the SWIG `-E` switch
   outputs results after the preprocessor has run.) Perhaps start by
   removing all SWIG functionality and gradually adding features. Find
   information on the web about how to fix the errors for your
   package.

#### f2c

[f2c](http://www.netlib.org/f2c/) is a tool that converts Fortran77
code to C/C++ code. The maintenance burden required to make f2c code
cross-platform is substantial. Since a fortran compiler (or emulator)
is required to install `R`, f2c is usually unnecessary. Several
packages use native fortran code without a problem.

Solution: If at all possible, remove the need for f2c. The only
recourse is to finesse makefiles to the point that each supported
platform more or less has a targeted makefile. Please write the
[bioc-devel mailing list](/help/mailing-list/) if you have trouble.
