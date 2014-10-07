# Coding Style

Creation Date: Antiquity.
Last Edit Date: 18 Feb 2014.

Using similar coding style helps people understand code.  These
guidelines are really just preferences; they are not enforced. 

* See also: [package guidelines](/developers/package-guidelines).

Indentation

* Use 4 spaces for indenting.  No tabs.
* No lines longer than 80 characters.

Variable names

* Use camelCaps: initial lowercase, then alternate case between words.

Function names

* Use camelCaps: initial lower case, then alternate case between words.
* Do not use '.' (in the S3 class system, `some(x)` where `x` is class
  `A` will dispatch to `some.A`).
  
Class names

* Use CamelCaps: initial upper case, then alternate case between
  words.

Use of space

* Always use space after a comma.  This: `a, b, c`.
* No space around "=" when using named arguments to functions. This:
  `somefunc(a=1, b=2)`
* Space around all binary operators: `a == b`.

Comments

* Use "##" to start full-line comments.
* Indent at the same level as surrounding code.

Namespaces

* Import all symbols used from packages other than "base". Except for
  default packages (base, graphics, stats, etc.) or when overly
  tedious, fully enumerate imports.
* Export all symbols useful to end users. Fully enumerate exports.
* For conditional package use, see
  [Package Guidelines](developers/package-guidelines/#dependencies)

Misc

* Use `<-` not `=` for assignment.
