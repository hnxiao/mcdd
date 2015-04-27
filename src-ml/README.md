April 27, 2015

--------------------------------------------------------
#*MathLink* interface to C-Libray *cddlib*
--------------------------------------------------------

1. This directory contains additional files necessary to create 
`cddml` and `cddmlgmp`. It does **NOT** follow the auto-config 
structure and the user must edit `Makefile` so that it conforms 
to the local *Mathematica* and *cddlib* installation.


2. *Mathematica MathLink* executables and library files must be 
installed properly. In *Mathematica* v10, the executables *mcc* 
and *mprep*, the libraries `libMLi3.a` and `libMLi4.a`, the header 
file `mathlink.h` are all installed in dictionary  
`/Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions`.


3. Before `cddml(gmp)` compilation, one must compile *cddlib* 
libraries `libcdd.a` and `libcddgmp.a` (and *GMP* libaray `libgmp.a`)
first. It is recommended to use *Homebrew* on *Macintosh* to compile 
*GMP* and *cddlib*, then all necessary libraries including `libgmp.a`,
`libcdd.a` and `libcddgmp.a` can be found in directory    
`/usr/local/lib`     
via symbolic links.


4. Examples for how to use *cddlib* can be found in subdirectory 
*examples* of *mma-cdd*.
