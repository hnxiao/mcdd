April 27, 2015

--------------------------------------------------------
#*WSTP* interface to C-Libray *cddlib* 
--------------------------------------------------------

1. This directory contains additional files necessary to create 
`cddwstp` and `cddwstpgmp`. It does **NOT** follow the auto-config 
structure and the user must edit `Makefile` so that it conforms 
to the local *Mathematica* and *cddlib* installation.


2. *Mathematica WSTP* executables and library files must be installed 
properly. In *Mathematica* v10, the executables `wscc` and `wsprep`, 
the libraries `libWSTPi3.a` and `libWSTPi4.a`, the header file `wstp.h`
are all installed in dictionary  
`/Applications/Mathematica.app/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions`.


3. Before `cddwstp(gmp)` compilation, one must compile *cddlib* libraries 
`libcdd.a` and `libcddgmp.a` (and *GMP* library `libgmp.a`) first. It is 
recommended to use *Homebrew* on *Macintosh* to compile *GMP* and *cddlib*, 
then all necessary libraries including `libgmp.a`, `libcdd.a` and 
`libcddgmp.a` can be found in directory    
`/usr/local/lib`     
via symbolic links.


4. Examples for how to use *cddlib* can be found in subdirectory 
*examples* of *mma-cdd*.
