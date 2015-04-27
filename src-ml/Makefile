# Makefile for cddml and cddmlgmp
# 
# see README
#
# make cddml	to produce binary cddml which uses floating-point arithmetic
# make cddmlgmp	to produce binary cddmlgmp which uses GMP rational exact arithmetic
# make all	to produce binaries both cddml and cddmlgmp
# make clean	to clean all objects
#
# Portions of this makefile require the use of GNU make.
# see http://www.gnu.org/software/make for more information.

MLINKDIR = /Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR}
    
LIBDIR = ${CADDSDIR}

MPREP = ${CADDSDIR}/mprep
RM = rm

EXTRA_LIBS = -stdlib=libstdc++ -lstdc++ -framework Foundation
# The following setting 
# EXTRA_LIBS = -lstdc++ -framework Foundation
# in Makefile for examples coming with Mma doesn’t work properly.

MATHLINK_LIB = -lMLi4

CC = /usr/bin/clang
# CXX = /usr/bin/clang++

BINARIES = cddml cddmlgmp

all : $(BINARIES)
	

cddml : cddmathlinktm.o cddmathlink.o cddmlio.o
	${CC} -I${INCDIR} cddmlio.o cddmathlinktm.o cddmathlink.o -L${LIBDIR} ${MATHLINK_LIB} ${EXTRA_LIBS} /usr/local/lib/libcdd.a -o cddml

cddmlgmp : cddmathlinktm.o cddmathlink_gmp.o cddmlio_gmp.o
	${CC} -I${INCDIR} cddmlio_gmp.o cddmathlinktm.o cddmathlink_gmp.o -L${LIBDIR} ${MATHLINK_LIB} ${EXTRA_LIBS} -lgmp /usr/local/lib/libcddgmp.a -o cddmlgmp

cddmathlinktm.o: cddmathlinktm.c
	$(CC) -c cddmathlinktm.c

cddmathlink.o: cddmathlink.c
	$(CC) -c cddmathlink.c

cddmathlink_gmp.o: cddmathlink.c
	$(CC) -DGMPRATIONAL -c cddmathlink.c -o cddmathlink_gmp.o

cddmlio.o: cddmlio.c
	$(CC) -c cddmlio.c

cddmlio_gmp.o: cddmlio.c
	$(CC) -DGMPRATIONAL -c cddmlio.c -o cddmlio_gmp.o

cddmathlinktm.c : cddmathlink.tm
	${MPREP} $? > $@

clean:
	${RM} -rf *.o *tm.c $(BINARIES)