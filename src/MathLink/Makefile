# Makefile for mcdd
#
# see README
#
# make              to produce binary cdd
# make clean        to clean all objects
#
# Portions of this makefile require the use of GNU make.
# see http://www.gnu.org/software/make for more information.

MATHLINKDIR = /Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit
SYS = MacOSX-x86-64
CADDSDIR = ${MATHLINKDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = ${CADDSDIR}/mprep
RM = rm

EXTRA_LIBS = -stdlib=libstdc++ -lstdc++ -framework Foundation

MATHLINK_LIB = -lMLi4

CC = /usr/bin/clang
CXX = /usr/bin/clang++


mcdd : mcddtm.o mcdd.o
	${CC} -I${INCDIR} mcddtm.o mcdd.o -L${LIBDIR} ${MATHLINK_LIB} ${EXTRA_LIBS} -lgmp /usr/local/lib/libcddgmp.a -o mcdd

mcddtm.o: mcddtm.c
	$(CC) -c mcddtm.c

mcdd.o: mcdd.c
	$(CC) -DGMPRATIONAL -c mcdd.c -o mcdd.o

mcddtm.c : mcdd.tm
	${MPREP} $? > $@

clean:
	${RM} -rf *.o *tm.c 
