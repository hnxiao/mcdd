# Makefile for mcdd
#
# see README
#
# make   	   to produce binary mcdd which uses GMP rational exact arithmetic
# make clean      to clean all objects
#
# Portions of this makefile require the use of GNU make.
# see http://www.gnu.org/software/make for more information.

WSTPDIR = /Applications/Mathematica.app/SystemFiles/Links/WSTP/DeveloperKit
SYS = MacOSX-x86-64
CADDSDIR = ${WSTPDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

WSPREP = ${CADDSDIR}/wsprep
RM = rm

EXTRA_LIBS = -stdlib=libstdc++ -lstdc++ -framework Foundation

WSTP_LIB = -lWSTPi4

CC = /usr/bin/clang
CXX = /usr/bin/clang++


mcdd : mcddtm.o mcdd.o
	${CC} -I${INCDIR} mcddtm.o mcdd.o -L${LIBDIR} ${WSTP_LIB} ${EXTRA_LIBS} -lgmp /usr/local/lib/libcddgmp.a -o mcdd

mcddtm.o: mcddtm.c
	$(CC) -c mcddtm.c

mcdd.o: mcdd.c
	$(CC) -DGMPRATIONAL -c mcdd.c -o mcdd.o

mcddtm.c : mcdd.tm
	${WSPREP} $? > $@

clean:
	${RM} -rf *.o *tm.c
