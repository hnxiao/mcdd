# Makefile for cddwstp and cddwstpgmp
#
# see README
#
# make cddwstp      to produce binary cddwstp which uses floating-point arithmetic
# make cddwstpgmp   to produce binary cddwstpgmp which uses GMP rational exact arithmetic
# make all          to produce binaries both cddwstp and cddwstpgmp
# make clean        to clean all objects
#
# Portions of this makefile require the use of GNU make.
# see http://www.gnu.org/software/make for more information.

WSTPLINKDIR = /Applications/Mathematica.app/SystemFiles/Links/WSTP/DeveloperKit
SYS = MacOSX-x86-64
CADDSDIR = ${WSTPLINKDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

WSPREP = ${CADDSDIR}/wsprep
RM = rm

EXTRA_LIBS = -stdlib=libstdc++ -lstdc++ -framework Foundation
# The following setting
# EXTRA_LIBS = -lstdc++ -framework Foundation
# in Makefile for examples coming with Mma doesn’t work properly.

WSTP_LIB = -lWSTPi4

CC = /usr/bin/clang
# CXX = /usr/bin/clang++

BINARIES = cddwstp cddwstpgmp

all : $(BINARIES)


cddwstp : cddwstptm.o cddwstp.o cddwsio.o
	${CC} -I${INCDIR} cddwsio.o cddwstptm.o cddwstp.o -L${LIBDIR} ${WSTP_LIB} ${EXTRA_LIBS} /usr/local/lib/libcdd.a -o cddwstp

cddwstpgmp : cddwstptm.o cddwstp_gmp.o cddwsio_gmp.o
	${CC} -I${INCDIR} cddwsio_gmp.o cddwstptm.o cddwstp_gmp.o -L${LIBDIR} ${WSTP_LIB} ${EXTRA_LIBS} -lgmp /usr/local/lib/libcddgmp.a -o cddwstpgmp

cddwstptm.o: cddwstptm.c
	$(CC) -c cddwstptm.c

cddwstp.o: cddwstp.c
	$(CC) -c cddwstp.c

cddwstp_gmp.o: cddwstp.c
	$(CC) -DGMPRATIONAL -c cddwstp.c -o cddwstp_gmp.o

cddwsio.o: cddwsio.c
	$(CC) -c cddwsio.c

cddwsio_gmp.o: cddwsio.c
	$(CC) -DGMPRATIONAL -c cddwsio.c -o cddwsio_gmp.o

cddwstptm.c : cddwstp.tm
	${WSPREP} $? > $@

clean:
	${RM} -rf *.o *tm.c $(BINARIES)
