/* cddwsio.h: Header file for WSTP/IO cddwsio.c,
   rewritten by Han Xiao based on cddmlio.h by Prof. Fukuda,
   Oct 18, 2014
*/

/* cddlib.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron
   P= {x :  b - A x >= 0}.
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#ifndef  __CDDWSIO_H
#define  __CDDWSIO_H
#endif  /* __CDDWSIO_H */

#ifndef  __CDD_H
#include "cdd.h"
#endif  /* __CDD_H */

/* ---------- FUNCTIONS MEANT TO BE PUBLIC ---------- */

/* basic IO */

void dd_WSWriteAmatrix(dd_Amatrix, long, long);
void dd_WSWriteMatrix(dd_MatrixPtr);
void dd_WSWriteSet(set_type);
void dd_WSWriteSetFamily(dd_SetFamilyPtr);
void dd_WSWriteError(dd_PolyhedraPtr);
void dd_WSSetMatrixWithString(dd_rowrange, dd_colrange, char *,dd_MatrixPtr);
char *dd_WSGetStrForNumber(mytype);


/* end of cddwsio.h */
