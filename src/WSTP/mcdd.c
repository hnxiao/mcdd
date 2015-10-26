/* cddwstp.c: Main program to call the cdd library cddlib
   from Mathematica using WSTP,
   rewritten by Han Xiao based on cddml.c by Prof. Fukuda,
   Oct 18, 2014
*/

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "setoper.h"
#include "cdd.h"
#include "wstp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

extern void allvertices(int m_input, int d_input, char *a_input);
extern void allvertices2(int m_input, int d_input, char *a_input);
extern void allfacets(int n_input, int d_input, char *g_input);
extern void allfacets2(int n_input, int d_input, char *g_input);

void dd_WSWriteAmatrix(dd_Amatrix, long, long);
void dd_WSWriteMatrix(dd_MatrixPtr);
void dd_WSWriteSet(set_type);
void dd_WSWriteSetFamily(dd_SetFamilyPtr);
void dd_WSWriteError(dd_PolyhedraPtr);
void dd_WSSetMatrixWithString(dd_rowrange, dd_colrange, char *,dd_MatrixPtr);
char *dd_WSGetStrForNumber(mytype);

/*
 * WSTP visible functions.
 */

void allvertices(int m_input, int d_input, char *a_input)
/* output vertices and incidences */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr GI=NULL;
  dd_rowrange i,m;
  dd_colrange j,d;
  dd_ErrorType err;

  m=(dd_rowrange)m_input; d=(dd_colrange)d_input;
  A=dd_CreateMatrix(m,d);
/*
  for (i=0; i<m; i++){
    for (j=0; j<d; j++) dd_set_d(A->matrix[i][j],a_input[i*d+j]);
  }
*/
  dd_WSSetMatrixWithString(m, d, a_input, A);

  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);
    /* compute the second (generator) representation */
  if (err==dd_NoError) {
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);

    WSPutFunction(stdlink,"List",2);
    dd_WSWriteMatrix(G);
    dd_WSWriteSetFamily(GI);
  } else {
    dd_WSWriteError(poly);
  }

  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(GI);
}


void allvertices2(int m_input, int d_input, char *a_input)
/* output vertices, incidences and adjacency */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr GI=NULL,GA=NULL;
  dd_SetFamilyPtr AI=NULL,AA=NULL;
  dd_rowrange i,m;
  dd_colrange j,d;
  dd_ErrorType err;

  m=(dd_rowrange)m_input; d=(dd_colrange)d_input;
  A=dd_CreateMatrix(m,d);
/*
  for (i=0; i<m; i++){
    for (j=0; j<d; j++) dd_set_d(A->matrix[i][j],a_input[i*d+j]);
  }
*/
  dd_WSSetMatrixWithString(m, d, a_input, A);

  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);
    /* compute the second (generator) representation */
  if (err==dd_NoError){
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);
    GA=dd_CopyAdjacency(poly);
    AI=dd_CopyInputIncidence(poly);
    AA=dd_CopyInputAdjacency(poly);

    WSPutFunction(stdlink,"List",5);
    dd_WSWriteMatrix(G);
    dd_WSWriteSetFamily(GI);
    dd_WSWriteSetFamily(GA);
    dd_WSWriteSetFamily(AI);
    dd_WSWriteSetFamily(AA);
  } else {
    dd_WSWriteError(poly);
  }
  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(GI);
  dd_FreeSetFamily(GA);
  dd_FreeSetFamily(AI);
  dd_FreeSetFamily(AA);
}


void allfacets(int n_input, int d_input, char *g_input)
/* output facets and incidences */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr AI=NULL;
  dd_rowrange i,n;
  dd_colrange j,d;
  dd_ErrorType err;

  n=(dd_rowrange)n_input; d=(dd_colrange)d_input;
  G=dd_CreateMatrix(n,d);
/*
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) dd_set_d(G->matrix[i][j],g_input[i*d+j]);
  }
*/
  dd_WSSetMatrixWithString(n, d, g_input, G);

  G->representation=dd_Generator;
  poly=dd_DDMatrix2Poly(G, &err);
    /* compute the second (inequality) representation */
  if (err==dd_NoError){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);

    WSPutFunction(stdlink,"List",2);
    dd_WSWriteMatrix(A);
    dd_WSWriteSetFamily(AI);
  } else {
    dd_WSWriteError(poly);
  }

  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(AI);
}


void allfacets2(int n_input, int d_input, char *g_input)
/* output facets, incidences and adjacency */
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A=NULL,G=NULL;
  dd_SetFamilyPtr AI=NULL, AA=NULL;
  dd_SetFamilyPtr GI=NULL, GA=NULL;
  dd_rowrange i,n;
  dd_colrange j,d;
  dd_ErrorType err;

  n=(dd_rowrange)n_input; d=(dd_colrange)d_input;
  G=dd_CreateMatrix(n,d);
/*
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) dd_set_d(G->matrix[i][j],g_input[i*d+j]);
  }
*/
  dd_WSSetMatrixWithString(n, d, g_input, G);

  G->representation=dd_Generator;
  poly=dd_DDMatrix2Poly(G, &err);
    /* compute the second (inequality) representation */
  if (err==dd_NoError){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);
    AA=dd_CopyAdjacency(poly);
    GI=dd_CopyInputIncidence(poly);
    GA=dd_CopyInputAdjacency(poly);

    WSPutFunction(stdlink,"List",5);
    dd_WSWriteMatrix(A);
    dd_WSWriteSetFamily(AI);
    dd_WSWriteSetFamily(AA);
    dd_WSWriteSetFamily(GI);
    dd_WSWriteSetFamily(GA);
  } else {
    dd_WSWriteError(poly);
  }

  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  dd_FreeSetFamily(AI);
  dd_FreeSetFamily(AA);
  dd_FreeSetFamily(GI);
  dd_FreeSetFamily(GA);
}

/*
 * Utilities.
 */

void dd_WSWriteAmatrix(dd_Amatrix A, long rowmax, long colmax)
{
    long i,j;
    double a;
    char *str=NULL;
    
    if (A==NULL){
        rowmax=0; colmax=0;
    }
    WSPutFunction(stdlink,"List",rowmax);
    for (i=0; i < rowmax; i++) {
        WSPutFunction(stdlink,"List",colmax);
        for (j=0; j < colmax; j++) {
#if defined GMPRATIONAL
            str=dd_WSGetStrForNumber(A[i][j]);
            WSPutString(stdlink, str);
            if (str!=NULL) free(str);
#else
            a=dd_get_d(A[i][j]);
            WSPutDouble(stdlink, a);
#endif
        }
    }
}


void dd_WSWriteMatrix(dd_MatrixPtr M)
{
    WSPutFunction(stdlink,"List",2);
    dd_WSWriteAmatrix(M->matrix, M->rowsize, M->colsize);
    dd_WSWriteSet(M->linset);
}


void dd_WSWriteSet(set_type S)
{
    long j;
    
    WSPutFunction(stdlink,"List",set_card(S));
    for (j=1; j <= S[0]; j++) {
        if (set_member(j, S)) WSPutLongInteger(stdlink, j);
   
    }
}

void dd_WSWriteSetFamily(dd_SetFamilyPtr F)
{
    long i,j;
    
    if (F!=NULL){
        WSPutFunction(stdlink,"List",F->famsize);
        for (i=0; i < F->famsize; i++) {
            WSPutFunction(stdlink,"List",set_card(F->set[i]));
            for (j=1; j <= F->setsize; j++) {
                if (set_member(j, F->set[i])) WSPutLongInteger(stdlink, j);
            }
        }
    }
}


void dd_WSWriteError(dd_PolyhedraPtr poly)
{
    WSPutFunction(stdlink,"List",3);
    WSPutFunction(stdlink,"List",0);
    WSPutFunction(stdlink,"List",1);
    WSPutString(stdlink,"Error occured: code");
    WSPutFunction(stdlink,"List",1);
    WSPutInteger(stdlink,poly->child->Error);
}


char *dd_WSGetStrForNumber(mytype x)
{
    /* This is to make a string of rational expression for GMP rational number x.
     It does nothing (return NULL) if x is double.  */
    
    char *sd=NULL,*sn=NULL,*st=NULL;
#if defined GMPRATIONAL
    mpz_t zn,zd;
    int len;
    
    mpz_init(zn); mpz_init(zd);
    mpq_canonicalize(x);
    mpq_get_num(zn,x);
    mpq_get_den(zd,x);
    if (mpz_sgn(zn)==0){
        st=(char *)malloc(2);
        strcpy(st," 0");
    } else if (mpz_cmp_ui(zd,1U)==0){
        st=mpz_get_str(st,10,zn);
    } else {
        sn=mpz_get_str(sn,10,zn);sd=mpz_get_str(sd,10,zd);
        len=strlen(sn)+strlen(sd)+2;
        st=(char *)malloc(len);
        strcpy(st,sn);
        strcat(st,"/");strcat(st,sd);
    }
    mpz_clear(zn); mpz_clear(zd);
    if (sd!=NULL) free(sd);
    if (sn!=NULL) free(sn);
#else
    /* do nothing */
#endif
    /* printf("String for Number =%s\n",st);  */
    return st;
}


void dd_WSSetMatrixWithString(dd_rowrange m, dd_colrange d, char line[], dd_MatrixPtr M)
{
    dd_rowrange i=0;
    dd_colrange j=0;
    char *next,*copy;
    const char ct[]=", {}\n";  /* allows separators ","," ","{","}". */
    mytype value;
    double dval;
    
    dd_init(value);
    copy=(char *)malloc(strlen(line) + 4048); /* some extra space as buffer.  Somehow linux version needs this */
    strcpy(copy,line);
    next=strtok(copy,ct);
    i=0; j=0;
    do {
#if defined GMPRATIONAL
        dd_sread_rational_value(next, value);
#else
        dval=atof(next);
        dd_set_d(value,dval);
#endif
        dd_WriteNumber(stderr,value);
        dd_set(M->matrix[i][j],value);
        j++;
        if (j == d) {
            i++; j=0;
        }
    } while ((next=strtok(NULL,ct))!=NULL);
    free(copy);
    dd_clear(value);
    return;
}

/*
 * Passes control to WSTP.
 */

#if MACINTOSH_WSTP

int main( int argc, char* argv[])
{
	/* Due to a bug in some standard C libraries that have shipped with
	 * MPW, zero is passed to WSMain below.  (If you build this program
	 * as an MPW tool, you can change the zero to argc.)
	 */
    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	argc = argc; /* suppress warning */
	return WSMain( 0, argv);
}

#elif WINDOWS_WSTP

#if __BORLANDC__
#pragma argsused
#endif

int PASCAL WinMain( HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;

    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	hinstPrevious = hinstPrevious; /* suppress warning */

	if( !WSInitializeIcon( hinstCurrent, nCmdShow)) return 1;
	WSScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
	return WSMain( argv_end - argv, argv);
}

#else

int main(argc, argv)
	int argc; char* argv[];
{
    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	return WSMain(argc, argv);
}

#endif

/* end of cddwstp.c */
