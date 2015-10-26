/* 
 * cddml.c: Main program to call the cdd library cddlib
 * from Mathematica using MATHLINK.
 */

/* 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "setoper.h"
#include "cdd.h"
#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

extern void allvertices(int m_input, int d_input, char *a_input);
extern void allvertices2(int m_input, int d_input, char *a_input);
extern void allfacets(int n_input, int d_input, char *g_input);
extern void allfacets2(int n_input, int d_input, char *g_input);

void dd_MLWriteAmatrix(dd_Amatrix, long, long);
void dd_MLWriteMatrix(dd_MatrixPtr);
void dd_MLWriteSet(set_type);
void dd_MLWriteSetFamily(dd_SetFamilyPtr);
void dd_MLWriteError(dd_PolyhedraPtr);
void dd_MLSetMatrixWithString(dd_rowrange, dd_colrange, char *,dd_MatrixPtr);
char *dd_MLGetStrForNumber(mytype);

/*
 * MathLink visible functions.
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
  dd_MLSetMatrixWithString(m, d, a_input, A);

  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);
    /* compute the second (generator) representation */
  if (err==dd_NoError) {
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);

    MLPutFunction(stdlink,"List",2);
    dd_MLWriteMatrix(G);
    dd_MLWriteSetFamily(GI);
  } else {
    dd_MLWriteError(poly);
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
  dd_MLSetMatrixWithString(m, d, a_input, A);

  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);
    /* compute the second (generator) representation */
  if (err==dd_NoError){
    G=dd_CopyGenerators(poly);
    GI=dd_CopyIncidence(poly);
    GA=dd_CopyAdjacency(poly);
    AI=dd_CopyInputIncidence(poly);
    AA=dd_CopyInputAdjacency(poly);

    MLPutFunction(stdlink,"List",5);
    dd_MLWriteMatrix(G);
    dd_MLWriteSetFamily(GI);
    dd_MLWriteSetFamily(GA);
    dd_MLWriteSetFamily(AI);
    dd_MLWriteSetFamily(AA);
  } else {
    dd_MLWriteError(poly);
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
  dd_MLSetMatrixWithString(n, d, g_input, G);

  G->representation=dd_Generator;
  poly=dd_DDMatrix2Poly(G, &err);
    /* compute the second (inequality) representation */
  if (err==dd_NoError){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);

    MLPutFunction(stdlink,"List",2);
    dd_MLWriteMatrix(A);
    dd_MLWriteSetFamily(AI);
  } else {
    dd_MLWriteError(poly);
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
  dd_MLSetMatrixWithString(n, d, g_input, G);

  G->representation=dd_Generator;
  poly=dd_DDMatrix2Poly(G, &err);
    /* compute the second (inequality) representation */
  if (err==dd_NoError){
    A=dd_CopyInequalities(poly);
    AI=dd_CopyIncidence(poly);
    AA=dd_CopyAdjacency(poly);
    GI=dd_CopyInputIncidence(poly);
    GA=dd_CopyInputAdjacency(poly);

    MLPutFunction(stdlink,"List",5);
    dd_MLWriteMatrix(A);
    dd_MLWriteSetFamily(AI);
    dd_MLWriteSetFamily(AA);
    dd_MLWriteSetFamily(GI);
    dd_MLWriteSetFamily(GA);
  } else {
    dd_MLWriteError(poly);
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

void dd_MLWriteAmatrix(dd_Amatrix A, long rowmax, long colmax)
{
    long i,j;
    double a;
    char *str=NULL;
    
    if (A==NULL){
        rowmax=0; colmax=0;
    }
    MLPutFunction(stdlink,"List",rowmax);
    for (i=0; i < rowmax; i++) {
        MLPutFunction(stdlink,"List",colmax);
        for (j=0; j < colmax; j++) {
#if defined GMPRATIONAL
            str=dd_MLGetStrForNumber(A[i][j]);
            MLPutString(stdlink, str);
            if (str!=NULL) free(str);
#else
            a=dd_get_d(A[i][j]);
            MLPutDouble(stdlink, a);
#endif
        }
    }
}


void dd_MLWriteMatrix(dd_MatrixPtr M)
{
    MLPutFunction(stdlink,"List",2);
    dd_MLWriteAmatrix(M->matrix, M->rowsize, M->colsize);
    dd_MLWriteSet(M->linset);
}


void dd_MLWriteSet(set_type S)
{
    long j;
    
    MLPutFunction(stdlink,"List",set_card(S));
    for (j=1; j <= S[0]; j++) {
        if (set_member(j, S)) MLPutLongInteger(stdlink, j);
   
    }
}

void dd_MLWriteSetFamily(dd_SetFamilyPtr F)
{
    long i,j;
    
    if (F!=NULL){
        MLPutFunction(stdlink,"List",F->famsize);
        for (i=0; i < F->famsize; i++) {
            MLPutFunction(stdlink,"List",set_card(F->set[i]));
            for (j=1; j <= F->setsize; j++) {
                if (set_member(j, F->set[i])) MLPutLongInteger(stdlink, j);
            }
        }
    }
}


void dd_MLWriteError(dd_PolyhedraPtr poly)
{
    MLPutFunction(stdlink,"List",3);
    MLPutFunction(stdlink,"List",0);
    MLPutFunction(stdlink,"List",1);
    MLPutString(stdlink,"Error occured: code");
    MLPutFunction(stdlink,"List",1);
    MLPutInteger(stdlink,poly->child->Error);
}


char *dd_MLGetStrForNumber(mytype x)
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


void dd_MLSetMatrixWithString(dd_rowrange m, dd_colrange d, char line[], dd_MatrixPtr M)
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
 * Passes control to MathLink.
 */

#if MACINTOSH_MATHLINK

int main( int argc, char* argv[])
{
	/* Due to a bug in some standard C libraries that have shipped with
	 * MPW, zero is passed to MLMain below.  (If you build this program
	 * as an MPW tool, you can change the zero to argc.)
	 */
    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	argc = argc; /* suppress warning */
	return MLMain( 0, argv);
}

#elif WINDOWS_MATHLINK

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

	if( !MLInitializeIcon( hinstCurrent, nCmdShow)) return 1;
	MLScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain( argv_end - argv, argv);
}

#else

int main(argc, argv)
	int argc; char* argv[];
{
    dd_set_global_constants();  /* First, this must be called to use cddlib. */

	return MLMain(argc, argv);
}

#endif

/* end of cddml.c */
