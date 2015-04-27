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
#include "cddwsio.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

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
