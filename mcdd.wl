(* ::Package:: *)

mcddBegin::usage = "mcddBegin[] initializes resources.";
mcddBegin[] := (
	Clear[cddlib];
	Off[General::spell1];
	Off[General::spell];
    cddlib = Install["/Users/HXiao/GitHub/mcdd/src/WSTP/mcdd"];
  );


mcddEnd::usage = "McddEnd[] closes and destroys resources.";
mcddEnd[] := (
    Uninstall[cddlib];
    Clear[cddlib];
  );


PolyhedronVertexList::usage="PolyhedronVertexList[A,b] returns the vertex list of the polyhedron defined by linear inequality system Ax<=b.";
PolyhedronVertexList[A_List,b_List]:=Module[{M,rown,coln,vl},
	M=MapThread[Insert,{-A,b,Table[1,{Length[b]}]}];
	{rown,coln}=Dimensions[M];
	vl=AllVertices[rown, coln, ToString[Flatten@M]][[1,1]];
	Return[ToExpression/@vl]];


mcddBegin[];
Print["Connected to cddlib."];
