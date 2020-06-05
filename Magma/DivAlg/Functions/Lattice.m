import "../../BasicData.m": K,BB,R,eps,dimsym;
import "../../Initialize.m": n,BBsym,IM,phi;

intrinsic Dagger(A::Any) -> AlgMatElt
 {The involution dagger.}
 return Transpose(A);
end intrinsic;

intrinsic InnerProduct(A::Mtrx,B::Mtrx) -> FldElt
 {Inner product of A and B.}
 return Evaluate(Trace(A*Dagger(B)),phi);
end intrinsic;

intrinsic PreInnerProduct(A::Mtrx,B::Mtrx) -> FldElt
 {Inner product of A and B in K.}
 return Trace(A*Dagger(B));
end intrinsic;

intrinsic Evaluate(A::Mtrx,x::Mtrx) -> Any
 {Evaluates the vector x at the matrix A via the inner product.}
 return InnerProduct(A,x*Dagger(x));
end intrinsic;

intrinsic PreEvaluate(A::Mtrx,x::Mtrx) -> Any
 {Value of the evaluation of x at A in the number field K.}
 return PreInnerProduct(A,x*Dagger(x));
end intrinsic;

intrinsic FormMinimum(A::DAForm) -> FldReElt
 {The minimum of the form A on the maximal order.}
 if assigned A`Minimum then return A`Minimum; end if;
 Gram:=MatrixRing(R,#IM)!0;
 for i in [1..#IM] do
  for j in [i..#IM] do
   Gram[i][j]:=InnerProduct(A`Matrix,IM[i]*Dagger(IM[j]));
   Gram[j][i]:=Gram[i][j];
  end for;
 end for;
 L:=LatticeWithGram(Gram);
 A`Minimum:=Minimum(L);
 return A`Minimum;
end intrinsic;

intrinsic PreFormMinimum(A::DAForm) -> FldReElt
 {The value of the minimum of the form A on the maximal order in the number field K.}
 if assigned A`PreMinimum then return A`PreMinimum; end if;
 v:=MinimalVectors(A)[1];
 value:=PreEvaluate(A`Matrix,v);
 A`PreMinimum:=value;
 return value;
end intrinsic;

intrinsic MinimalVectors(A::DAForm) -> SeqEnum
 {Returns the shortest vectors of the form A in M.}
 if assigned A`MinimalVectors then return A`MinimalVectors; end if;
 Gram:=MatrixRing(R,#IM)!0;
 for i in [1..#IM] do
  for j in [i..#IM] do
   Gram[i][j]:=InnerProduct(A`Matrix,IM[i]*Dagger(IM[j]));
   Gram[j][i]:=Gram[i][j];
  end for;
 end for;
 L:=LatticeWithGram(Gram);
 m:=Minimum(L);
 S:=ChangeRing(ShortVectorsMatrix(L,m*(1-eps),m*(1+eps)),Integers());

 vecs:=[];
 for i in [1..NumberOfRows(S)] do
  v:=Parent(IM[1])!0;
  for j in [1..NumberOfColumns(S)] do
   v+:=S[i][j]*IM[j];
  end for;
  Append(~vecs,v);
 end for;
 A`MinimalVectors:=vecs;
 v:=vecs[1];
 A`Minimum:=Evaluate(A`Matrix,v);
 return A`MinimalVectors;
end intrinsic;

intrinsic FindPerp(S::SeqEnum) -> AlgMatElt
 {Find perpendicular vector to the projections in the input list S.}
 SS:=S;
 T:=KMatrixSpace(K,#BBsym,#SS)!0;
 for i in [1..#BBsym] do
  for j in [1..#SS] do
   T[i][j]:=PreInnerProduct(BBsym[i],SS[j]);
  end for;
 end for;
 KM:=KernelMatrix(T);
 for i in [1..NumberOfRows(KM)] do
  km:=MatrixRing(K,n)!0;
  for j in [1..NumberOfColumns(KM)] do
   X:=BBsym[1];
   km+:=(K!(KM[i][j])*BBsym[j]);
  end for;
  if km ne 0 then
   return km;
  end if;
 end for;
 return false;
end intrinsic;

intrinsic PerfectionRank(A::DAForm) -> RngIntElt
 {Returns the perfection rank of the input form A.}
 if assigned A`PerfectionRank then return A`PerfectionRank; end if;
 L:=[x*Dagger(x) : x in MinimalVectors(A)];
 LL:=[Eltseq(L[i]) : i in [1..#L]];
 MM:=Matrix(LL);
 A`PerfectionRank:=Rank(MM);
 return A`PerfectionRank;
end intrinsic;

intrinsic PerfectionRankDAList(S::SeqEnum) -> Any
 {Returns the perfection rank of the projections onto the vectors contained in the input list S.}
 L:=[x*Dagger(x) : x in S];
 LL:=[Vector(&cat [Eltseq(x) : x in Eltseq(L[i])]) : i in [1..#L]];
 MM:=Matrix(LL);
 return Rank(MM);
end intrinsic;

intrinsic IsPositiveForm(A::DAForm) -> Any
 {Checks whether the input form A is positive definite on the maximal order M.}
 if assigned A`IsPositive then return A`IsPositive; end if;
 Gram:=MatrixRing(R,#IM)!0;
 for i in [1..#IM] do
  for j in [i..#IM] do
   Gram[i][j]:=InnerProduct(A`Matrix,IM[i]*Dagger(IM[j]));
   Gram[j][i]:=Gram[i][j];
  end for;
 end for;
 A`IsPositive:=IsPositiveDefinite(Gram);
 return A`IsPositive;
end intrinsic;

intrinsic IsIntegral(A::Mtrx) -> Any
 {Checks if an element of the rational algebra is contained in the maximal order specified via the basis BB.}
 V:=KMatrixSpaceWithBasis([KMatrixSpace(K,n,n)!x : x in BB]);
 c:=Eltseq(Coordinates(V,V!A));
 for i in [1..#c] do
  if not IsIntegral(c[i]) then
   return false;
  end if;
 end for;
 return true;
end intrinsic;

intrinsic IsIntegralUnit(A::Mtrx) -> Any
 {Checks if an element of the rational algebra is a unit in the maximal order specified via the basis BB.}
 return IsIntegral(A) and IsIntegral(A^-1);
end intrinsic;

intrinsic TestIsometry(A::DAForm,B::DAForm) -> Any
 {Tests whether the two forms A and B are isometric.}
 if Determinant(A`Matrix) ne Determinant(B`Matrix) then
  return false,false;
 end if;
 if #AutomorphismGroup(A) ne #AutomorphismGroup(B) then
  return false,false;
 end if;
 SA:=MinimalVectors(A); SB:=MinimalVectors(B);
 if #SA ne #SB then
  return false,false;
 end if;
 bi:=SB[1]^-1;
 for a in SA do
  x:=a*bi;
  if IsIntegralUnit(x) and B`Matrix eq Dagger(x)*A`Matrix*x then
   return true,x;
  end if;
 end for;
 return false,false;
end intrinsic;

intrinsic SymmetricCoordinates(A::Mtrx) -> Any
 {Returns the coordinates of a symmetric element w.r.t. the basis BBsym.}
 require A eq Dagger(A): "A must be dagger-symmetric";
 L:=BBsym;
 LL:=[Eltseq(L[i]) : i in [1..#L]];
 MM:=Matrix(LL);
 x:=Vector(Eltseq(A));
 a,b:=Solution(MM,x);
 return a;
end intrinsic;

intrinsic SymmetricCoordinatesToMatrix(L::SeqEnum) -> Any
 {Returns the symmetric element specified by the coordinates in the input list L.}
 return &+[L[i]*BBsym[i] : i in [1..#L]];
end intrinsic;

intrinsic IsFace(L::SeqEnum) -> Any
 {Checks the dimension of a ist of possible facet vectors.}
 return PerfectionRankDAList(L) eq dimsym-1;
end intrinsic;
