import "BasicData.m" : BB,K;

n:=NumberOfRows(BB[1]);
IM:=BB;

function Dagger(A)
 return Transpose(A);
end function;

L:=[x+Dagger(x) : x in IM];
LL:=[Eltseq(L[i]) : i in [1..#L]];
MM:=Transpose(Matrix(LL));
MM:=EchelonForm(MM);
temp:=[];
for i in [1..NumberOfRows(MM)] do
 for j in [1..NumberOfColumns(MM)] do
  if MM[i][j] eq 1 then
   Append(~temp,L[j]);
   break j;
  end if;
 end for;
end for;
BBsym:=temp;

//We compute the discriminant
F:=MatrixRing(Rationals(),#BB)!0;
for i in [1..#BB] do
 for j in [1..#BB] do
  F[i][j]:=Trace(BB[i]*BB[j]);
 end for;
end for;
MyRamifiedPrimes:=[f[1]: f in Factorization(Integers()!Determinant(F))];

phi:=RealPlaces(K)[1];
