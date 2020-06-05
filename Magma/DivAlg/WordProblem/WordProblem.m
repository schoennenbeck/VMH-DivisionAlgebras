import "../../BasicData.m": K,BB,dimsym,eps;
import "../../Initialize.m": n,CP,KP,RK,TT,psi,RR,IB,IM,BBsym,Q2,RQ2K,RKPQ2,phi;

ERank:=dimsym-2;
FRank:=dimsym-1;

intrinsic ComputeBarycenters(V::VData)
 {Compute and save the attribute Barycenters in V.}
 BZ:=[];
 PL:=V`PerfectList;
 for i in [1..#PL] do
  M:=MinimalVectors(PL[i]);
  b:=1/#M * &+[(x*Dagger(x))/Trace(x*Dagger(x)) : x in M];
  Append(~BZ,b);
 end for;
 V`Barycenters:=BZ;
end intrinsic;

intrinsic ComputeEdges(V::VData)
 {Compute and save the attribute Edges in V.}
 if not assigned V`Barycenters then ComputeBarycenters(V); end if;
 Edges:=[* *];
 FL:=V`FacesList;
 PL:=V`PerfectList;
 for i in [1..#FL] do
  edges:=[* *];
  F:=FL[i];
  for j in [1..#F] do
   for k in [j+1..#F] do
    s:=SetToSequence(F[j] meet F[k]);
    M:=MinimalVectors(PL[i]);
    s:=[M[x]*Dagger(M[x]) : x in s];
    s:=[x/Trace(x) : x in s];
    d:=KMatrixSpace(K,#s,n^2)!0;
    for ii in [1..#s] do
     v:=Eltseq(s[ii]);
     for jj in [1..n^2] do
      d[ii][jj]:=v[jj];
     end for;
    end for;
    if Rank(d) eq ERank then
     Append(~edges,s);
    end if;
   end for;
  end for;
  Append(~Edges,edges);
 end for;
 V`Edges:=Edges;
end intrinsic;

intrinsic ComputeStabilizers(V::VData)
 {Compute and assign the attribute Stabilizers to V.}
 temp:=[ [x : x in AutomorphismGroup(P)] : P in V`PerfectList ];
 V`Stabilizers:=&cat temp;
end intrinsic;

intrinsic ComputePolytopeVectors(V::VData)
 {Computes a list of vectors in the symmetric space such that for each facet of the Voronoi domain of each perfect form in PerfectList the facet is given by a system of inequalities involving the scalar products on these vectors. Consider this an internal procedure.}
 if not assigned V`Edges then ComputeEdges(V); end if;
 if not assigned V`Stabilizers then ComputeStabilizers(V); end if;
 PV:=[];
 PL:=V`PerfectList;
 FL:=V`FacesList;
 EL:=V`Edges;
 for i in [1..#PL] do
  PV[i]:=[];
  Mi:=MinimalVectors(PL[i]);
  for j in [1..#FL[i]] do
   temp:=[];
   SFace:={(Mi[x]*Dagger(Mi[x]))/Trace(Mi[x]*Dagger(Mi[x])) : x in FL[i][j]};
   FB:=1/#SFace*&+SFace; //Barycenter of the face
   edgestemp:=[e : e in EL[i] | SequenceToSet(e) subset SFace];
   faceperp:=FindPerp(Setseq(SFace));
   edgestemp:=[e cat [faceperp] : e in edgestemp];
   temp:=[FindPerp(e) : e in edgestemp];
   for x in [1..#temp] do
    vz:=Sign(InnerProduct(temp[x],FB));
    temp[x]:=vz*temp[x];
   end for;
   Append(~PV[i],temp);
  end for;
 end for;
 V`PolytopeVectors:=PV;
end intrinsic;

intrinsic GeodesicIntersection(M1::Mtrx,M2::Mtrx,Face::Any,V::VData : ReturnSolution:=false) -> BoolElt
 {Checks if the geodesic from M1 to M2 intersects the face "Face" given as a pair [i,j] corresponding to the j-th face of the i-th perfect form. This function returns 0 if the geodesic does not intersect the face, 1 if it intersects the relative interior of the face and 2 if it intersects the boundary of the face.}
 M1:=M1/Trace(M1);
 M2:=M2/Trace(M2);
 i:=Face[1];
 j:=Face[2];
 M:=MinimalVectors(V`PerfectList[i]);
 PV:=V`PolytopeVectors[i][j];
 FaceGens:=[(M[x]*Dagger(M[x])/Trace(M[x]*Dagger(M[x]))) : x in V`FacesList[i][j]];
 Translations:=[FaceGens[i]-FaceGens[1] : i in [2..#FaceGens]];
 A:=KMatrixSpace(K,#FaceGens,n^2) ! 0;
 for i in [1..#Translations] do
  a:=Eltseq(Translations[i]);
  for j in [1..#a] do
   A[i][j]:=a[j];
  end for;
 end for;
 a:=Eltseq(M2-M1);
 for i in [1..#a] do
  A[#FaceGens][i]:=a[i];
 end for;
 LS:=Eltseq(M1-FaceGens[1]);
 a,b,c:=IsConsistent(A,Vector(LS));
 if not a then 
  return 0,0;
 end if;
 t:=b[#FaceGens];
 et:=-Evaluate(t,phi);
 if not (-eps lt et and et lt 1+eps) then
  return 0,0;
 end if;

 solution:=M1-t*(M2-M1);
 if solution eq M1 then 
  return 0,0; 
 end if;
 for v in PV do
  vs:=PreInnerProduct(v,solution);
  if vs eq 0 then
   return 2,solution;
  end if;
  if Evaluate(vs,phi) lt 0 then
   return 0,0;
  end if;
 end for;
 if ReturnSolution then
  return 1,solution;
 else
  return 1;
 end if;
end intrinsic;

intrinsic ConstructWord(x::Mtrx,V::VData) -> SeqEnum
 {Other method.}
 if not assigned V`Barycenters then ComputeBarycenters(V); end if;
 if not assigned V`Stabilizers then ComputeStabilizers(V); end if;
 if not assigned V`PolytopeVectors then ComputePolytopeVectors(V); end if;
 Barycenters:=V`Barycenters;
 Stabs:=V`Stabilizers;
 FL:=V`FacesList;
 FacVec:=V`FacetVectors;
 FTL:=V`FaceTrafoList;
 NL:=V`NeighbourList;
 OKGens:=V`MultFreeList;
 One:=Parent(x)!1;
 B1:=Barycenters[1];
 B:=B1;
 B2:=x^-1*B1*Transpose(x^-1);
 cur:=1; //current index
 Word:=[];
 found:=true;
 while not x in Stabs and found do
  found:=false;
  for j in [1..#FL[cur]] do
   a,b:=GeodesicIntersection(B1,B2,[cur,j],V : ReturnSolution:=true);
   if a eq 2 then
    //this is the case in which the geodesic intersects an edge of the fundamental domain
    B1/:=Trace(B1);
    B1:=(1/2)*(B1+b);
    B1pert:=B1;
    interior:=false;
    N:=1;
    while not interior do
     N*:=1/10;
     for i in [1..n] do
      randoms:=[1/Random([N,2*N]) : i in [1..Degree(K)]];
      for k in [1..Degree(K)] do
       B1pert[i][i]+:=k eq 1 select randoms[k]*K!1 else randoms[k]*(K.(k-1));
      end for;
      for j in [i+1..n] do
       randoms:=[1/Random([N,2*N]) : i in [1..Degree(K)]];
       for k in [1..Degree(K)] do
        B1pert[i][j]+:=k eq 1 select randoms[k]*K!1 else randoms[k]*(K.(k-1));
        B1pert[j][i]+:=k eq 1 select randoms[k]*K!1 else randoms[k]*(K.(k-1));
       end for;
      end for;
     end for;
     B1pert/:=Trace(B1pert);
     interior:=&and [InnerProduct(B1pert,FacVec[cur][j]) gt 0 : j in [1..#FL[cur]]];
     if not interior then
      B1pert:=B1;
     else
      B1:=B1pert;
     end if;
    end while; 
    found:=true; break j; //check all faces again
   end if;
   if a eq 1 then
    found:=true;
    g:=FTL[cur][j];
    cur:=NL[cur][j];
    x:=x*g;
    if g ne One then 
     p:=0; e:=0;
     if g in OKGens then p:=Position(OKGens,g); e:=1; end if;
     Append(~Word,[p,e]);
    end if;
    B1:=g^-1*b*Transpose(g^-1);
    B2:=g^-1*B2*Transpose(g^-1);
    break;
   end if;
  end for;
 end while;
 Append(~Word,[0,Position(Stabs,x)]);
 //At this point the output is to be read as follows:
 //x*w[1]^exp[1]*w[2]^exp[2]*....*w[r]^exp[r] = Stabs[l]
 //The stabilizer element is found in the last entry of the list "Word".

 //For convenience, we change this, so that now we have x=s*w[2]^exp[2]*...*w[r]^exp[r] and s is the stabilizer element, which is now found in the first entry
 temp:=[ Word[#Word] ];
 for i in [1..#Word-1] do
  Word[#Word-i][2]*:=-1;
  Append(~temp,Word[#Word-i]);
 end for;
 Word:=temp;
 return found,Word,x;
end intrinsic;

intrinsic ConstructWordViaBarycenters(x::Mtrx,V::VData) -> SeqEnum
 {Barycenter method - this is *not* the method described in the article. It may not terminate, but experience shows taht if it does, it produces slightly shorter words than ConstrucWord().}
 if not assigned V`Barycenters then ComputeBarycenters(V); end if;
 if not assigned V`Stabilizers then ComputeStabilizers(V); end if;
 if not assigned V`PolytopeVectors then ComputePolytopeVectors(V); end if;
 Barycenters:=V`Barycenters;
 Stabs:=V`Stabilizers;
 FL:=V`FacesList;
 FTL:=V`FaceTrafoList;
 NL:=V`NeighbourList;
 OKGens:=V`MultFreeList;
 One:=Parent(x)!1;
 B1:=Barycenters[1];
 B:=B1;
 B2:=x^-1*B1*Transpose(x^-1);
 cur:=1; 
 Word:=[];
 found:=true;
 while not x in Stabs and found do
  found:=false;
  for j in [1..#FL[cur]] do
   if GeodesicIntersection(B1,B2,[cur,j],V) in {1,2} then
    found:=true;
    g:=FTL[cur][j];
    cur:=NL[cur][j];
    x:=x*g;
    if g ne One then  
     p:=0; e:=0;
     if g in OKGens then p:=Position(OKGens,g); e:=1; end if;
     Append(~Word,[p,e]);
    end if;
    B1:=Barycenters[cur];
    B2:=x^-1*B2*Transpose(x^-1);
    break;
   end if;
  end for;
 end while;
 Append(~Word,[0,Position(Stabs,x)]);
 //At this point the output is to be read as follows:
 //x*w[1]^exp[1]*w[2]^exp[2]*....*w[r]^exp[r] = Stabs[l]
 //The stabilizer element is found in the last entry of the list "Word".

 //For convenience, we change this, so that now we have x=s*w[2]^exp[2]*...*w[r]^exp[r] and s is the stabilizer element, which is now found in the first entry
 temp:=[ Word[#Word] ];
 for i in [1..#Word-1] do
  Word[#Word-i][2]*:=-1;
  Append(~temp,Word[#Word-i]);
 end for;
 Word:=temp;
 return found,Word,x;
end intrinsic;


intrinsic SolveWordProblem(x::Mtrx,V::VData)->GrpFpElt
 {}
 found,Word,xx:=ConstructWord(x,V);
 if not found then
  error "This element does not belong to the computed group";
 end if;
 
 G,mathom:=Presentation(V);
 simpl:=V`GNonsimplifiedToSimplifiedHom;
 Gnon:=V`GNonsimplified;
 res:=G!1;
 for i in [2..#Word] do
  res:=res*simpl(Gnon.Word[i][1]^Word[i][2]);
 end for;

 s:=V`Stabilizers[Word[1][2]];
 pos:=Position([P[1]: P in V`StabilizerPreimages],s);
 res:=simpl(V`StabilizerPreimages[pos][2])*res;
 assert V`GSimplifiedHom(res) eq x;
 
 return res;
end intrinsic;
