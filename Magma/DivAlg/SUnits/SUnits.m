import "../../BasicData.m": BB,K,R,dimsym,eps;
import "../../Initialize.m": n,IM,BBsym,phi,MyRamifiedPrimes;

VV:=KMatrixSpaceWithBasis([KMatrixSpace(K,n,n)!x : x in BB]);  //The regular (or equivalently unique simple) module
coords:=func<x| Vector([Rationals()!y: y in Eltseq(Coordinates(VV,VV!x))])>;

intrinsic OrbitAndStabilizer(start::Any,Gens::SeqEnum,act::Any)->SeqEnum,SeqEnum
 {Computes the orbit of start under the action of the group generated by Gens via act}
 Orbit:=[start];
 StabWords:=[];
 OrbitPlus:=[<start,[0]>];
 ToDo:=OrbitPlus;
 while  #ToDo gt 0 do
  t:=ToDo[1];
  Remove(~ToDo,1);
  for i in [1..#Gens] do
   tg:=act(t[1],Gens[i]);
   if tg in Orbit then
    pos:=Position(Orbit,tg);
    Append(~StabWords,t[2] cat [i] cat [-y: y in Reverse(OrbitPlus[pos][2])]);
   else
    Append(~Orbit,tg);
    Append(~OrbitPlus,<tg,t[2] cat [i]>);
    Append(~ToDo,<tg,t[2] cat [i]>);
   end if;
  end for;
 end while;
 return OrbitPlus,[[x: x in y | x ne 0]: y in StabWords];
end intrinsic;

intrinsic LatAct(L::Lat, g::Mtrx)->Lat
 {The lattice Lg where L lives in the regular module VV}
 return ZLatticeFromGenerators([v*g: v in ZLatticeToGenerators(L)]);
end intrinsic;
 

intrinsic SolveNormEquations(S::SeqEnum) -> SeqEnum
 {Computes elements of norm p (for p in S) in the maximal order corresponding to V}
 for i in [1..#S] do
  if S[i] lt 0 then
   S[i]:=-S[i];
   if S[i] in MyRamifiedPrimes or not IsPrime(S[i]) then
    error "This only makes sense for non-ramified primes.";
   end if;
  end if;
 end for;
 //We do something very stupid at this point:
 max:=1;
 elts:=[Parent(BB[1])!0: i in [1..#S]];
 counter:=1;
 while not [p: p in S] eq [Determinant(x): x in elts] do
  if counter mod 1000 eq 0 then
   max+:=1;
  end if;
  v:=&+[Random([-max..max])*BB[i]: i in [1..#BB]];
  if Determinant(v) in S then
   pos:=Position(S,Determinant(v));
   if Determinant(elts[pos]) eq 0 then
    elts[pos]:=v;
   end if;
  end if;
  counter+:=1;
 end while;
 return elts;
end intrinsic;

intrinsic ZLatticeFromGenerators(L::SeqEnum)->Lat
 {The Z-lattice generated by the elements in L as a subspace of the regular module VV}
 return LatticeWithBasis(Matrix([coords(v): v in L]));
end intrinsic;

intrinsic ZLatticeToGenerators(L::Lat) ->SeqEnum
 {The inverse}
 return [&+[v[i]*BB[i]: i in [1..#BB]]: v in Basis(L)];
end intrinsic;


intrinsic SUnits(V::VData, S::SeqEnum:projective:=false,rewriting:=false, Generators:=[]) ->Any
 {Computes a presentation for the SUnit-group of the maximal order described in V at the set of (rational) primes S}
 if n gt 3 then
  error "Not implemented for division algebras of this index yet.";
 end if;
 Gens:=SolveNormEquations(S);
 for i in [1..#Gens] do
  pos:=Position([Determinant(x): x in Generators], Determinant(Gens[i]));
  if pos ne 0 then
   Gens[i]:=Generators[pos];
  end if;
 end for;
 Stab,stab:=Presentation(V);
 if rewriting then
  StabRW:=RWSGroup(Stab);
 else
  StabRW:=Stab;
 end if;
 StabGens:=[stab(Stab.i): i in [1..Ngens(Stab)]];
 Gens:=StabGens cat SolveNormEquations(S);
 G:=FreeGroup(#Gens);
 Rel:=[G!Eltseq(LHS(r))=G!Eltseq(RHS(r)): r in Relations(Stab)];

 solver:=func<x|G!Eltseq(Stab!Eltseq(StabRW!Eltseq(SolveWordProblem(x,V))))>;
 StabCenterGens:=[MatrixRing(K,n)!p: p in S];
 StabCenterGensWords:=[];

 L0:=ZLatticeFromGenerators(BB);


 //If n=2 we have edges of minus-type and need to find the correct transforming elements.
 if n eq 2 then 
  for i in [1..#S] do
  //We compute the orbit of the stabilizer on Lp to find s1 such that Lp*s1 supset pL0 gp^-1
   p:=S[i];
   gp:=Gens[Ngens(Stab)+i];
   Lp:=LatAct(L0,gp);
   Orbit:=[Lp];
   ToDo:=[<Lp,[0]>];
   target:=p*LatAct(L0,gp^-1);
   found:=target subset Lp;
   s1:=Parent(gp)!1;
   s1word:=G!1;
   while not found do
    t:=ToDo[1];
    Remove(~ToDo,1);
    for j in [1..#StabGens] do
     tg:=LatAct(t[1],StabGens[j]);
     if target subset tg then
      s1word:=G![y: y in t[2] | y ne 0]*G.j;
      found:=true;
      s1:=stab(Stab!([y: y in t[2] | y ne 0] cat [j]));
      break j;
     elif tg in Orbit then
      continue j;
     else
      Append(~Orbit,tg);
      Append(~ToDo,<tg,t[2] cat [j]>);
     end if;
    end for;
   end while;
   Gens[Ngens(Stab)+i]:=gp*s1;
   gp:=gp*s1;
   zpword:=G.(Ngens(Stab)+i)^2*solver((gp*gp)^-1*StabCenterGens[i]);
   Append(~StabCenterGensWords,zpword);
   for j in [1..#Gens] do
    Append(~Rel,zpword*G.j*zpword^-1*G.j^-1 = G!1);
   end for;
  end for;
 end if;

 for i in [1..#S] do
  //We start with the relation coming from the single edge.
  p:=S[i];
  print "Now computing at prime ", p;
  gp:=Gens[Ngens(Stab)+i];
  if n eq 2 then
   zpword:=StabCenterGensWords[i];
   zp:=StabCenterGens[i];
  end if;
  Lp:=LatAct(L0,gp);
  Type1Neighbours,EdgeStab:=OrbitAndStabilizer(Lp,StabGens,LatAct);
  print "Edge orbit computed. There are ", #EdgeStab, " generators to consider";
  if n eq 2 then //This is the minus-type situation and we have already chosen the "correct" generator
   for w in EdgeStab cat [[]] do
    ww:=stab(w);
    Append(~Rel,zpword^-1*G.(Ngens(Stab)+i)*G!w*G.(Ngens(Stab)+i) = solver(zp^-1*gp*ww*gp));
   end for;
  else
   for w in EdgeStab do
    ww:=stab(w);
    Append(~Rel,G.(Ngens(Stab)+i)*G!w*G.(Ngens(Stab)+i)^-1 = solver(gp*ww*gp^-1));
   end for;
  end if;
  print "Edge relations done.";

  if n eq 3 then 
   //We now find an expression for p*Id as a "power" of gp. We find the corresponding cycle along the way.
   //First we compute the orbit of the stabilizer on Lp to find s1 such that Lp*s1 supset pL0 gp^-1
   Orbit:=[Lp];
   ToDo:=[<Lp,[0]>];
   target:=p*LatAct(L0,gp^-1);
   found:=target subset Lp;
   s1:=Parent(gp)!1;
   s1word:=G!1;
   while not found do
    t:=ToDo[1];
    Remove(~ToDo,1);
    for j in [1..#StabGens] do
     tg:=LatAct(t[1],StabGens[j]);
     if target subset tg then
      s1word:=G![y: y in t[2] | y ne 0]*G.j;
      found:=true;
      s1:=stab(Stab!([y: y in t[2] | y ne 0] cat [j]));
      break j;
     elif tg in Orbit then
      continue j;
     else
      Append(~Orbit,tg);
      Append(~ToDo,<tg,t[2] cat [j]>);
     end if;
    end for;
   end while;
   print "Second edge in 2-cell found.";
   //We have to do the same again with the target this time being pL0gp^-1s1^-1gp^-1
   Orbit:=[Lp];
   ToDo:=[<Lp,[0]>];
   target:=p*LatAct(L0,gp^-1*s1^-1*gp^-1);
   found:=target subset Lp;
   s2:=Parent(gp)!1;
   s2word:=G!1;
   while not found do
    t:=ToDo[1];
     Remove(~ToDo,1);
    for j in [1..#StabGens] do
     tg:=LatAct(t[1],StabGens[j]);
     if target subset tg then
      s2word:=G![y: y in t[2] | y ne 0]*G.j;
      found:=true;
      s2:=stab(Stab!([y: y in t[2] | y ne 0] cat [j]));
      break j;
      elif tg in Orbit then
      continue j;
     else
      Append(~Orbit,tg);
      Append(~ToDo,<tg,t[2] cat [j]>);
     end if;
    end for;
   end while;
   print "Third edge in 2-cell found.";
   zpword:=G.(Ngens(Stab)+i)*s2word*G.(Ngens(Stab)+i)*s1word*G.(Ngens(Stab)+i)*solver((gp*s2*gp*s1*gp)^-1*StabCenterGens[i]);
   Append(~StabCenterGensWords,zpword);
  //zpword is now a word for the central element p
   for j in [1..#Gens] do
    Append(~Rel,zpword*G.j*zpword^-1*G.j^-1 = G!1);
   end for;
   print "Stabilizer center added.";
  end if; 
  //Now L0gp s2 gp s1 gp = pL0 

  //Now we add some additional stabilizer relations (or equivalently the relations coming from the mixed 2-cells in the product of buildings.
  for j in [1..i-1] do
   gq:=Gens[j+Ngens(Stab)];
   Lpq:=LatAct(L0,gp*gq);
   Lqp:=LatAct(L0,gq*gp);
   Orbit:=[Lpq];
   ToDo:=[<Lpq,[0]>];
   target:=Lqp;
   found:=target eq Lpq;
   s:=Parent(gp)!1;
   sword:=G!1;
   while not found do
    t:=ToDo[1];
    Remove(~ToDo,1);
    for k in [1..#StabGens] do
     if #Orbit mod 200 eq 0 then print #Orbit; end if; //DEBUG PRINT
     tg:=LatAct(t[1],StabGens[k]);
     if target eq tg then
      sword:=G![y: y in t[2] | y ne 0]*G.k;
      found:=true;
      s:=stab(Stab!([y: y in t[2] | y ne 0] cat [k]));
      break k;
      elif tg in Orbit then
      continue k;
     else
      Append(~Orbit,tg);
      Append(~ToDo,<tg,t[2] cat [k]>);
     end if;
    end for;
   end while;
   assert found;
   print "Mixed relation found.";
   commrel:=gp*gq*s*gp^-1*gq^-1;
   Append(~Rel,G.(i+Ngens(Stab))*G.(j+Ngens(Stab))*sword*G.(i+Ngens(Stab))^-1*G.(j+Ngens(Stab))^-1 = solver(commrel));
  end for;
  print "Mixed relations done.";
 end for;
 assert IsSatisfied(Rel,[GL(n,K)!x: x in Gens]);
 G:=quo<G|Rel>;
 h:=hom<G->GL(n,K)|Gens>;
 if not projective then
  return G,h,_;
 else 
  return G,h,StabCenterGensWords cat [solver(MatrixRing(K,n)!-1)];
 end if;
end intrinsic;


intrinsic MyTest(V::VData, S::SeqEnum:projective:=false) ->Any
 {Computes a presentation for the SUnit-group of the maximal order described in V at the set of (rational) primes S}
 if n ge 3 or #S gt 1 then
  error "Not implemented for division algebras of this index yet.";
 end if;
 Gens:=SolveNormEquations(S);
 Stab,stab:=Presentation(V);
 StabGens:=[stab(Stab.i): i in [1..Ngens(Stab)]];
 SNE:=SolveNormEquations(S);
 ss:=SNE[1];
 Gens:=StabGens cat [ss^-1*g*ss: g in StabGens];
 StabGens2:=[ss^-1*g*ss: g in StabGens];
 G:=FreeGroup(#Gens);
 Rel:=[G!Eltseq(LHS(r))=G!Eltseq(RHS(r)): r in Relations(Stab)] cat [G![i gt 0 select i + Ngens(Stab) else i - Ngens(Stab): i in Eltseq(LHS(r))] = G![i gt 0 select i + Ngens(Stab) else i - Ngens(Stab): i in Eltseq(RHS(r))]: r in Relations(Stab)];

 solver:=func<x|G!Eltseq(SolveWordProblem(x,V))>;

 L0:=ZLatticeFromGenerators(BB);


 //If n=2 we have edges of minus-type and need to find the correct transforming elements.


 for i in [1..#S] do
  //We start with the relation coming from the single edge.
  p:=S[i];
  print "Now computing at prime ", p;
  gp:=ss;
  Lp:=LatAct(L0,gp);
  Type1Neighbours,EdgeStab:=OrbitAndStabilizer(L0,[ss^-1*g*ss: g in StabGens],LatAct);
  print "Edge orbit computed. There are ", #EdgeStab, " generators to consider";
  if n eq 2 then //This is the minus-type situation and we have already chosen the "correct" generator
   for w in EdgeStab cat [[]] do
    ww:=ss^-1*stab(w)*ss;
    assert LatAct(L0,ww) eq L0 and LatAct(Lp,ww) eq Lp; 
    www:=G![j gt 0 select j + Ngens(Stab) else j - Ngens(Stab): j in w];
    Append(~Rel,www = solver(ww));
   end for;
  
  end if;
  print "Edge relations done.";

  
 
 end for;
 assert IsSatisfied(Rel,[GL(n,K)!x: x in Gens]);
 G:=quo<G|Rel>;
 g:=hom<G->GL(n,K)|Gens>;
 H,h:=Simplify(G);
 
 hh:=hom<H->GL(n,K)|[g(G!H.i): i in [1..Ngens(H)]]>;
 if not projective then
  return H,hh,_;
 else 
  return H,hh,[H!(G!Eltseq(solver(MatrixRing(K,n)!-1)))];
 end if;
end intrinsic;