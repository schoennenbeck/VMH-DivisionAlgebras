import "../../BasicData.m": C,K,BB;
import "../../Initialize.m": n,CP,KP,RK,TT,psi,RR,eps,IB,IM,BBsym,Q2,RQ2K,RKPQ2;
import "../Functions/Lattice.m": Q2sym;

intrinsic AutomorphismGroup(A::DAForm) -> Grp
 {Calculates the automorphism group of the form A.}
 if assigned A`Aut then return A`Aut; end if;
 S:=MinimalVectors(A);
 bi:=S[1]^-1;
 Gens:=[];
 for i in [1..#S] do
  a:=S[i];
  x:=a*bi;
  if IsIntegralUnit(x) and Dagger(x)*A`Matrix*x eq A`Matrix then
   Append(~Gens,x);
  end if;
 end for;
 Gens cat:= [MatrixRing(K,n)!(-1)];
 A`Aut:=sub<GL(n,K)|Gens>;
 return A`Aut;
end intrinsic;
