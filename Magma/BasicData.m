R:=RealField(50);
eps:=R ! 0.0000000001;
P<x>:=PolynomialRing(Rationals());


//Q_2,3
a:=2; b:=3;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[2,0]];
BB:=[eins , (1/2)*(eins+J+I*J) , (1/2)*(eins+J-I*J) ,
(1/2)*(I+I*J) ];

/*
//Q_2,5
a:=-3; b:=10;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/6,0,1/3],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_3,5
a:=2; b:=15;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,1/2],[1/2,0,1/2,-1/2],[0,1/4,1/2,1/4]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,7
a:=-2; b:=7;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,1/2],[1/2,0,1/2,-1/2],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_3,7
a:=-1; b:=21;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,0],[0,0,0,1],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;
*/
/*
//Q_5,7
a:=-2; b:=35;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,1/2],[1/2,0,1/2,-1/2],[0,1/4,1/2,1/4]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

/*
//Q_2,11
a:=-3; b:=22;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/6,0,1/3],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,3,5,7
a:=-42; b:=5;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,0],[0,0,0,1],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,3,5,7,11,13
a:=-42; b:=5;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/134,0,26/67],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,3,5,11
a:=-165; b:=2;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/2,0,1/2],[1/2,-1/2,0,1/2],[0,0,1/2,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,3,7,11
a:=-1; b:=231;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/2,1/2,1/2],[1/2,-1/2,1/2,-1/2],[1/2,-1/2,-1/2,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,13
a:=-5; b:=26;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/10,0,1/10],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_3,11
a:=-1; b:=33;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,0],[0,0,0,1],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,17
a:=-3; b:=34;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/6,0,1/3],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;
*/

/*
//Q_19,37
a:=19; b:=37;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,0],[0,0,0,1],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

/*
//Q_2,19
a:=-19; b:=2;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/2,0,0],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;
*/

/*
//Q_3,13
a:=39; b:=2;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/4,3/4],[1/2,0,1/4,-1/4],[0,1/2,1/4,1/4]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;
*/

/*
//Q_2,23
a:=-3; b:=46;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/6,0,1/3],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_3,17
a:=51; b:=5;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/10,2/5],[0,0,0,1],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_5,11
a:=55; b:=2;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/4,3/4],[1/2,0,1/4,-1/4],[0,1/2,1/4,1/4]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_3,19
a:=-1; b:=57;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,0],[0,0,0,1],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,29
a:=-3; b:=58;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,1/6,0,1/3],[0,0,1/2,1/2],[0,0,0,1]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_2,31
a:=-2; b:=31;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,1/2],[1/2,0,1/2,-1/2],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

//Q_31,37
a:=31; b:=37;
f:=x^2-b;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
J:=Mat![[theta,0],[0,-theta]]; I:=Mat![[0,1],[a,0]];
bb:=[Mat ! 1,I,J,I*J];
bm:=MatrixRing(Rationals(),4)![[1,0,0,0],[1/2,0,1/2,0],[0,0,0,1],[0,1/2,0,1/2]];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;

/*
//Q_19,37
f:=x^2-37;
K<theta>:=NumberField(f);
dimsym:=3;
Mat:=MatrixRing(K,2);
eins:=Mat ! 1;
I:=Mat![[0,1],[19,0]]; J:=Mat![[theta,0],[0,-theta]];
BB:=[eins , (1/2)*(eins+J) , I*J ,
(1/2)*(I+I*J) ];
*/
/*
//D_2,3 (Degree 3 example from paper)
f:=x^3-3*x+1;
K<theta>:=NumberField(f);
dimsym:=6;
Mat:=MatrixRing(K,3);
eins:=Mat ! 1;
A:=Mat ! DiagonalMatrix(K,[theta,theta^2-2,(theta^2-2)^2-2]);
X:=Mat ![ 0,1,0,0,0,1,2,0,0];
bb:=[A^i*X^j : i in [0..2], j in [0..2]];
BB:=[
 (1/3)*(bb[1]+2*bb[2]+bb[3]+2*bb[4]+bb[5]+2*bb[6]+bb[7]+2*bb[8]+bb[9]) , 
 (1/3)*(bb[2]+bb[3]+2*bb[5]+2*bb[6]+bb[8]+bb[9]) ,
 bb[3] ,
 (1/3)*(bb[4]+2*bb[5]+bb[6]+bb[7]+2*bb[8]+bb[9]) ,
 bb[5] ,
 bb[6] ,
 bb[7] ,
 bb[8] ,
 bb[9]
];
*/
/*
//D_3,5
f:=x^3-3*x+1;
K<theta>:=NumberField(f);
dimsym:=6;
Mat:=MatrixRing(K,3);
eins:=Mat ! 1;
A:=Mat ! DiagonalMatrix(K,[theta,theta^2-2,(theta^2-2)^2-2]);
X:=Mat ![ 0,1,0,0,0,1,5,0,0];
bb:=[eins,A,A^2,X,A*X,A^2*X,X^2,A*X^2,A^2*X^2];
bm:=(1/3)*MatrixRing(Rationals(),9)![ 1,2,1,2,1,2,1,2,1,0,1,1,0,2,2,0,1,1,0,0,3,0,0,0,0,0,0,0,0,0,1,2,1,1,2,1,0,0,0,0,3,0,0,0,0 , 0,0,0,0,0 , 3 ,0,0,0 , 0,0,0,0,0,0,3,0,0 , 0,0,0,0,0,0,0,3,0 , 0,0,0,0,0,0,0,0,3 ];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;
*/
/*
//D_2,5
f:=x^3-3*x+1;
K<theta>:=NumberField(f);
dimsym:=6;
Mat:=MatrixRing(K,3);
eins:=Mat ! 1;
A:=Mat ! DiagonalMatrix(K,[theta,theta^2-2,(theta^2-2)^2-2]);
X:=Mat ![ 0,1,0,0,0,1,10,0,0];
bb:=[eins,A,A^2,X,A*X,A^2*X,X^2,A*X^2,A^2*X^2];
bm:=(1/9)*MatrixRing(Rationals(),9)![ 1,2,1,1,2,7,7,2,7,0,3,3,0,3,3,0,3,3,0,0,3,0,3,6,0,6,0,0,0,0,3,6,3,6,3,6,0,0,0,0,3,3,6,
0,3,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9];
BB:=[];
for i in [1..#bb] do 
   m:=Mat ! 0;
   for j in [1..#bb] do m +:= K ! bm[i][j]*bb[j]; end for;
   Append(~BB,m);
end for;
*/
/*
//D_2,7
f:=x^3+x^2-2*x-1;
K<theta>:=NumberField(f);
dimsym:=6;
Mat:=MatrixRing(K,3);
eins:=Mat ! 1;
G,_,map:=AutomorphismGroup(K);
A:=Mat ! DiagonalMatrix(K,[theta,map(G.1)(theta),map(G.1^2)(theta)]);
X:=Mat ![ 0,1,0,0,0,1,2,0,0];
bb:=[Mat!1,A,A^2,X,A*X,A^2*X,X^2,A*X^2,A^2*X^2];
BB:=bb;
*/
