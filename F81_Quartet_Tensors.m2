restart
--Ring definition for QUARTETS with identity at the leaves and inner egde matrix
-- on p1..p4 and EIGENVALUES l1..l2
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(5,1),l_(5,2)]
gens R

--Diagonalizing change of basis
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})

M=H*diagonalMatrix(R,4,4,{l_(5,1),l_(5,2),l_(5,2),l_(5,2)})*inverse(H)


--H^t=A^(-1), H=A^(-t), H^(-1)=A^t
--Test with different basis
AG=transpose matrix{{p_1,p_2,p_3,p_4},{-1,1,0,0},{-p_1,-p_2,p_1+p_2,0},{-p_1,-p_2,-p_3,p_1+p_2+p_3}}

inverse AG

HG=matrix {{1, 1, 1,
      1}, {(-p_2)/(p_1+p_2), p_1/(p_1+p_2), 0, 0},
      {(-p_3)/(p_1^2+2*p_1*p_2+p_2^2+p_1*p_3+p_2*p_3),
      (-p_3)/(p_1^2+2*p_1*p_2+p_2^2+p_1*p_3+p_2*p_3), 1/(p_1+p_2+p_3), 0},
      {(-p_4)/(p_1^2+2*p_1*p_2+p_2^2+2*p_1*p_3+2*p_2*p_3+p_3^2+p_1*p_4+p_2*p_4+p_3*p_4),
      (-p_4)/(p_1^2+2*p_1*p_2+p_2^2+2*p_1*p_3+2*p_2*p_3+p_3^2+p_1*p_4+p_2*p_4+p_3*p_4),
      (-p_4)/(p_1^2+2*p_1*p_2+p_2^2+2*p_1*p_3+2*p_2*p_3+p_3^2+p_1*p_4+p_2*p_4+p_3*p_4),
      1}}

M=(transpose inverse AG)*diagonalMatrix(R,4,4,{l_(5,1),l_(5,2),l_(5,2),l_(5,2)})*(transpose AG)


--Identity at the leaves
M1=id_(R^4)
M2=M1
M3=M1
M4=M1

st=toList(1..4)
S=sort elements (set st)^**4/splice/splice
netList S

--Tensor for configuration 12|34
qq=mutableMatrix(R,256,1)
for i to 255 do (
	qq_(i,0)=sum flatten toList apply(st,k->apply(st,kk->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))*M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))*M3_(position(st,l->l==kk),position(st,l->l==(S_i)_2))*M4_(position(st,l->l==kk),position(st,l->l==(S_i)_3))*M_(position(st,l->l==k),position(st,l->l==kk))))
) 
qq=matrix qq;
netList (flatten entries qq)

--Tensor with change of basis
H4=(transpose H)**(transpose H)**(transpose H)**(transpose H);
qbar=time H4*qq;
-- used 2.8125 seconds

--do not rerun unless in need of overwriting
"F81_Quartet_qbar_1234.txt" << toString qbar << endl << close


qbar_(position(S,j->j==(4,4,4,4)),0)
qbar_(position(S,j->j==(3,3,3,3)),0)
qbar_(position(S,j->j==(1,1,1,1)),0)

qbar_(position(S,j->j==(1,1,1,1)),0)*qbar_(position(S,j->j==(4,4,4,4)),0)
qbar_(position(S,j->j==(1,1,4,4)),0)*qbar_(position(S,j->j==(4,4,1,1)),0)


--DIFFERENT BASIS
--Tensor with change of basis
G4=time (inverse AG)**(inverse AG)**(inverse AG)**(inverse AG);
-- used 12.6563 seconds

qbar=time G4*qq;
-- used 2.8125 seconds

length unique select(flatten entries qbar,i->i!=0)
---------------------------------
---------------------------------
---------------------------------
--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms qbar_(position(S,j->j==i),0)>1))
netList nonMonomial
length nonMonomial --9
netList apply(nonMonomial,i->support qbar_(position(S,j->j==i),0))

nonZeroEntries=S_(positions(flatten entries qbar,i->i!=0))
length nonZeroEntries 
monomialNonZeroEntries=select(nonZeroEntries,i->(length terms qbar_(position(S,j->j==i),0)==1))
length monomialNonZeroEntries
netList monomialNonZeroEntries
---------------------------------
---------------------------------
---------------------------------

--Flattening 12|34 for identity at the leaves
flattq=mutableMatrix(R,16,16)
for i to 15 do (for j to 15 do 
    (if(j%4==j//4 and i%4==i//4) then flattq_(i,j)=p_(st_(i%4))*M_(i%4,j%4));      
    );
flattq=matrix(flattq);
flattq
--Change of basis of the flattening
flattQ=time (transpose(H)**transpose(H))*flattq*(H**H); 


--Quasi-block form (reordering according to s in next paragraph of code)
blockQ=flattQ_{0,3,12,7,13,15,5,1,4,10,2,8,6,9,11,14}^{0,3,12,7,13,15,5,1,4,10,2,8,6,9,11,14}
rank blockQ_{1,2,3,4}--1
rank blockQ_{14,15}--0
rank blockQ_{10,11,12,13}--1
rank blockQ_{7,8}--1

vflattQ = sub(blockQ, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
vflattQ

--rk 0
rank vflattQ_{position(s,i->i==(3,4)),position(s,i->i==(4,3))}--0
--rk 1
rank vflattQ_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,3))}--1
rank vflattQ_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(3,2)),position(s,i->i==(2,3))}--1
rank vflattQ_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}--1
--rk 2
rank vflattQ_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}--2
--rk 3
rank vflattQ_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}--3
rank vflattQ_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}--3


--Ring definition for flattening 12|34 general
Rgeneral=K[l_(1,1)..l_(5,4)]
gens Rgeneral
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
-- Moving from identity at the leaves to general case
flattQ=sub(blockQ,Rgeneral);

flattP=matrix toList apply(0..15,i->toList apply(0..15,j->l_(1,(s_i)_0)*l_(2,(s_i)_1)*l_(3,(s_j)_0)*l_(4,(s_j)_1)*flattQ_(i,j)));
aux=time sub(flattP,flatten toList apply(1..5,i->{l_(i,4)=>l_(i,2),l_(i,3)=>l_(i,2)}));
RG=K[l_(1,1)..l_(5,2)]
flattP=sub(aux,RG);
--rk 0
rank flattP_{position(s,i->i==(3,4)),position(s,i->i==(4,3))}--0
--rk 1
rank flattP_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,3))}--1
rank flattP_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(3,2)),position(s,i->i==(2,3))}--1
rank flattP_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}--1
--rk 2
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}--2
--rk 3
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}--3
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}--3
--rk 4: all
--column generators:
rank flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}--4
--note that it forms a diagonal submatrix with non-zero entries in the diagonal as long as pi is generic and all eigenvalues are non-zero
flattP_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}^{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(1,4))}



vflattP = sub(flattP, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
vflattP

toExternalString vflattP

--Tensor in general case
qbar=sub(qbar,Rgeneral);
qbar_(position(S,j->j==(2,2,3,3)),0)
use Rgeneral
pbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*qbar_(position(S,j->j==i),0))};
auxpbar=time sub(pbar,flatten toList apply(1..5,i->{l_(i,4)=>l_(i,2),l_(i,3)=>l_(i,2)}));
pbar=sub(auxpbar,RG);
--test
pbar_(position(S,j->j==(2,2,3,3)),0)
pbar_(position(S,j->j==(2,2,4,4)),0)
pbar_(position(S,j->j==(2,2,3,3)),0)
pbar_(position(S,j->j==(2,2,4,4)),0)
"F81_Quartet_pbar_1234.txt" << toString pbar << endl << close






