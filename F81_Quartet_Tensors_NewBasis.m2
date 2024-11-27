------------------------------------------------------
--- Computation of tensors in new basis and ptilde
------------------------------------------------------

-- We can jump directly to section 4 unless we want to rerun something
 
-------------------------------------------------------------------------------------
-- 1. qbar
-------------------------------------------------------------------------------------

restart
--Ring definition for QUARTETS with identity at the leaves and inner egde matrix
-- on p1..p4 and EIGENVALUES l1..l2
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(5,1),l_(5,2)]
gens R

AG=transpose matrix{{p_1,p_2,p_3,p_4},{p_1,-p_1,0,0},{p_1,p_2,-p_1-p_2,0},{p_1,p_2,p_3,-p_1-p_2-p_3}}
AG=AG_{0}|AG_{3}|AG_{2}|AG_{1}

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

--DIFFERENT BASIS
--Tensor with change of basis
G4=time (inverse AG)**(inverse AG)**(inverse AG)**(inverse AG);
-- used 12.6563 seconds

qbar=time G4*qq;
 -- used 86.7656 seconds

--DO NOT RERUN UNLESS YOU REALLY WANT TO EDIT THE FILE
"F81_Quartet_qbar_1234_NewBasis.txt" << toString qbar << endl << close

-----------------------------------------------------------------
-- 2. pbar, p0 and ptilde
-----------------------------------------------------------------

restart
K=frac(QQ[p_1..p_4])
Rgeneral=K[l_(1,1)..l_(5,4)]
gens Rgeneral
qbar=value get "F81_Quartet_qbar_1234_NewBasis.txt";
st=toList(1..4);
S=sort elements (set st)^**4/splice/splice;
auxpbar=transpose matrix{toList apply(S,i->l_(1,i_0)*l_(2,i_1)*l_(3,i_2)*l_(4,i_3)*qbar_(position(S,j->j==i),0))};
pbar=time sub(auxpbar,flatten toList apply(1..5,i->{l_(i,4)=>l_(i,2),l_(i,3)=>l_(i,2)}));
length select(flatten entries pbar,i->i!=0)--112

--checks
pbar_(position(S,j->j==(1,2,3,3)),0)==l_(1,1)*l_(2,2)*l_(3,2)*l_(4,2)*qbar_(position(S,j->j==(1,2,3,3)),0)
pbar_(position(S,j->j==(2,2,3,3)),0)==l_(1,2)*l_(2,2)*l_(3,2)*l_(4,2)*qbar_(position(S,j->j==(2,2,3,3)),0)

--DO NOT RERUN
"F81_Quartet_pbar_1234_NewBasis.txt" << toString pbar << endl << close


--No evolution point
p0=transpose matrix{toList apply(flatten entries pbar,i->sub(i,matrix{toList apply(1..20, i->1_K)}))};
length select(flatten entries p0,i->i!=0)--112

pbar_(position(S,j->j==(2,3,2,3)),0)
p0_(position(S,j->j==(2,3,2,3)),0)

--DO NOT RERUN
"F81_Quartet_p0_NewBasis.txt" << toString p0 << endl << close

--pTilde
pTilde=transpose matrix{toList apply(0..255,i->if(p0_(i,0)!=0) then (1/p0_(i,0))*pbar_(i,0) else pbar_(i,0))};
length select(flatten entries pTilde,i->i!=0)--112

pbar_(position(S,j->j==(2,3,2,3)),0)
pTilde_(position(S,j->j==(2,3,2,3)),0)

pbar_(position(S,j->j==(2,2,3,3)),0)
pTilde_(position(S,j->j==(2,2,3,3)),0)

pbar_(position(S,j->j==(3,3,2,2)),0)
pTilde_(position(S,j->j==(3,3,2,2)),0)

pbar_(position(S,j->j==(2,2,3,3)),0)==pbar_(position(S,j->j==(3,3,2,2)),0)--true
pTilde_(position(S,j->j==(2,2,3,3)),0)==pTilde_(position(S,j->j==(3,3,2,2)),0)--true


--DO NOT RERUN
"F81_Quartet_pTilde_NewBasis.txt" << toString pTilde << endl << close

-----------------------------------------------------------------
-- 3. Non-zero entries
-----------------------------------------------------------------

nonZeroEntries=S_(positions(flatten entries p0,i->i!=0))
length nonZeroEntries  --112 

nonMonomial=select(S,i->(length terms pTilde_(position(S,j->j==i),0)>1))
length nonMonomial   --9
netList nonMonomial

--DO NOT RERUN
"F81_Quartet_NonZeroEntries_NewBasis.txt" << toString nonZeroEntries << endl << close

-----------------------------------------------------------------
-- 4. Parametrization with ptilde
-----------------------------------------------------------------

restart
K=frac(QQ[p_1..p_4])
R=K[l_(1,1)..l_(5,2)]
pTilde=value get "F81_Quartet_pTilde_NewBasis.txt";
st=toList(1..4);
S=sort elements (set st)^**4/splice/splice;

pTilde_(position(S,j->j==(3,3,2,2)),0)
pTilde_(position(S,j->j==(4,4,2,2)),0)
pTilde_(position(S,j->j==(3,3,4,4)),0)


pTilde_(position(S,j->j==(2,2,2,2)),0)
pTilde_(position(S,j->j==(3,3,3,3)),0)
pTilde_(position(S,j->j==(4,4,4,4)),0)
