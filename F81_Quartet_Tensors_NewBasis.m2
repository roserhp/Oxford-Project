restart
--Ring definition for QUARTETS with identity at the leaves and inner egde matrix
-- on p1..p4 and EIGENVALUES l1..l2
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(5,1),l_(5,2)]
gens R


--H^t=A^(-1), H=A^(-t), H^(-1)=A^t
--Test with different basis

AG=transpose matrix{{p_1,p_2,p_3,p_4},{p_1,-p_1,0,0},{p_1,p_2,-p_1-p_2,0},{p_1,p_2,p_3,-p_1-p_2-p_3}}
--AG=AG_{1}|AG_{2}|AG_{3}|AG_{0}
AG=AG_{0}|AG_{3}|AG_{2}|AG_{1}

inverse AG

--Careful with position of the eigenvalues!!!
M=(transpose inverse AG)*diagonalMatrix(R,4,4,{l_(5,1),l_(5,2),l_(5,2),l_(5,2)})*(transpose AG)
--M=(transpose inverse AG)*diagonalMatrix(R,4,4,{l_(5,2),l_(5,2),l_(5,2),l_(5,1)})*(transpose AG)


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

qbar_(position(S,j->j==(4,4,4,4)),0)
qbar_(position(S,j->j==(3,3,3,3)),0)
qbar_(position(S,j->j==(1,1,1,1)),0)

qbar_(position(S,j->j==(1,1,3,3)),0)==qbar_(position(S,j->j==(3,3,1,1)),0)
qbar_(position(S,j->j==(1,1,1,3)),0)
qbar_(position(S,j->j==(1,1,1,2)),0)


---------------------------------
---------------------------------
---------------------------------
--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms qbar_(position(S,j->j==i),0)>1));
netList nonMonomial
length nonMonomial --9

qbarList=flatten entries qbar;
nonZeroEntries=S_(positions(qbarList,i->i!=0));
length nonZeroEntries --112
monomialNonZeroEntries=select(nonZeroEntries,i->(length terms qbar_(position(S,j->j==i),0)==1))
length monomialNonZeroEntries
netList monomialNonZeroEntries


aux=unique select(qbarList,i->i!=0);
length aux --26
netList aux
L=time apply(aux,i->positions(qbarList,j->j==i));
 -- used 23.2344 seconds
netList apply(L,i->S_i)

zeroEntries=S_(positions(qbarList,i->i==0));
length zeroEntries --144
netList zeroEntries

member((1,2,2,3),zeroEntries)
member((1,2,3,3),zeroEntries)
member((2,3,3,3),zeroEntries)
member((4,3,3,4),zeroEntries)

qbar_(position(S,j->j==(1,1,1,2)),0)
qbar_(position(S,j->j==(1,1,1,3)),0)
qbar_(position(S,j->j==(1,1,1,4)),0)


qbar_(position(S,j->j==(1,2,2,2)),0)
qbar_(position(S,j->j==(1,3,3,3)),0)

---------------------------------
---------------------------------
---------------------------------


