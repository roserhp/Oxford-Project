restart

--Ring declaration. Warning: for some operations it will be needed to consider a ring on QQ including parameters pi
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K; --We don't need to take eigenvalues for variables in the case of star trees for no evolution points (identity at the leaves, no interior edges)

--Matrices with identity at the leaves
M1=id_(R^4)
M2=M1
M3=M1

--Building tensor 
st={1,2,3,4}
S=sort elements (set st)^**3/splice
netList S

qq=mutableMatrix(R,64,1)
for i to 63 do (
	qq_(i,0)=sum apply({1,2,3,4},k->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))
	                           *M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))
		                   *M3_(position(st,l->l==k),position(st,l->l==(S_i)_2)))) 
qq=matrix qq;
--How to check specific entries?
qq_(position(S,i->i==(2,2,2)),0)
qq_(position(S,i->i==(2,2,3)),0)

--Basis change matrices
AG=transpose matrix{{p_1,p_2,p_3,p_4},{p_1,-p_1,0,0},{p_1,p_2,-p_1-p_2,0},{p_1,p_2,p_3,-p_1-p_2-p_3}}
AG=AG_{0}|AG_{3}|AG_{2}|AG_{1}

--Change of basis in the tensor
qbar=((inverse AG)**(inverse AG)**(inverse AG))*qq;



--Do not rerun unless you want to overwrite this file
"Tripod_qbar_NewBasis.txt" << toString qbar<< endl << close


--Different types of entries in the tensor in the new basis
qbarList=flatten entries qbar;
nonZeroEntries=S_(positions(qbarList,i->i!=0));
length nonZeroEntries --22

aux=unique select(qbarList,i->i!=0);
length aux --10
netList aux
L=time apply(aux,i->positions(qbarList,j->j==i));
 -- used 23.2344 seconds
netList apply(L,i->S_i)

zeroEntries=S_(positions(qbarList,i->i==0));
length zeroEntries --42
netList zeroEntries
