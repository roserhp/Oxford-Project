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
H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})
A=inverse(transpose H) --this is the one in the paper (not exactly the same matrix because we are using p1+p2+p3+p4=1)

--Test with different basis
B=matrix {{p_1, p_1*p_3+p_1*p_4, p_1*p_2*p_3+p_1*p_2*p_4, 1}, {p_2, p_2*p_3+p_2*p_4, -p_1*p_2*p_3-p_1*p_2*p_4, -1}, {p_3,-p_1*p_3-p_2*p_3, p_1*p_3*p_4+p_2*p_3*p_4, -1},
      {p_4, -p_1*p_4-p_2*p_4, -p_1*p_3*p_4-p_2*p_3*p_4, 1}}
G=transpose matrix{{p_1,p_2,p_3,p_4},{-1,1,0,0},{-p_1,-p_2,p_1+p_2,0},{-p_1,-p_2,-p_3,1-p_4}}

inv=inverse G;
qbar=time (inv**inv**inv)*qq;
netList S_(positions(flatten entries qbar,i->i==0))
length S_(positions(flatten entries qbar,i->i!=0))
length S_(positions(flatten entries qbar,i->i==0))

length unique select(flatten entries qbar,i->i!=0)
netList unique select(flatten entries qbar,i->i!=0)


--Change of basis in the tensor
qbar=((transpose H)**(transpose H)**(transpose H))*qq


--Do not rerun unless you want to overwrite this file
"Tripod_qbar.txt" << toString qbar<< endl << close
