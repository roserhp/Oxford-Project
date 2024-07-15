restart

--Ring declaration. Warning: for some operations it will be needed to consider a ring on QQ including parameters pi
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K; --We don't need to take eigenvalues for variables in the case of star trees for no evolution points (identity at the leaves, no interior edges)


--Basis change matrices
A=transpose(1/2*matrix{{2*p_1,2*p_2,2*p_3,2*p_4},
	               {p_1*(p_3+p_4),p_2*(p_3+p_4),-p_3*(p_1+p_2),-p_4*(p_1+p_2)},
		       {p_1*(p_2+p_4),-p_2*(p_1+p_3),p_3*(p_2+p_4),-p_4*(p_1+p_3)},
		       {p_1*(p_2+p_3),-p_2*(p_1+p_4),-p_3*(p_1+p_4),p_4*(p_2+p_3)}})
H1=matrix{{1,1/p_1,1/p_1,1/p_1},{1,1/p_2,-1/p_2,-1/p_2},{1,-1/p_3,1/p_3,-1/p_3},{1,-1/p_4,-1/p_4,1/p_4}}

H=transpose(matrix{{1,1,1,1},{1/(p_1+p_2),1/(p_1+p_2),-1/(p_3+p_4),-1/(p_3+p_4)},
	{0,0,1/p_3,-1/p_4},{1/p_1,-1/p_2,0,0}})

--are the columns of A an orthogonal eigenbasis?
M=A
for i to 2 do (for j from i+1 to 3 do (print sum toList apply(0..3,k->(1/p_(k+1))*M_(k,i)*M_(k,j))))


--are the columns of A an orthonormal eigenbasis?
for i to 3 do (print sum toList apply(0..3,k->(1/p_(k+1))*M_(k,i)*M_(k,i)))

u=transpose matrix{{p_1,p_2,p_3,p_4}}
v=transpose matrix{{p_1*(p_3+p_4),p_2*(p_3+p_4),-p_3*(p_1+p_2),-p_4*(p_1+p_2)}}

numcols u!=1
numrows u

--scalar product
--Input: two column vectors
--Output: scalar product
piScalarProd=(u,v)->(
    n=numrows u;
    if(numrows v!=n) then error( "Vectors must have the same size");
    if(numcols v!=1) then error( "Input must be column vectors");
    if(numcols u!=1) then error( "Input must be column vectors");
    return sum toList apply(0..(n-1),k->(1/p_(k+1))*u_(k,0)*v_(k,0))
    )
--Orthogonal basis
--Input: a matrix where column vectors form a basis
--Output: orthogonality true or false
orthoBasis=(M)->(
    n=numrows M;
    aux=0;
    for i to n-2 do (for j from i+1 to n-1 do (aux=piScalarProd(M_{i},M_{j})));
    if aux==0 then return("true");
    if aux!=0 then return("false");
    );

piScalarProd(A_0,A_1)
orthoBasis(A)
orthoBasis(H)
orthoBasis(H1)

B=transpose matrix{{p_1,p_2,p_3,p_4},{-1,0,0,1},{-1,1,0,0},{-1,0,1,0}}
orthoBasis(B)
M=B

u_2=M_{2}-(piScalarProd(M_{1},M_{2})/piScalarProd(M_{1},M_{1}))*M_{1}
u_3=M_{3}-(piScalarProd(M_{1},M_{3})/piScalarProd(M_{1},M_{1}))*M_{1}-(piScalarProd(u_2,M_{3})/piScalarProd(u_2,u_2))*u_2

piScalarProd(M_{1},u_2)
piScalarProd(M_{1},u_3)
piScalarProd(u_2,u_3)


MM=M_{0}|M_{1}|u_2|u_3
orthoBasis(MM)

B=matrix {{p_1, p_1*p_3+p_1*p_4, p_1*p_2*p_3+p_1*p_2*p_4, 1}, {p_2, p_2*p_3+p_2*p_4, -p_1*p_2*p_3-p_1*p_2*p_4, -1}, {p_3,-p_1*p_3-p_2*p_3, p_1*p_3*p_4+p_2*p_3*p_4, -1},
      {p_4, -p_1*p_4-p_2*p_4, -p_1*p_3*p_4-p_2*p_3*p_4, 1}}

orthoBasis(B)

G=transpose matrix{{p_1,p_2,p_3,p_4},{-1,1,0,0},{-p_1,-p_2,p_1+p_2,0},{-p_1,-p_2,-p_3,1-p_4}}
orthoBasis(G)
----------------------------
R=K[x_2..x_4]

M=transpose matrix{{p_1,p_2,p_3,p_4},{p_2,p_1*x_2,p_4*x_3,p_3*x_4},{p_3,p_4*x_4,p_1*x_2,p_2*x_3},{p_4,p_3*x_3,p_2*x_4,p_1*x_2}}
M=A

I=ideal{piScalarProd(M_{0},M_{1}),
piScalarProd(M_{0},M_{2}),
piScalarProd(M_{0},M_{3}),
piScalarProd(M_{1},M_{2}),
piScalarProd(M_{1},M_{3}),
piScalarProd(M_{2},M_{3})}
netList I_*
dim I, degree I


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
H=matrix{{1,1/p_1,1/p_1,1/p_1},{1,1/p_2,-1/p_2,-1/p_2},{1,-1/p_3,1/p_3,-1/p_3},{1,-1/p_4,-1/p_4,1/p_4}}

--Change of basis in the tensor
qbar=((transpose H)**(transpose H)**(transpose H))*qq

netList flatten entries qbar

netList S_(positions(flatten entries qbar,i->i==0))
netList toList apply(0..63,i->{S_i,qbar_(i,0)})
length unique S_(positions(flatten entries qbar,i->i!=0))

netList toList apply(positions(flatten entries qbar,i->i!=0),j->{S_j,qbar_(j,0)})

".txt" << toString qbar<< endl << close
