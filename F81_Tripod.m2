restart

--Ring declaration. Warning: for some operations it will be needed to consider a ring on QQ including parameters pi
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(3,2)];
gens R

--Retrieve tensor for no evolution point in the new basis
qbar=value get "Tripod_qbar.txt";

--Index set
st={1,2,3,4}
S=sort elements (set st)^**3/splice;

--Matrices with anything at the leaves
D1=matrix{{l_(1,1),0,0,0},{0,l_(1,2),0,0},{0,0,l_(1,2),0},{0,0,0,l_(1,2)}}
D2=matrix{{l_(2,1),0,0,0},{0,l_(2,2),0,0},{0,0,l_(2,2),0},{0,0,0,l_(2,2)}}
D3=matrix{{l_(3,1),0,0,0},{0,l_(3,2),0,0},{0,0,l_(3,2),0},{0,0,0,l_(3,2)}}

--From identity at the leaves to general case
pbar=(D1**D2**D3)*qbar;
