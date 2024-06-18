------------------------------
-- MODEL DEFINITION ----------
-----------------------------

--JC69

restart
R=QQ[b_1,b_2,b_3];
m1=mutableMatrix{{0,b_1,b_1,b_1},{b_1,0,b_1,b_1},{b_1,b_1,0,b_1},{b_1,b_1,b_1,0}}
apply(0..3,i->m1_(i,i)=1-sum flatten entries (matrix m1)^{i});
M1=matrix m1

m2=mutableMatrix{{0,b_2,b_2,b_2},{b_2,0,b_2,b_2},{b_2,b_2,0,b_2},{b_2,b_2,b_2,0}}
apply(0..3,i->m2_(i,i)=1-sum flatten entries (matrix m2)^{i});
M2=matrix m2

m3=mutableMatrix{{0,b_3,b_3,b_3},{b_3,0,b_3,b_3},{b_3,b_3,0,b_3},{b_3,b_3,b_3,0}}
apply(0..3,i->m3_(i,i)=1-sum flatten entries (matrix m3)^{i});
M3=matrix m3

---------------------------------
--BUILDING TENSOR----------------
---------------------------------

--Notation: 1-A, 2-C, 3-G, 4-T
st={1,2,3,4}
S=sort elements (set st)^**3/splice

pp=mutableMatrix(R,64,1)
for i to 63 do (
	pp_(i,0)=sum apply({1,2,3,4},k->(1/4)*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))
	                           *M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))
		                   *M3_(position(st,l->l==k),position(st,l->l==(S_i)_2)))) 
p=matrix pp;
--Display entries of tensor p
netList flatten entries p
--How to check specific entries?
p_(position(S,i->i==(2,2,2)),0)
p_(position(S,i->i==(2,2,3)),0)

--Fourier change of basis: Hadamard matrix
H=sub(matrix{{1,1,1,1},{1,-1,1,-1},{1,1,-1,-1},{1,-1,-1,1}},R)
(inverse H)*M1*H
(inverse H)*M2*H
(inverse H)*M3*H

pbar=((transpose H)**(transpose H)**(transpose H))*p;
--Display entries of tensor p with new coordinates
netList flatten entries pbar

--Check specific entries. Note that after the change of vars, 
-- 1,2,3,4 do not correspond to nucleotides
pbar_(position(S,i->i==(2,2,2)),0)
pbar_(position(S,i->i==(2,2,3)),0)


--------------------------------------------------------------
--VANISHING IDEAL---------------------------------------------
--------------------------------------------------------------

--Build ring with p variables (19 nonvanishing tensor coordinates)
nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0))
PBAR=flatten entries pbar^(positions(flatten entries pbar,i->i!=0));
netList PBAR -- check nonzero entries
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=QQ[varp]; 
gens Rp -- check ring variables
--Compute vanishing ideal as the kernel of the map
f=map(R,Rp,PBAR);
I=time trim kernel f; -- used 0.00760187 seconds
betti I --check degrees of ideal generators: 12 linear, 1 cubic
netList I_* --display ideal generators

