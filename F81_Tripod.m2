restart

------------------------------
-- MODEL DEFINITION ---------
-----------------------------

--Ring declaration. Warning: for some operations it will be needed to consider a ring on QQ including parameters pi
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(3,2)];
gens R

--Retrieve tensor for no evolution point in the new basis
qbar=value get "Tripod_qbar.txt";

--Matrices with anything at the leaves
D1=matrix{{l_(1,1),0,0,0},{0,l_(1,2),0,0},{0,0,l_(1,2),0},{0,0,0,l_(1,2)}}
D2=matrix{{l_(2,1),0,0,0},{0,l_(2,2),0,0},{0,0,l_(2,2),0},{0,0,0,l_(2,2)}}
D3=matrix{{l_(3,1),0,0,0},{0,l_(3,2),0,0},{0,0,l_(3,2),0},{0,0,0,l_(3,2)}}


--------------------------------------
--BUILDING AND UNDERSTANDING TENSOR --
--------------------------------------

--From identity at the leaves to general case
pbar=(D1**D2**D3)*qbar;

--Index set
st={1,2,3,4}
S=sort elements (set st)^**3/splice;

--Double-checking relation between identity at the leaves and general case
qbar_(position(S,i->i==(2,2,2)),0)*l_(1,2)*l_(2,2)*l_(3,2)==pbar_(position(S,i->i==(2,2,2)),0)

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
I=time kernel f; 
betti I --check degrees of ideal generators
netList I_* --display ideal generators
