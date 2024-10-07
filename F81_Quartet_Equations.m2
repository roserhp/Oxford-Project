restart
needsPackage "MultigradedImplicitization";
--needsPackage "Polyhedra";
--Ring definitions: 
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,2)];
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=K[var];

--Retrieve tensor on our desired basis
use R;
--choose coordinates
pbar=value get "F81_4leaves_tensor.txt";
--pbar=value get "F81_Quartet_pbar_1234.txt";

--pbar=value get "F81_4leaves_tensor_ptilde.txt"; (parametrization as monic as possible)
--does not give nicer linear equations

F = map(R,Rx,transpose pbar);

-- substitute in a random root distribution
pbar = sub(pbar, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
nonZeroPBar = delete(null, apply(flatten entries pbar, m -> if m != 0 then m));
length nonZeroPBar 

--substitute special root distributions
--pbar= sub(pbar, {p_1 => 1/4, p_2 => 1/4, p_3 => 1/4, p_4 => 1/4});
--nonZeroPBar = delete(null, apply(flatten entries pbar, m -> if m != 0 then m));
--length nonZeroPBar

--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms pbar_(position(S,j->j==i),0)>1));
netList nonMonomial
netList apply(nonMonomial,i->pbar_(position(S,j->j==i),0))

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0));
length nonZeroEntries
netList nonZeroEntries
"TN93_1234_nonZeroEntries.txt" << toString nonZeroEntries << endl << close

monomialNonZeroEntries=select(nonZeroEntries,i->(length terms pbar_(position(S,j->j==i),0)==1));
netList monomialNonZeroEntries


--Parametrization: almost monomial map

S = QQ[l_(1,1)..l_(5,2)];
T = QQ[apply(nonZeroEntries, ind -> x_ind)];
f = map(S,T, apply(nonZeroPBar, i -> sub(i,S)));

--grading
D = maxGrading f;
D

--dimension of grading
rank(D) --5

--compute linear invariants
G1=time componentsOfKernel(1,f);
   -- used 62.3125 seconds
peek G1
L1=flatten toList apply(keys G1,i->G1#i);
length L1
netList L1
"F81_4leaves_linear.txt" << toString L << endl << close

--compute the quadratics
G2 = time componentsOfKernel(2,f);
L2=flatten toList apply(keys G2,i->G2#i);
length L2
netList L2
length select(L2,i->(degree sub(i,QQ[support i]))=={2})

H = (flatten delete({},values(G))) / (g -> sub(g, source f))

--double check in the ideal
(unique (H / f)) == {0}

G3 = time componentsOfKernel(3,f);
     -- used 85.4375 seconds
L3=flatten toList apply(keys G3,i->G3#i);
length L3
netList L3
length select(L3,i->(degree sub(i,QQ[support i]))=={3})

L3T=time apply(L3,i->sub(i,T));
-- used 107.844 seconds

I=ideal L3T;
betti I
I=time trim I;
-- used 1.3125 seconds
betti I
numgens I

G4 = time componentsOfKernel(4,f); --not yet ever computed
L4=flatten toList apply(keys G4,i->G4#i);
length L4
netList L4


"JC_quartet_cubics.txt" << toString L3T << endl << close
"JC_quartet_cubics.txt" << toString L3T << endl << close


--stick in a file
fileName = "TN93_quartet_quadrics" << "";

for g in H do (
    fileName << g << endl
);

fileName << close;

------------------------------------------
------------------------------------------
restart
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=QQ[var];
L=value get "F81_4leaves_linear.txt";
IL=ideal L;
betti IL
256-176-67 --13
netList apply(L,i->support i)
length select(L,i->length (terms i)==2) --45
netList apply(select(L,i->length (terms i)==2),j->support j)
length select(L,i->length (terms i)==3) --22
netList apply(select(L,i->length (terms i)==3),j->support j)
netList L

--------------------------------------------
--------------------------------------------
--------------------------------------------


gb(M,ChangeMatrix=>true)
