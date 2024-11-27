-------------------------------------------------------------------------
-- Linear equations using parametrization with new basis in ptilde
-------------------------------------------------------------------------

restart
needsPackage "MultigradedImplicitization";

--------------------------------------------------------------------------------------
-- SETUP
--------------------------------------------------------------------------------------

--Define ring with generic root distribution
K=frac(QQ[p_1..p_4]);
R=K[l_(1,1)..l_(5,2)];

--Choose a random root distribution
a=0_QQ
b=0_QQ
c=0_QQ
d=-1_QQ
while(d<=0 or a==0 or b==0 or c==0) do (
a=random(0,100)/100;
b=random(0,100)/100;
c=random(0,100)/100;
d=1-a-b-c)
r=(a,b,c,d)
sum toList r

--Parametrization of ptilde for a random root distribution
pTilde=value get "F81_Quartet_pTilde_NewBasis.txt";
pTilde = sub(pTilde, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3});

--Define ring with random root distribution
S = QQ[l_(1,1)..l_(5,2)];
pTilde=sub(pTilde,S);
nonZeroPTilde = delete(null, apply(flatten entries pTilde, m -> if m != 0 then m));
length nonZeroPTilde --112

--Define ring with variables non-zero ptildes.

nonZeroEntries=value get "F81_Quartet_NonZeroEntries_NewBasis.txt"; --list of the 112 non-zero entries for 12|34
-- Note that here they are x's
T = QQ[apply(nonZeroEntries, ind -> x_ind)];

--Define parametrization map
f = map(S,T, apply(nonZeroPTilde, i -> sub(i,S)));

--Linear equations
G1=time componentsOfKernel(1,f);
-- used 82.0937 seconds
L=flatten toList apply(keys G1,i->G1#i); --list of linear equations
LT=apply(L,l->sub(l,T));
I=trim sub(ideal L,T); --ideal of linear equations
dim I,degree I
betti I --99
netList I_*

--Binomial equations
I2=ideal select(flatten entries gens I,i->length (terms i)==2);
betti (trim I2)  --91
netList I2_*

(x_(2,2,3,3)-x_(3,3,2,2)) % gb I --0
(x_(2,2,4,4)-x_(4,4,2,2)) % gb I --0
(x_(4,4,3,3)-x_(3,3,4,4)) % gb I --0
(x_(4,4,2,2)-x_(3,3,2,2)) % gb I --0

(x_(2,2,3,3)-x_(3,3,2,2)) % gb I2 --not 0
(x_(2,2,4,4)-x_(4,4,2,2)) % gb I2 --not 0
(x_(4,4,3,3)-x_(3,3,4,4)) % gb I2 --not 0
(x_(4,4,2,2)-x_(3,3,2,2)) % gb I2 --not 0

(x_(4,4,3,3)-x_(4,4,2,2)) % gb I --not 0
(x_(4,4,3,3)-x_(3,3,2,2)) % gb I --not 0

I2plus=trim(I2+ideal{x_(4,4,2,2)-x_(3,3,2,2),x_(2,2,3,3)-x_(3,3,2,2),x_(2,2,4,4)-x_(4,4,2,2),x_(4,4,3,3)-x_(3,3,4,4)});
betti I2plus --95

TI2=T/I2plus
II2=trim sub(I,TI2);
netList II2_*
support II2
length support II2
