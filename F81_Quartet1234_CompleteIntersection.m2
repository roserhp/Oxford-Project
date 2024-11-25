--Compute linear equations for ptilde for the quartet 12|34
--where ptilde is obtained by rescaling pbar with the
--no evolution point whenever possible (i.e. when the corresponding coordinate is not zero)
--Note that computations are done with a random root distribution, 
--since the package MultigradedImplicitization doesn't work with fields of fractions
--Recall that we are fixing a topology, so some of this invariants will only hold for 12|34

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
pTilde=value get "F81_4leaves_tensor_ptilde.txt";
pTilde = sub(pTilde, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3});

--Define ring with random root distribution
S = QQ[l_(1,1)..l_(5,2)];
pTilde=sub(pTilde,S);
nonZeroPTilde = delete(null, apply(flatten entries pTilde, m -> if m != 0 then m));
length nonZeroPTilde

--Define ring with variables non-zero ptildes.
nonZeroEntries=value get "TN93_1234_nonZeroEntries.txt"; --list of the 80 non-zero entries for 12|34
-- Note that here they are x's
T = QQ[apply(nonZeroEntries, ind -> x_ind)];

--Define parametrization map
f = map(S,T, apply(nonZeroPTilde, i -> sub(i,S)));

---------------------------------------------------------------
--Linear equations
---------------------------------------------------------------

G1=time componentsOfKernel(1,f); -- used 9.56805 seconds
L=flatten toList apply(keys G1,i->G1#i); --list of linear equations
LT=apply(L,l->sub(l,T));
I=trim sub(ideal L,T); --ideal of linear equations
dim I,degree I
betti I
netList I_*

length select(flatten entries gens I,i->length (terms i)==2) --59
I2=ideal select(flatten entries gens I,i->length (terms i)==2);
betti (trim I2)
length select(flatten entries gens I,i->length (terms i)==3) --8
I3=ideal select(flatten entries gens I,i->length (terms i)==3) 
netList I3_*
support I3
(x_(4,4,3,3)-x_(3,3,4,4)) % gb I
(x_(4,4,2,2)-x_(2,2,4,4)) % gb I
(x_(3,3,2,2)-x_(2,2,3,3)) % gb I

I2plus=I2+ideal{x_(4,4,3,3)-x_(3,3,4,4),x_(4,4,2,2)-x_(2,2,4,4),x_(3,3,2,2)-x_(2,2,3,3)};
betti (trim I2plus) --62

Rp=K[gens T]
p12=p_1+p_2
p34=p_3+p_4
p1144=p12/(p_1*p_2)
p1122=1/(p12*p34)
p1133=p34/(p_3*p_4)
p1222=(p34-p12)/(p12^2*p34^2)
p1233=-1/(p_3*p_4)
p1244=1/(p_1*p_2)
p1333=p34*(p_4-p_3)/(p_3^2*p_4^2)
p1444=p12*(p_2-p_1)/(p_1^2*p_2^2)
p2222=(p12^2-p12*p34+p34^2)/(p12^3*p34^3)
p2233=1/(p_3*p_4*p34)
p2244=1/(p_1*p_2*p12)
p2333=(p_3-p_4)/(p_3^2*p_4^2)
p2444=(p_2-p_1)/(p_1^2*p_2^2)
p3333=(p_3^3+p_4^3)/(p_3^3*p_4^3)
p4444=(p_1^3+p_2^3)/(p_1^3*p_2^3)

T3=matrix{{p1144,p2244*x_(2,2,4,4),p1244},
          {p1122,p2222*x_(2,2,2,2),p1222},
	  {p1133,p2233*x_(2,2,3,3),p1233}};
t3=det T3;

T4=matrix{{p1144,p1244,x_(3,3,4,4),0},
          {p1122,p1222,p2233*x_(2,2,3,3),0},
	  {p1133,p1233,p3333*x_(3,3,3,3),p1333},
	  {0,0,p2333*x_(4,4,4,2),p1233}};
t4=det T4;

T5=matrix{{p1144,p1244,p4444*x_(4,4,4,4),p1444},
          {p1122,p1222,p2244*x_(2,2,4,4),0},
	  {p1133,p1233,x_(3,3,4,4),0},
	  {0,0,p2444*x_(4,4,4,2),p1244}};
t5=det T5;

T61=matrix{{p1122,p1222},{p1244*x_(4,4,4,2),p2244*x_(2,2,4,4)}}
T62=matrix{{p1122,p1222},{p1233*x_(4,4,4,2),p2233*x_(2,2,3,3)}}
t6=p1133*det T61-p1144*det T62

T71=matrix{{p1122,p1244},{p1233*x_(4,4,4,2),x_(3,3,4,4)}}
T72=matrix{{p1122,p1244},{p1244*x_(4,4,4,2),p4444*x_(4,4,4,4)}}
t7=(p1144)^2*det T71-p1133*p1144*det T72+p1133*p1444^2*p1122*x_(4,4,4,2)


m3=sub(t3, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})
m4=sub(t4, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})
m5=sub(t5, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})
m6=sub(t6, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})
m7=sub(t7, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})

m3=sub(m3,T) 
m4=sub(m4,T)  
m5=sub(m5,T) 
m6=sub(m6,T)
m7=sub(m7,T)

m4 % gb I
m6 % gb I
m7 % gb I

support m3
support m4
support m5
support m6
support m7

Im=trim ideal{m3,m4,m5,m6,m7}
betti Im

isSubset(Im,I) --true

I2plus+Im==I --true

---------------------------------------------------------------
--Binomial quadratic equations
---------------------------------------------------------------

use T
q=matrix{{x_(1,2,1,2),x_(1,2,2,1),x_(1,2,2,2)},
       {x_(2,1,1,2),x_(2,1,2,1),x_(2,1,2,2)},
       {x_(2,2,1,2),x_(2,2,2,1),x_(4,4,4,2)}}

Q=minors(2,q)
betti Q
codim Q, degree Q

JBin=ideal{x_(2,1,2,2)*x_(2,2,2,1)-x_(2,1,2,1)*x_(4,4,4,2),
x_(1,2,2,2)*x_(2,2,2,1)-x_(1,2,2,1)*x_(4,4,4,2),
x_(2,1,2,2)*x_(2,2,1,2)-x_(2,1,1,2)*x_(4,4,4,2),
x_(2,1,2,1)*x_(2,2,1,2)-x_(2,1,1,2)*x_(2,2,2,1),
x_(1,2,2,2)*x_(2,2,1,2)-x_(1,2,1,2)*x_(4,4,4,2),
x_(1,2,2,1)*x_(2,2,1,2)-x_(1,2,1,2)*x_(2,2,2,1),
x_(1,2,2,2)*x_(2,1,2,1)-x_(1,2,2,1)*x_(2,1,2,2),
x_(1,2,2,2)*x_(2,1,1,2)- x_(1,2,1,2)*x_(2,1,2,2),
x_(1,2,2,1)*x_(2,1,1,2)- x_(1,2,1,2)*x_(2,1,2,1)}
--Binomial quadrics from Jennifer's computations
betti trim JBin
codim JBin, degree JBin --(4,6)

Q==JBin --true

qlocal=ideal{det q_{0,1}^{0,1},det q_{0,2}^{0,1},det q_{0,1}^{0,2},det q_{0,2}^{0,2}}
codim qlocal, degree qlocal --(3,2) 

saturate(qlocal,ideal{q_(0,0)})==Q --true

JacQLocal=sub(jacobian qlocal,apply(support qlocal,i->i=>1)) --jacobian evaluated at no evolution point
rank JacQLocal --4

--qlocal is NOT a complete intersection but it is locally at the no evolution point

PD=time primaryDecomposition qlocal;
 -- used 3.98437 seconds
netList PD
PD_2==Q --true

-------------------------------------------------------------
-- Non-binomial quadratic equations
------------------------------------------------------------


