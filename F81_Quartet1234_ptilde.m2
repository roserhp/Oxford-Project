--Compute linear equations for ptilde for the quartet 12|34
--where ptilde is obtained by rescaling pbar with the
--no evolution point whenever possible (i.e. when the corresponding coordinate is not zero)
--Note that computations are done with a random root distribution, 
--since the package MultigradedImplicitization doesn't work with fields of fractions
--Recall that we are fixing a topology, so some of this invariants will only hold for 12|34

restart
needsPackage "MultigradedImplicitization";

--Retrieve values of ptilde in terms of the eigenvalues 
--for generic root distribution
K=frac(QQ[p_1..p_4]);
R=K[l_(1,1)..l_(5,2)];
pTilde=value get "F81_4leaves_tensor_ptilde.txt";

a=random(0,100)/100;
b=random(0,100)/100;
c=random(0,100)/100;
d=1-a-b-c;
r=(a,b,c,d)
sum toList r

--Choose a random root distribution
-- r=(1/9,2/63,11/21,1/3)
pTilde = sub(pTilde, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3});

S = QQ[l_(1,1)..l_(5,2)];
pTilde=sub(pTilde,S);
nonZeroPTilde = delete(null, apply(flatten entries pTilde, m -> if m != 0 then m));

--Define ring with variables non-zero ptildes. Note that here they are x's
nonZeroEntries=value get "TN93_1234_nonZeroEntries.txt"; --list of the 80 non-zero entries for 12|34
T = QQ[apply(nonZeroEntries, ind -> x_ind)];
--Note that the ordering of variables affects the resulting system of generators of the ideal
--However, it seems that the original ordering is the one that maximizes the number of binomial equations
--T = QQ[reverse(apply(nonZeroEntries, ind -> x_ind))]; 

--Define parametrization map
f = map(S,T, apply(nonZeroPTilde, i -> sub(i,S)));
--f = map(S,T,reverse(apply(nonZeroPTilde, i -> sub(i,S))));

--Compute linear equations of the kernel of the map
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


m3=sub(t3, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})
m4=sub(t4, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})
m5=sub(t5, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3})

m3=sub(m3,T) 
m4=sub(m4,T)  
m5=sub(m5,T) 

Im=ideal{m3,m4,m5}

isSubset(Im,I) --true

Tquot=T/I2plus
gens Tquot
Iquot=trim sub(I,Tquot)
betti Iquot
netList Iquot_*
m4 % gb Iquot

Iaux=Iquot+Im
Iquot==Iaux
betti (trim Iaux)

Tquotquot=Tquot/sub(Im,Tquot)
Iquotquot=trim sub(Iquot,Tquotquot);
betti Iquotquot
netList Iquotquot_*



--Equations from Niharika:
use Rp
n1=(p_3-p_4)*p_3*p_4*p12^2*x_(4,4,4,2)+
         p_3*p_4*p12^2*p34*x_(2,2,3,3)-
       p12^2*(p_3^3+p_4^3)*x_(3,3,3,3)+
 p_1*p_2*p_3^2*p_4^2*p34^2*x_(3,3,4,4)
n1T=sub(sub(n1, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3}),T)

n1T % gb I --not in the ideal!!!

n2=(p_2-p_1)*p_1*p_2*p34^2*x_(4,4,4,2)+
         p_1*p_2*p34^2*p12*x_(2,2,4,4)-
       p34^2*(p_1^3+p_2^3)*x_(4,4,4,4)+
 p_3*p_4*p_1^2*p_2^2*p12^2*x_(3,3,4,4)
n2T=sub(sub(n2, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3}),T)

n2T % gb I  --not in the ideal!

n2aux=(p_2-p_1)*p_1*p_2*p34^2*x_(2,4,4,4)+
         p_1*p_2*p34^2*p12*x_(4,4,2,2)-
       p34^2*(p_1^3+p_2^3)*x_(4,4,4,4)+
 p_3*p_4*p_1^2*p_2^2*p12^2*x_(3,3,4,4)
n2auxT=sub(sub(n2aux, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3}),T)

n2auxT % gb I --not in the ideal!!!

(n2auxT % gb I)==(n2T % gb I) --true

--Check model equations from ATR - model selection
use Rp
M1=(p12/p34)*x_(2,2,3,3)-(p34/p12)*x_(2,2,4,4)-(p12/p34-p34/p12)*x_(2,3,3,3)
M2=p_1*p_2*p_3*p_4*x_(3,3,4,4)-p12^2*x_(2,2,3,3)+p12^2*x_(2,3,3,3)
M3=(p12*p34+(p12-p34)^2)/(p12*p34)*x_(2,2,2,2)-(p_3^3+p_4^3)/(p_3*p_4*p34^2)*x_(3,3,3,3)-(p_3*p_4*(p12-p34)^2-p12^2*p_3*p_4-p12*(p_3-p_4)^2)/(p12*p34*p_3*p_4)*x_(2,3,3,3)
M4=(p12*p34+(p12-p34)^2)/(p12*p34)*x_(2,2,2,2)-(p_1^3+p_2^3)/(p_1*p_2*p12^2)*x_(4,4,4,4)-(p_1*p_2*(p12-p34)^2-p34^2*p_1*p_2-p34*(p_1-p_2)^2)/(p12*p34*p_1*p_2)*x_(2,3,3,3)
M5=(p_3^3+p_4^3)/(p_3*p_4*p34^2)*x_(3,3,3,3)-(p_1^3+p_2^3)/(p_1*p_2*p12^2)*x_(4,4,4,4)-((p12^2-p34^2)*p_1*p_2*p_3*p_4+p12*p_1*p_2*(p_3-p_4)^2-p34*p_3*p_4*(p_1-p_2)^2)/(p12*p34*p_1*p_2*p_3*p_4)*x_(2,3,3,3)

M5==M4-M3

M1T=sub(sub(M1, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
M2T=sub(sub(M2, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
M3T=sub(sub(M3, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
M4T=sub(sub(M4, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)

M1T % gb I --0
M2T % gb I --0
M3T % gb I --0
M4T % gb I --0

M=ideal{M1T,M2T,M3T,M4T}
Mquot=trim sub(M,Tquot);
betti Mquot

Rmodel=Tquot/Mquot
Imodel=trim sub(I,Rmodel) 
support Imodel            --x_4444 and x_4442 don't depend on the topology but
                          --but x_3344 vanishes or not depending on it
--MODEL EQUATION: x_3344+x_3434+x_3443=alpha x_2333+ beta x_2222
--where beta=coef(x_3344)/coef(lambda1,x_33)

use Rp
A=p12*p34/(p_1*p_2*p_3*p_4)
B=p12*p34/(p12*p34+(p12-p34)^2)
beta=sub(A/B,Rp)
alpha=sub(-A-beta*(1-B),Rp)
M5=x_(3,3,4,4)-alpha*x_(2,3,3,3)-beta*x_(2,2,2,2)
M5T=sub(sub(M5, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)

M5T % gb I --0

M=ideal{M1T,M2T,M3T,M4T,M5T}
Mquot=trim sub(M,Tquot);
betti Mquot
Rmodel=Tquot/Mquot
Imodel=trim sub(I,Rmodel) 
