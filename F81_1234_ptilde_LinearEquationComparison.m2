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

--Linear equations
G1=time componentsOfKernel(1,f); -- used 9.56805 seconds
L=flatten toList apply(keys G1,i->G1#i); --list of linear equations
LT=apply(L,l->sub(l,T));
I=trim sub(ideal L,T); --ideal of linear equations
dim I,degree I
betti I
netList I_*

--Binomial equations
I2=ideal select(flatten entries gens I,i->length (terms i)==2);
betti (trim I2)  --59
I2plus=I2+ideal{x_(4,4,3,3)-x_(3,3,4,4),x_(4,4,2,2)-x_(2,2,4,4),x_(3,3,2,2)-x_(2,2,3,3)};
betti (trim I2plus) --62

--Ring on non-zero ptildes with generic root distribution
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

--------------------------------------------------------------------------------------
-- EQUATIONS FROM RANK CONDITIONS
--------------------------------------------------------------------------------------
use Rp

--Rang 2 conditions (table 3 in ATR paper)
T3=matrix{{p1144,p2244*x_(2,2,4,4),p1244},
          {p1122,p2222*x_(2,2,2,2),p1222},
	  {p1133,p2233*x_(2,2,3,3),p1233}};
t3=det T3;

tt3=-((p_1^2+2*p_1*p_2+p_2^2-p_1*p_3-p_2*p_3+p_3^2-p_1*p_4-p_2*p_4+2*p_3*p_4+p_4^2
)/(p_1^4*p_2*p_3^4*p_4+3*p_1^3*p_2^2*p_3^4*p_4+3*p_1^2*p_2^3*p_3^4*p_4+p_1*p_2^4*p_3^4*p_4+3*p_1^4*p_2*p_3^3*p_4^2+9*p_
       1^3*p_2^2*p_3^3*p_4^2+9*p_1^2*p_2^3*p_3^3*p_4^2+3*p_1*p_2^4*p_3^3*p_4^2+3*p_1^4*p_2*p_3^2*p_4^3+9*p_1^3*p_2^2*p_3^2*p_4^3+9*p_1^2*p_2^3*p_3^2*p_4^3+3*p_1*p_2^4*p_3^2*p_4^3+p_1^4*p_2*p_3*p_4^4+3
       *p_1^3*p_2^2*p_3*p_4^4+3*p_1^2*p_2^3*p_3*p_4^4+p_1*p_2^4*p_3*p_4^4))*x_(2,2,2,2)+(1/(p_1*p_2*p_3^4*p_4+3*p_1*p_2*p_3^3*p_4^2+3*p_1*p_2*p_3^2*p_4^3+p_1*p_2*p_3*p_4^4))*x_(2,2,3,3)+(1/(p_1^4*p_2*
       p_3*p_4+3*p_1^3*p_2^2*p_3*p_4+3*p_1^2*p_2^3*p_3*p_4+p_1*p_2^4*p_3*p_4))*x_(2,2,4,4)



--Rang 3 condition (table 4 in ATR paper)
T4=matrix{{p1144,p1244,x_(3,3,4,4),0},
          {p1122,p1222,p2233*x_(2,2,3,3),0},
	  {p1133,p1233,p3333*x_(3,3,3,3),p1333},
	  {0,0,p2333*x_(4,4,4,2),p1233}};
t4=det T4;

tt4=(-1/(p_1*p_2*p_3^4*p_4^3+p_1*p_2*p_3^3*p_4^4))*x_(2,2,3,3)+((p_3^2-p_3*p_4+p_4^2)/(p_1*p_2*p_3^5*p_4^4+p_1*p_2*p_3^4*p_4^5))*x_(3,3,3,3)+(1/(p_1^2*p_3^2*p_4^2+2*p_1*p_2*p_3^2*p_
       4^2+p_2^2*p_3^2*p_4^2))*x_(3,3,4,4)+((-p_3^2+2*p_3*p_4-p_4^2)/(p_1*p_2*p_3^5*p_4^4+p_1*p_2*p_3^4*p_4^5))*x_(4,4,4,2)


--Rang 3 condition (table 5 in ATR paper)
T5=matrix{{p1144,p1244,p4444*x_(4,4,4,4),p1444},
          {p1122,p1222,p2244*x_(2,2,4,4),0},
	  {p1133,p1233,x_(3,3,4,4),0},
	  {0,0,p2444*x_(4,4,4,2),p1244}};
t5=det T5;

tt5=(1/(p_1^4*p_2^3*p_3*p_4+p_1^3*p_2^4*p_3*p_4))*x_(2,2,4,4)+((-1)/(p_1^2*p_2^2*p_3^2+2*p_1^2*p_2^2*p_3*p_4+p_1^2*p_2^2*p_4^2))*x_(3,3,4,4)+((p_1^2-2*p_1*p_2+p_2^2)/(p_1^5*p_2^4*p_
       3*p_4+p_1^4*p_2^5*p_3*p_4))*x_(4,4,4,2)+((-p_1^2+p_1*p_2-p_2^2)/(p_1^5*p_2^4*p_3*p_4+p_1^4*p_2^5*p_3*p_4))*x_(4,4,4,4)

--Equation combining ranks 1 and 2 
T61=matrix{{p1122,p1222},{p1244*x_(4,4,4,2),p2244*x_(2,2,4,4)}}
T62=matrix{{p1122,p1222},{p1233*x_(4,4,4,2),p2233*x_(2,2,3,3)}}
t6=p1133*det T61-p1144*det T62;

--Equation combining ranks 1 and 3 
T71=matrix{{p1122,p1244},{p1233*x_(4,4,4,2),x_(3,3,4,4)}}
T72=matrix{{p1122,p1244},{p1244*x_(4,4,4,2),p4444*x_(4,4,4,4)}}
T73=matrix{{p1122,p1244},{0,p1444}}
t7=(p1144)^2*det T71-p1133*p1144*det T72+p1133*p1444*det T73*x_(4,4,4,2);

--Rang equations specialised to random root distribution
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

Im=ideal{m3,m4,m5,m6,m7}
betti (trim Im)

It=ideal{tt3,tt4,tt5,t6,t7}
betti (trim It)


I2plus+Im==I --true

--------------------------------------------------------------------------------------
-- EQUATIONS FROM PARAMETRIZATION
--------------------------------------------------------------------------------------

use Rp
--Model equations that do not yield binomial equations when specializing to tree 12|34
M1=p34^2*x_(2,2,4,4)-p12^2*x_(2,2,3,3)-(p34-p12)*x_(4,4,4,2)
M4=p_1*p_2*p_3*p_4*x_(3,3,4,4)-p12^2*x_(2,2,3,3)+p12^2*x_(4,4,4,2)
M7=(p12*p34+(p12-p34)^2)/(p12*p34)*x_(2,2,2,2)-(p_3^3+p_4^3)/(p_3*p_4*p34^2)*x_(3,3,3,3)-(p_3*p_4*(p12-p34)^2-p12^2*p_3*p_4-p12*(p_3-p_4)^2)/(p12*p34*p_3*p_4)*x_(4,4,4,2)
M8=(p12*p34+(p12-p34)^2)/(p12*p34)*x_(2,2,2,2)-(p_1^3+p_2^3)/(p_1*p_2*p12^2)*x_(4,4,4,4)-(p_1*p_2*(p12-p34)^2-p34^2*p_1*p_2-p34*(p_1-p_2)^2)/(p12*p34*p_1*p_2)*x_(4,4,4,2)
M9=p_1*p_2*p_3*p_4*x_(3,3,4,4)-(p12*p34+(p12-p34)^2)*x_(2,2,2,2)+(p12*p34+(p12-p34)^2)*x_(4,4,4,2)

--Parametrization equations specialised to random root distribution
M1T=sub(sub(M1, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
M4T=sub(sub(M4, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
M7T=sub(sub(M7, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
M8T=sub(sub(M8, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
M9T=sub(sub(M9, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)

IM=ideal{M1T,M4T,M7T,M8T,M9T};
betti (trim IM)

Im==IM --true
--Rang equations and model 

netList Im_* --rang equations
netList IM_* --equations from parametrization

--------------------------------------------------------------------------------------
-- EQUATIONS FROM ANNACHIARA
--------------------------------------------------------------------------------------
use Rp

--Equations from Jupiter notebook
Eq2Rk3=-x_(2,2,4,4)+
p_1*p_2*p_3*p_4*(p_1+p_2)*((p_1+p_2)-(p_3+p_4))/((p_3+p_4)^2*((p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4)))*x_(3,3,4,4)+
(p_1^3+p_2^3)*(p_3+p_4)/((p_1+p_2)*p_1*p_2*(p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4))*x_(4,4,4,4)+
(p_1*p_2*(p_1+p_2)^2+2*p_1*p_2*(p_3+p_4)^2-(p_3+p_4)*p_1^2-(p_3+p_4)*p_2^2)/(p_1*p_2*(p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4))*x_(4,4,4,2)

Eq1Rk3 = -x_(3,3,4,4)-
p_1*p_2*p_3*p_4/((p_3+p_4)^2+(p_1+p_2)^2-(p_1+p_2)*(p_3+p_4))*x_(3,3,4,4)+
((p_1+p_2)*(p_3^3+p_4^3)/(((p_3+p_4)^2+(p_1+p_2)^2-(p_1+p_2)*(p_3+p_4))*p_3*p_4*(p_3+p_4)))*x_(3,3,3,3) +
((p_1+p_2)^2/p_1*p_2*p_3*p_4)*x_(2,2,3,3)-
x_(4,4,4,2)-
((p_1+p_2)^2/p_1*p_2*p_3*p_4)*x_(4,4,4,2)+
(p_3*p_4*(p_3+p_4)^2-(p_1+p_2)*p_3^2+2*p_3*p_4*(p_1+p_2)^2-(p_1+p_2)*p_4^2)/(p_3*p_4*((p_3+p_4)^2+(p_1+p_2)^2-(p_1+p_2)*(p_3+p_4))) *x_(4,4,4,2)

AK7=-(p_1^3+p_2^3)*(p_3+p_4)/((p_1+p_2)*p_1*p_2*((p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4)))*x_(4,4,4,4)+
(p_1*p_2*p_3*p_4)/((p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4))*x_(3,3,4,4)-
(p_1*p_2*((p_1+p_2)-(p_3+p_4))^2-(p_3+p_4)^2*p_1*p_2-(p_3+p_4)*(p_1-p_2)^2)/(p_1*p_2*((p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4)))*x_(4,4,4,2)+
x_(4,4,4,2)

--Equations from 3.6 (using chatgpt to translate to M2)
AK5= ((p_1 * p_2 * p_3 * p_4 * p12 * (p12 - p34)) / (p34^2 * (p12^2 + p34^2 - p12 * p34))) * x_(3,3,4,4) - ((p_1^3 + p_2^3) * p34) / (p12 * p_1 * p_2 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,4) + (p_1 * p_2 * p12^2 + 2 * p_1 * p_2 * p34^2 - p34 * p_1^2 - p34 * p_2^2) / (p_1 * p_2 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,2)-x_(2,2,4,4)
AK4=- ((p12^2 + p34^2 - p12 * p34 + p_1 * p_2 * p_3 * p_4) / (p12^2 + p34^2 - p12 * p34)) * x_(3,3,4,4) + (p12 * (p_3^3 + p_4^3)) / ((p12^2 + p34^2 - p12 * p34) * p_3 * p_4 * p34) * x_(3,3,3,3) + (p12^2) / (p_1 * p_2 * p_3 * p_4) * x_(2,2,3,3) + (p12 * (-p12^3 - p12 * p34^2 + p12^2 * p34 - p_1 * p_2 * p_3^2 - p_1 * p_2 * p_4^2 + p_1 * p_2 * p_3 * p_4 * p34)) / (p_1 * p_2 * p_3 * p_4 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,2)

--Corrected coeff of p2333 by Niharika (manual translation and chatgpt)
AK4=- ((p12^2 + p34^2 - p12 * p34 + p_1 * p_2 * p_3 * p_4) / (p12^2 + p34^2 - p12 * p34)) * x_(3,3,4,4) + 
(p12 * (p_3^3 + p_4^3)) / ((p12^2 + p34^2 - p12 * p34) * p_3 * p_4 * p34) * x_(3,3,3,3) + 
(p12^2) / (p_1 * p_2 * p_3 * p_4) * x_(2,2,3,3) + 
(p_3*p_4*p34^3-2*p_3*p_4*p12^2*p34-(p_3^2+p_4^2)*p12*p34-p_3*p_4*p12^2+p_3*p_4*p34^2-p_1*p_2*p_3*p12) / (p_3 * p_4 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,2)
 

--CORRECT Equations from 3.6 (using chatgpt to translate to M2) WITH A PLUS INSTEAD OF A MINUS
AK5= ((p_1 * p_2 * p_3 * p_4 * p12 * (p12 - p34)) / (p34^2 * (p12^2 + p34^2 - p12 * p34))) * x_(3,3,4,4) + ((p_1^3 + p_2^3) * p34) / (p12 * p_1 * p_2 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,4) + (p_1 * p_2 * p12^2 + 2 * p_1 * p_2 * p34^2 - p34 * p_1^2 - p34 * p_2^2) / (p_1 * p_2 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,2)-x_(2,2,4,4)

--New equations from AK
AK4=- ((p34 * p_1 * p_2 * p_3 * p_4 * (p34 - p12))) / (p12^2 * (p12^2 + p34^2 - p12 * p34)) * x_(3,3,4,4)+ 
       ((p12 * (p_3^3 + p_4^3)))/ ((p12^2 + p34^2 - p12 * p34) * p_3 * p_4 * p34) * x_(3,3,3,3)- 
       x_(2,2,3,3)+ 
       (((p_3 * p_4 * (p34^2 + 2 * p12^2) - p12 * (p_3^2 + p_4^2)))) / (p_3 * p_4 * (p12^2 + p34^2 - p12 * p34)) * x_(2,3,3,3)
--last version
AK4=((p34 * p_1 * p_2 * p_3 * p_4 * (p34 - p12)) / (p12^2 * (p12^2 + p34^2 - p12 * p34))) * x_(3,3,4,4)+ 
       ((p12 * (p_3^3 + p_4^3))/ ((p12^2 + p34^2 - p12 * p34) * p_3 * p_4 * p34)) * x_(3,3,3,3)- 
       x_(2,2,3,3)+ 
       (((p_3 * p_4 * (p34^2 + 2 * p12^2) - p12 * (p_3^2 + p_4^2))) / (p_3 * p_4 * (p12^2 + p34^2 - p12 * p34))) * x_(4,4,4,2)



ideal{Eq2Rk3}==ideal{AK5} --false
ideal{Eq1Rk3}==ideal{AK4} --false


--Parametrization equations specialised to random root distribution
AK5T=sub(sub(AK5, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
AK4T=sub(sub(AK4, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
AK7T=sub(sub(AK7, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)

AK5T % gb I -- 0
AK4T % gb I -- 0
AK7T % gb I -- 0


AK5T % gb IM -- 0
AK4T % gb IM -- 0
AK7T % gb IM -- 0

support AK5T
support AK4T
support AK7T

IAK=ideal{AK5T,AK4T,AK7T};
betti (trim IAK)

ideal{AK7T}==ideal{m7} --true
ideal{AK5T}==ideal{m5} --false
support AK5T==support m5 --true


----------------------------------------------------
-- Proposition 3.6 in main.tex
----------------------------------------------------

use Rp

Eq1=-(p_1^3+p_2^3)*(p_3+p_4)/((p_1+p_2)*p_1*p_2*((p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4)))*x_(4,4,4,4)+
(p_1*p_2*p_3*p_4)/((p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4))*x_(3,3,4,4)-
(p_1*p_2*((p_1+p_2)-(p_3+p_4))^2-(p_3+p_4)^2*p_1*p_2-(p_3+p_4)*(p_1-p_2)^2)/(p_1*p_2*((p_1+p_2)^2+(p_3+p_4)^2-(p_1+p_2)*(p_3+p_4)))*x_(4,4,4,2)+
x_(4,4,4,2)

Eq3=((p34 * p_1 * p_2 * p_3 * p_4 * (p34 - p12)) / (p12^2 * (p12^2 + p34^2 - p12 * p34))) * x_(3,3,4,4)+ 
       ((p12 * (p_3^3 + p_4^3))/ ((p12^2 + p34^2 - p12 * p34) * p_3 * p_4 * p34)) * x_(3,3,3,3)- 
       x_(2,2,3,3)+ 
       (((p_3 * p_4 * (p34^2 + 2 * p12^2) - p12 * (p_3^2 + p_4^2))) / (p_3 * p_4 * (p12^2 + p34^2 - p12 * p34))) * x_(4,4,4,2)

Eq2= ((p_1 * p_2 * p_3 * p_4 * p12 * (p12 - p34)) / (p34^2 * (p12^2 + p34^2 - p12 * p34))) * x_(3,3,4,4) + ((p_1^3 + p_2^3) * p34) / (p12 * p_1 * p_2 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,4) + (p_1 * p_2 * p12^2 + 2 * p_1 * p_2 * p34^2 - p34 * p_1^2 - p34 * p_2^2) / (p_1 * p_2 * (p12^2 + p34^2 - p12 * p34)) * x_(4,4,4,2)-x_(2,2,4,4)

Eq4=(1/p12^3)*x_(2,2,4,4) - ((p12^2 + p34^2 - p12*p34) / (p12^3 * p34^3)) * x_(2,2,2,2) + (1/p34^3)*x_(2,2,3,3)
Eq5=(p12^2/p34^2)*x_(2,2,3,3)+((p34^2-p12^2)/p34^2)*x_(4,4,4,2)-x_(2,2,4,4)


Eq1T=sub(sub(Eq1, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
Eq2T=sub(sub(Eq2, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
Eq3T=sub(sub(Eq3, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
Eq4T=sub(sub(Eq4, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)
Eq5T=sub(sub(Eq5, {p_1 => r_0, p_2 => r_1, p_3 => r_2, p_4 => r_3}),T)

Eq1T % gb I -- 0
Eq2T % gb I -- 0
Eq3T % gb I -- 0
Eq4T % gb I -- 0
Eq5T % gb I -- 0

Eq1T % gb IM -- 0
Eq2T % gb IM -- 0
Eq3T % gb IM -- 0
Eq4T % gb IM -- 0
Eq5T % gb IM -- 0

IAK=ideal{Eq1T,Eq2T,Eq3T,Eq4T,Eq5T};
betti (trim IAK)

IM==IAK --true
Im==IAK --true

ideal{m3}==ideal{Eq4T} --true: rank 2 condition
m3
Eq4T
ideal{m6}==ideal{Eq5T} --true: ranks 1 and 2 combined
m6
Eq5T
ideal{m7}==ideal{Eq1T} --true: ranks 1 and 3 combined
m7
Eq1T
ideal{m4,m5}==ideal{Eq2T,Eq3T} --false
support ideal{m4,m5}==support ideal{Eq2T,Eq3T} --true

m4
Eq3T
ideal{m4}==ideal{Eq3T} --false
support m4==support Eq3T --true
ideal{m5}==ideal{Eq2T} --false
support m5==support Eq2T  --true

ideal{m3,m6,m7}==ideal{Eq1T,Eq4T,Eq5T} --true
ideal{m3,m4,m5,m7}==ideal{Eq1T,Eq2T,Eq3T,Eq4T} --false
support ideal{m3,m4,m5,m7}==support ideal{Eq1T,Eq2T,Eq3T,Eq4T} --true
ideal{m3,m4,m5,m6,m7}==ideal{Eq1T,Eq2T,Eq3T,Eq4T,Eq5T} --true


ideal{M4T,M7T,M8T,M9T}==ideal{Eq1T,Eq2T,Eq3T,Eq4T} --false
ideal{M1T,M4T,M7T,M8T,M9T}==ideal{Eq1T,Eq2T,Eq3T,Eq4T,Eq5T} --true



use Rp
IEq=ideal{Eq1,Eq2,Eq3,Eq4,Eq5};
It==IEq

Eq1 % gb It -- 0
Eq2 % gb It -- not 0
Eq3 % gb It -- not 0
Eq4 % gb It -- 0
Eq5 % gb It -- 0


tt3 % gb IEq -- 0
tt4 % gb IEq -- not 0
tt5 % gb IEq -- not 0
t6 % gb IEq -- 0
t7 % gb IEq -- 0

ideal{tt3,t6,t7}==ideal{Eq1,Eq4,Eq5} --true

ideal{t6}==ideal{Eq5} --true
ideal{t7}==ideal{Eq1} --true
ideal{tt3}==ideal{Eq4} --false
netList terms t3
netList terms Eq4
netList apply(terms t3,i->coefficients i)


flatten entries (coefficients sub(p2444*x_(2,3,3,3),Rp))_1
numerator sub(((coefficients sub(p2444*x_(2,3,3,3),Rp))_1)_(0,0),K)

numerator (p12/p34)

help numerator

netList apply(terms t3,i->(coefficients i))
netList apply(terms Eq4,i->(coefficients i))

netList apply(terms t3,i->factor numerator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms t3,i->factor denominator sub(((coefficients i)_1)_(0,0),K))

netList apply(terms Eq4,i->factor numerator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms Eq4,i->factor denominator sub(((coefficients i)_1)_(0,0),K))

toExternalString (terms t3)_0

toExternalString factor numerator sub(((coefficients (terms t3)_0)_1)_(0,0),K)
p_1^2+2*p_1*p_2+p_2^2-p_1*p_3-p_2*p_3+p_3^2-p_1*p_4-p_2*p_4+2*p_3*p_4+p_4^2

((p_1^2+2*p_1*p_2+p_2^2-p_1*p_3-p_2*p_3+p_3^2-p_1*p_4-p_2*p_4+2*p_3*p_4+p_4^2
)/(p_1^4*p_2*p_3^4*p_4+3*p_1^3*p_2^2*p_3^4*p_4+3*p_1^2*p_2^3*p_3^4*p_4+p_1*p_2^4*p_3^4*p_4+3*p_1^4*p_2*p_3^3*p_4^2+9*p_
       1^3*p_2^2*p_3^3*p_4^2+9*p_1^2*p_2^3*p_3^3*p_4^2+3*p_1*p_2^4*p_3^3*p_4^2+3*p_1^4*p_2*p_3^2*p_4^3+9*p_1^3*p_2^2*p_3^2*p_4^3+9*p_1^2*p_2^3*p_3^2*p_4^3+3*p_1*p_2^4*p_3^2*p_4^3+p_1^4*p_2*p_3*p_4^4+3
       *p_1^3*p_2^2*p_3*p_4^4+3*p_1^2*p_2^3*p_3*p_4^4+p_1*p_2^4*p_3*p_4^4))*x_(2,2,2,2)


ideal{tt3}==ideal{Eq4} --true


netList apply(terms t4,i->factor numerator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms t4,i->factor denominator sub(((coefficients i)_1)_(0,0),K))

netList apply(terms t5,i->factor numerator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms t5,i->factor denominator sub(((coefficients i)_1)_(0,0),K))

netList apply(terms t6,i->factor numerator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms t6,i->factor denominator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms Eq5,i->factor numerator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms Eq5,i->factor denominator sub(((coefficients i)_1)_(0,0),K))

netList apply(terms t7,i->factor numerator sub(((coefficients i)_1)_(0,0),K))
netList apply(terms t7,i->factor denominator sub(((coefficients i)_1)_(0,0),K))
