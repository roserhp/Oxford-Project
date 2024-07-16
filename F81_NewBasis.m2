restart
--Field declaration
K=frac(QQ[p_1,p_2,p_3,p_4]);

---------------------------------
--Choose for TRIPOD or QUARTET---
---------------------------------
n=3
R=K[l_(1,1)..l_(3,2)];
qbar=sub(value get "Tripod_qbar_NewBasis.txt",R);
---------------------------------
n=4
R=K[l_(1,1)..l_(5,2)];
qbar=value get "F81_Quartet_qbar_1234_NewBasis.txt";
---------------------------------

------------------------------------------
--Analysis of different types of entries--
------------------------------------------
st={1,2,3,4}
S=sort elements (set st)^**n/splice;
qbarList=flatten entries qbar;

--Zero entries
zeroEntries=S_(positions(qbarList,i->i==0));
length zeroEntries 
netList zeroEntries

--Non-zero entries
uniqueNonZero=unique select(qbarList,i->i!=0);
length uniqueNonZero
L=time apply(uniqueNonZero,i->positions(qbarList,j->j==i)); --takes a bit for quartets
netList apply(L,i->S_i)

--Non-monomial entries
nonMonomial=select(S,i->(length terms qbar_(position(S,j->j==i),0)>1));
netList nonMonomial
length nonMonomial 

