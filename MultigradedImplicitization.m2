newPackage(
  "MultigradedImplicitization",
  Version => "0.1",
  Date => "October 26, 2023",
  Authors => {
    {Name => "Joseph Cummings",
    Email => "jcummin7@nd.edu",
    HomePage => "https://sites.google.com/view/josephcummingsuky/home"},
    {Name => "Benjamin Hollering",
    Email => "benhollering@gmail.com",
    HomePage => "https://sites.google.com/view/benhollering"}
  },
  Headline => "A package for levaraging multigradings to solve implicitization problems",
  DebuggingMode => true,
  PackageImports => {"gfanInterface"}
)

--------------------
--Exports
--------------------

export {
  -- Methods
  "maxGrading",
  "findBasisInDegree",
  "componentOfKernel",
  "componentsOfKernel",
  -- Options
  "Grading", "PreviousGens", "ReturnTargetGrading"
}



---------------------
---- maxGrading -----
---------------------
maxGrading = method(Options => {ReturnTargetGrading => false});
maxGrading RingMap := Matrix => opts -> F -> (

  dom := source F;
  codom := target F;
  elimRing := dom ** codom;
  X := vars dom;
  n := numgens dom;
  elimIdeal := ideal(sub(X, elimRing) - sub(F(X), elimRing));
  
  if opts.ReturnTargetGrading then (transpose linealitySpace(gfanHomogeneitySpace(elimIdeal))) else (transpose linealitySpace(gfanHomogeneitySpace(elimIdeal)))_(toList(0..n-1))
  )


TEST ///
A = matrix {{1,1,1,0,0,0,0,0,0}, {0,0,0,1,1,1,0,0,0}, {0,0,0,0,0,0,1,1,1}, {1,0,0,1,0,0,1,0,0}, {0,1,0,0,1,0,0,1,0}};
R = QQ[x_1..x_(numcols A)];
S = QQ[t_1..t_(numrows A)];
F = map(S, R, apply(numcols(A), i -> S_(flatten entries A_i)));
assert(ker(A) == ker(maxGrading(F)));
///



----------------------------
----- findBasisInDegree ----
----------------------------
findBasisInDegree = method();
findBasisInDegree (List, Ring, List, MutableHashTable) := Matrix => (deg, dom, G, basisHash) -> (

  if #G == 0 then (
      return basis(deg, dom);
  );

  -- otherwise, we shift G in all possible ways to land in R_deg

  L := apply(G, g -> (
          checkDegree := deg - degree(g);
          if basisHash#?checkDegree then (
              g*basisHash#checkDegree
          ) else (
              g*basis(checkDegree, dom)
          )
      )
  );
  
  -- stick em all in a matrix
  mat := L#0;
  scan(1..#L-1, i -> mat = mat | L#i);

  -- and collect coefficients.
  (mons, coeffs) := coefficients(mat);

  -- find the independent linear relations
  coeffs = mingens(image sub(coeffs, QQ));

  -- remove monomials corresponding to pivots
  badMonomials := apply(pivots coeffs, i -> mons_(0,i#0));
  monomialBasis := flatten entries basis(deg, dom);
  scan(badMonomials, m -> monomialBasis = delete(m, monomialBasis));

  matrix{monomialBasis}
  )

findBasisInDegree (List, Ring, MutableHashTable) := Matrix => (deg, dom, basisHash) -> findBasisInDegree(deg, dom, {}, basisHash)

TEST ///
A = matrix {{1,1,1,0,0,0,0,0,0}, {0,0,0,1,1,1,0,0,0}, {0,0,0,0,0,0,1,1,1}, {1,0,0,1,0,0,1,0,0}, {0,1,0,0,1,0,0,1,0}};
R = QQ[x_1..x_(numcols A)];
S = QQ[t_1..t_(numrows A)];
F = map(S, R, apply(numcols(A), i -> S_(flatten entries A_i)));
dom = newRing(R, Degrees => A);
basisHash = new MutableHashTable from apply(gens(dom), i -> degree(i) => i);
assert(findBasisInDegree({1,1,1,1,0}, dom, {}, basisHash) == matrix {{x_(1,1)*x_(2,2), x_(1,2)*x_(2,1)}});
B = basis(2, source F);
lats = unique apply(flatten entries B, i -> degree(sub(i, dom)));
scan(lats, deg -> basisHash#deg = findBasisInDegree(deg, dom, {}, basisHash));
assert(findBasisInDegree({2,1,1,1,1}, deg,  {x_(1,2)*x_(2,1)-x_(1,1)*x_(2,2), x_(1,3)*x_(2,1)-x_(1,1)*x_(2,3), x_(1,3)*x_(2,2)-x_(1,2)*x_(2,3)}, basisHash) == matrix {{x_(1,2)*x_(1,3)*x_(2,1)}});
///


----------------------------
----- componentOfKernel ----
----------------------------
componentOfKernel = method(Options => {PreviousGens => {}});
componentOfKernel (List, Ring, RingMap, Matrix) := List => opts -> (deg, dom, F, monomialBasis) -> (


  -- collect coefficients into a matrix
  (mons, coeffs) := coefficients(F(sub(monomialBasis, source F)));

  -- find the linear relations among coefficients
  K := gens ker sub(coeffs,QQ);

  newGens := flatten entries (monomialBasis * K);

  
  newGens
  )


componentOfKernel (List, Ring, RingMap, MutableHashTable) := List => opts ->  (deg, dom, F, basisHash) -> (

      monomialBasis := if basisHash#?deg then basisHash#deg else findBasisInDegree(deg, dom, opts.PreviousGens, basisHash);

      componentOfKernel(deg, dom, F, monomialBasis)   
  )


componentOfKernel (List, Ring, RingMap) := List => opts -> (deg, dom, F) -> (


  monomialBasis := basis(deg, dom);

  componentOfKernel(deg, dom, F, monomialBasis)
  )



TEST ///
A = matrix {{1,1,1,0,0,0,0,0,0}, {0,0,0,1,1,1,0,0,0}, {0,0,0,0,0,0,1,1,1}, {1,0,0,1,0,0,1,0,0}, {0,1,0,0,1,0,0,1,0}};
R = QQ[x_1..x_(numcols A)];
S = QQ[t_1..t_(numrows A)];
F = map(S, R, apply(numcols(A), i -> S_(flatten entries A_i)));
dom = newRing(R, Degrees => A);
assert(componentOfKernel({1,1,1,1,0}, dom, F, matrix {{x_(1,1)*x_(2,2), x_(1,2)*x_(2,1)}}) == {x_(1,2)*x_(2,1)-x_(1,1)*x_(2,2)});

----------------------------
----- componentsOfKernel ----
--------------///
--------------
componentsOfKernel = method(Options => {Grading => null});
componentsOfKernel (Number, RingMap) := MutableHashTable => opts -> (d, F) -> (

--  A := if opts.Grading == null then maxGrading(F) else opts.Grading;
  A := if opts.Grading == null then (toList apply(0..(numgens source F -1),i->flatten entries (maxGrading(F))_i)) else opts.Grading;
--  print(A);
  dom := newRing(source F, Degrees => A);
  basisHash := new MutableHashTable;
  gensHash := new MutableHashTable;
  
  -- assumes homogeneous with normal Z-grading
  for i in 1..d do (

      B := sub(basis(i, source F), dom);
      lats := unique apply(flatten entries B, m -> degree m);
      
      G := flatten(values(gensHash));

      for deg in lats do (
        
        basisHash#deg = findBasisInDegree(deg, dom, G, basisHash);
        gensHash#deg = if numcols(basisHash#deg) == 1 then {} else componentOfKernel(deg, dom, F, basisHash#deg);
      );
  );
  
  gensHash
  )

TEST ///
A = matrix {{1,1,1,0,0,0,0,0,0}, {0,0,0,1,1,1,0,0,0}, {0,0,0,0,0,0,1,1,1}, {1,0,0,1,0,0,1,0,0}, {0,1,0,0,1,0,0,1,0}};
R = QQ[x_1..x_(numcols A)];
S = QQ[t_1..t_(numrows A)];
F = map(S, R, apply(numcols(A), i -> S_(flatten entries A_i)));
toList apply(0..(numgens source F -1),i->flatten entries (maxGrading(F))_i)
dom = newRing(R, Degrees => A);
B=toList apply(0..(numgens source F-1),i->flatten entries (maxGrading(F))_i)
dom=newRing(R, Degrees => B);
degree x_1
componentsOfKernel(2,F)
peek componentsOfKernel(2,F)
assert(componentsOfKernel(2, F) == {x_(1,2)*x_(2,1)-x_(1,1)*x_(2,2), x_(1,3)*x_(2,1)-x_(1,1)*x_(2,3), x_(1,3)*x_(2,2)-x_(1,2)*x_(2,3)});
///

-- Documentation below

beginDocumentation()

-- template for function documentation
--doc ///
--Key
--Headline
--Usage
--Inputs
--Outputs
--Consequences
--  Item
--Description
--  Text
--  Code
--  Pre
--  Example
--  CannedExample
--Subnodes
--Caveat
--SeeAlso
--///

doc ///
Key
  MultigradedImplicitization
Headline
  Package for levaraging multigradings to solve implicitization problems
Description
  Text
    The MultigradedImplicitization package provides methods for computing the maximal $\mathbb{Z}^k$ grading in which the 
    kernel of a polynomial map $F$ is homogeneous and exploiting it to find generators of $\ker(F)$. This package is
    particularly useful for problems from algebraic statistics which often involve highly structured maps {\tt F} which are
    often naturally homogeneous in a larger multigrading than the standard $\mathbb{Z}$-multigrading given by total degree.
  Text
    References:
    [1] INSERT OUR new PAPER\break
    [2] INSERT Joe's Multigraded Macaulay Dual Spaces Paper
    [3] INSERT our networks paper
    [4] Seth Sullivant, Algebraic Statistics, American Mathematical Soc.
///


doc ///
Key
  maxGrading
  (maxGrading, RingMap)
  [maxGrading, ReturnTargetGrading]
Headline
  computes the maximal $\mathbb{Z}^k$ grading such that $\ker(F)$ is homogeneous
Usage
  maxGrading(F)
Inputs
  F:RingMap
  ReturnTargetGrading => Boolean
    a boolean which encodes whether or not to also return the grading on {\tt target(F)} which induces the grading on $\ker(F)$
Outputs
  :Matrix
    the maximal $\mathbb{Z}^k$ grading in which $\ker(F)$ is homogeneous
--Consequences
--  asd
Description
  Text
    Computes the maximal $\mathbb{Z}^k$ grading such that the $\ker(F)$ is homogeneous. 
    The columns of the output matrix are the degrees of the corresponding variables in {\tt source(F)}.
    For example, the snippet below shows that the maximal grading of a toric ideal
    is exactly the integer matrix which encodes the monomial map parameterizing the ideal.
  Example
    A = matrix {{1,1,1,0,0,0}, {0,0,0,1,1,1}, {1,0,0,1,0,0}, {0,1,0,0,1,0}, {0,0,1,0,0,1}}
    R = QQ[x_(1,1)..x_(2,3)];
    S = QQ[t_1..t_2, s_1..s_3];
    F = map(S, R, {t_1*s_1, t_1*s_2, t_1*s_3, t_2*s_1, t_2*s_2, t_2*s_3})
    maxGrading(F)
  Text
    The option {\tt ReturnTargetGrading} returns a matrix which also gives the corresponding
    grading on the target ring of $F$ which induces the grading on $\ker(F)$. This option is {\tt false}
    by default. 
  Example
    maxGrading(F, ReturnTargetGrading => true)
///


doc ///
Key
  findBasisInDegree
  (findBasisInDegree, List, Ring, List, MutableHashTable)
  (findBasisInDegree, List, Ring, MutableHashTable)
Headline
  Finds a basis for the homogeneous component of a graded ring but removes basis elements which correspond to previously computed generators. 
Usage
  findBasisInDegree(deg, dom, G, B)
  findBasisInDegree(deg, dom, B)
Inputs
  deg:List
    the degree of the homogeneous component to compute
  dom:Ring
    a graded ring which is the source of a homogeneous ring map $F$
  G:List
    a list of previously computed generators of $\ker(F)$
  B:MutableHashTable
    a mutable hashtable which contains all bases of homogeneous components which correspond to lower total degrees than {\tt deg}
Outputs
  :Matrix
    A monomial basis for the homogeneous component of degree {\tt deg} of {\tt dom} with any monomials which cannot be involved in new generators of $\ker(F)$ removed.  
--Consequences
--  asd
Description
  Text
    Computes a monomial basis for the homogeneous component of degree {\tt deg} of the graded ring {\tt dom} which is the source of a ring map $F$. 
    Monomials which correspond to previously computed relations which are in {\tt G} are automatically removed since they will not yield new generators
    in $\ker(F)$ when applying @TO2{componentOfKernel, "componentOfKernel"}@ to this basis.  
  Example
    A = matrix {{1,1,1,0,0,0}, {0,0,0,1,1,1}, {1,0,0,1,0,0}, {0,1,0,0,1,0}, {0,0,1,0,0,1}}
    R = QQ[x_(1,1)..x_(2,3)];
    S = QQ[t_1..t_2, s_1..s_3];
    F = map(S, R, {t_1*s_1, t_1*s_2, t_1*s_3, t_2*s_1, t_2*s_2, t_2*s_3})
    dom = newRing(R, Degrees => A);
    B = new MutableHashTable from apply(gens(dom), i -> degree(i) => i);
    findBasisInDegree({1,1,1,1,0}, dom, B)
///


doc ///
Key
  componentsOfKernel
  (componentsOfKernel, Number, RingMap)
  [componentsOfKernel, Grading]
Headline
  Finds all minimal generators up to a given total degree in the kernel of a ring map 
Usage
  componentsOfKernel(d, F)
Inputs
  d:List
    the total degree at which to stop computing generators
  F:RingMap
  Grading => Matrix
    a matrix whose columns give a homogeneous multigrading on $\ker(F)$
Outputs
  :MutableHashTable
    A mutable hashtable whose keys correspond to all homogeneous components of $\ker(F)$ and values correspond to generators in $\ker(F)$ with those components
--Consequences
--  asd
Description
  Text
    Computes all minimal generators of $\ker(F)$ which are of total degree at most {\tt d}. 
  Example
    A = matrix {{1,1,1,0,0,0,0,0,0}, {0,0,0,1,1,1,0,0,0}, {0,0,0,0,0,0,1,1,1}, {1,0,0,1,0,0,1,0,0}, {0,1,0,0,1,0,0,1,0}};
    R = QQ[x_1..x_(numcols A)];
    S = QQ[t_1..t_(numrows A)];
    F = map(S, R, apply(numcols(A), i -> S_(flatten entries A_i)));
    peek componentsOfKernel(2, F)
  Text
    If a grading in which $\ker(F)$ is homogeneous is already known or a specific grading is desired then the option {\tt Grading} can be used to specify this.
    In this case the columns of the matrix {\tt Grading} are automatically used to grade the source of $F$. 
--  Code
--    todo
--  Pre
--    todo

--  todo
--SeeAlso
--  todo
///


doc ///
Key
  componentOfKernel
  (componentOfKernel, List, Ring, RingMap, Matrix)
  (componentOfKernel, List, Ring, RingMap, MutableHashTable)
  (componentOfKernel, List, Ring, RingMap)
  [componentOfKernel, PreviousGens]
Headline
  Finds all minimal generators of a given degree in the kernel of a ring map 
Usage
  componentOfKernel(deg, dom, F, M)
  componentOfKernel(deg, dom, F, B)
  componentOfKernel(deg, dom, F)
Inputs
  deg:List
    the degree of the homogeneous component to compute
  dom:Ring
    a graded ring which is the source of a homogeneous ring map $F$
  F:RingMap
    a map whose kernel is homogeneous in the grading of {\tt dom}
  M:Matrix
    a monomial basis for the homogeneous component of {\tt deg} with degree {\tt deg}
  B:MutableHashTable
    a mutable hashtable which contains all bases of homogeneous components which correspond to lower total degrees than {\tt deg}
  PreviousGens => List
    a list of generators of the kernel which have lower total degree
Outputs
  :Matrix
    A list of minimal generators for $\ker(F)$ which are in the homogeneous component of degree {\tt deg}
--Consequences
--  asd
Description
  Text
    Computes all minimal generators of $\ker(F)$ which are in the homogeneous component of degree {\tt deg}
  Example
    A = matrix {{1,1,1,0,0,0,0,0,0}, {0,0,0,1,1,1,0,0,0}, {0,0,0,0,0,0,1,1,1}, {1,0,0,1,0,0,1,0,0}, {0,1,0,0,1,0,0,1,0}};
    R = QQ[x_1..x_(numcols A)];
    S = QQ[t_1..t_(numrows A)];
    F = map(S, R, apply(numcols(A), i -> S_(flatten entries A_i)));
    dom = newRing(R, Degrees => A);
    componentOfKernel({1,1,1,1,0}, dom, F)
  Text
      The option {\tt PreviousGens} can be used to specify a set of previously computed generators. 
      In the case that a monomial basis or hash table of monomial bases is not given then @TO2{findBasisInDegree, "findBasisInDegree"}@
      will be used to compute a monomial basis and {\tt PreviousGens} will be used to trim this basis. 
--  Code
--    todo
--  Pre
--    todo

--  todo
--SeeAlso
--  todo
///


