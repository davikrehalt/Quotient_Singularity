--Warning!!! This code is WRONG! because the Nullspace function in the Macaulay2 library is bugged!
dotProd = (n, a, b) -> (
    s = 0;
    for i from 0 to n-1 do (
        s = s + a#i * b#i
    );
    s
)

sumTuples = (a, b) -> (
    apply(0..#a-1, i -> a#i + b#i)
)
grid = (dimension, groupsize, weightlist) -> (
    outputlist = {};
    for i from 1 to (groupsize+1)^dimension - 1 do (
        coordinates = apply(0..dimension-1, d -> (i//(groupsize+1)^d)%(groupsize+1));
        i;
        if (dotProd(dimension,coordinates,weightlist)% groupsize == 0) then (
            outputlist = append(outputlist, coordinates);
        )
    );
    outputlist
)

trimGenerators = (preGenerators) -> (
    previousStep = sort preGenerators;
    currentStep = set preGenerators;
    counter = 0;
    while (counter < #previousStep) do (
        for j from counter to #previousStep - 1 do (
            toDiscard = sumTuples(previousStep#counter, previousStep#j);
            currentStep = currentStep - set {toDiscard};
        );
        
        previousStep = sort toList currentStep;
        counter = counter + 1;
    );

    return previousStep;
);

blockDiag = (m1, m2) -> (
    rows1 = #m1;
    if rows1 == 0 then return m2;

    cols1 = #(m1#0);
    if cols1 == 0 then return m2;

    rows2 = #m2;
    cols2 = #(m2#0);

    -- Create the new list with the required dimensions
    newMatrix = {};
    for i from 0 to (rows1 + rows2 - 1) do (
        newRow = {};
        for j from 0 to (cols1 + cols2 - 1) do (
            if i < rows1 then (
                if j < cols1 then (
                    newRow = append(newRow, m1#i#j);
                ) else (
                    newRow = append(newRow, 0);
                );
            ) else (
                if j < cols1 then (
                    newRow = append(newRow, 0);
                ) else (
                    newRow = append(newRow, m2#(i-rows1)#(j-cols1));
                );
            );
        );
        newMatrix = append(newMatrix, newRow);
    );
    return newMatrix;
);
transposeList = (m) -> (
    rows = #m;
    if rows == 0 then return {};

    cols = #(m#0);
    if cols == 0 then return {};

    newMatrix = {};
    for j from 0 to (cols - 1) do (
        newRow = {};
        for i from 0 to (rows - 1) do (
            newRow = append(newRow, m#i#j);
        );
        newMatrix = append(newMatrix, newRow);
    );
    return newMatrix;
);

sublist = (lst, indicess) -> (
    newList = {};
    for i in indicess do (
        newList = append(newList, lst#i);
    );
    return newList;
);

hstack = (m1, m2) -> (
    if #m1 == 0 then return m2;
    if #m2 == 0 then return m1;
    return toList apply(0 .. (#m1 - 1), i -> m1#i | m2#i)
);

generatePairs = (n) -> (
    pairslist = {};
    for i from 0 to n-2 do (
        for j from i+1 to n-1 do (
            pairslist = append(pairslist , {i, j});
        )
    );
    return pairslist;
);

matrixRank = (m) -> (
    nrows = #m; 
    if nrows == 0 then return 0; 
    ncols = #m#0; 
    if ncols == 0 then return 0; 
    rring = RR; 
     
    mutableM = mutableMatrix(rring, nrows, ncols); 
         
    for i from 0 to nrows - 1 do ( 
        for j from 0 to ncols - 1 do ( 
            mutableM_(i,j) = m#i#j; 
        ) 
    ); 
     
    return rank(matrix mutableM); 
     
);  


nullsspace = (m) -> (
    nrows = #m;
    if nrows == 0 then return toList {};
    ncols = #m#0;
    if ncols == 0 then return toList {};
    rring = RR;
    
    mutableM = mutableMatrix(rring, nrows, ncols);
    
    for i from 0 to nrows - 1 do (
        for j from 0 to ncols - 1 do (
            mutableM_(i,j) = m#i#j*1.0;
        )
    );

    nullss = nullSpace(mutableM);
    nullsss= toList apply(0..numRows(nullss)-1, i -> toList apply(0..numColumns(nullss)-1,j -> nullss_(i,j)));

    return nullsss;
);

getT1 = (R, ggenerators) -> (
  N = #ggenerators;
  dimension = #ggenerators#0;
  
  ECoordinateDiag = {};
  ECoordinates = {};
  E1 = new MutableHashTable;
  indexLookup = new MutableHashTable;

  for d from 0 to dimension - 1 do (
    checkCondition = (n, d) -> (ggenerators#n#d) < R#d;
    E1#d = select(toList(0 .. N-1), n -> checkCondition(n, d));

    if #E1#d > 0 then (
      counter = 0;
      for point in E1#d do (
        indexLookup#(d, point) = if (#ECoordinateDiag == 0) then (counter) else (counter + #ECoordinateDiag#0);
        counter = counter + 1;
      );

      currentInclusionMatrix = transposeList(apply(E1#d, i -> ggenerators#i));
      ECoordinateDiag = blockDiag(ECoordinateDiag, currentInclusionMatrix);
      ECoordinates = hstack(ECoordinates, currentInclusionMatrix);
    );
  );

  -- Compute E2 and SecondDifferential
  E2 = new MutableHashTable;
  SecondDifferential = {};
  pairsList = generatePairs(dimension);

  for pair in pairsList do (
    d1 = pair#0;
    d2 = pair#1;
    E2#pair = sort toList(set(E1#d1)*set(E1#d2));

    if ((#E2#pair > 0) and (#ECoordinateDiag > 0)) then (
      newCols = toList apply(E2#pair, i -> (
        toList apply(#ECoordinateDiag#0, j ->
          if j == indexLookup#(d1, i) then (
            1.0
          )
          else if j == indexLookup#(d2, i) then (
            -1.0
          )
          else 0.0
        )
      ));
      SecondDifferential = hstack(SecondDifferential, transposeList(newCols));
    );
  );

  -- Compute KERNELDIM
  KERNELDIM = (if #ECoordinates == 0 then 0 else #ECoordinates#0 - matrixRank(ECoordinates));

  -- Compute IMAGEDIM
  NullSpaceOfECoordinates = nullsspace(ECoordinateDiag);
  CombinedMatrix = hstack(NullSpaceOfECoordinates, SecondDifferential); 
  if #CombinedMatrix > 0 then (
    IMAGEDIM = matrixRank(CombinedMatrix);
  ) else (
    IMAGEDIM = 0;
  );
  
  return KERNELDIM - IMAGEDIM;
);

hashTableToLists = h -> (
    keyList = toList keys h;
    toList apply(keyList, k -> {k,h#k})
);
compute2d = (groupsize, weightlist) -> (
    if groupsize == 1 then return {};
    Grid = grid(2, groupsize, weightlist);
    ggenerators = trimGenerators(Grid);
    output = new MutableHashTable;
    scan(0..#Grid-1, i -> (
        T1 = getT1(Grid#i, ggenerators);
        if T1 != 0 then output#(Grid#i)=T1;
    ));
    return hashTableToLists(output);
);
checkInterior = (dimension, groupSize, weightList) -> (
    if dimension != 3 then (
        return "NotImplementedError";
    );

    -- Pre-check
    if gcd(weightList#0, weightList#1, groupSize) != 1 then (
        error "ValueError";
    );
    if gcd(weightList#0, weightList#2, groupSize) != 1 then (
        error "ValueError";
    );
    if gcd(weightList#1, weightList#2, groupSize) != 1 then (
        error "ValueError";
    );

    -- Check if more than two prime factors
    if groupSize < 6 then (
        return new MutableHashTable;
    );
    factors = factor groupSize;
    if #factors < 2 then (
        return new MutableHashTable;
    );

    gcdList = apply(0 .. 2, i -> gcd(weightList#i, groupSize));
    sumValue = 0;
    for i from 0 to 2 do (
        sumValue = sumValue + if gcdList#i > 1 then 1 else 0;
    );
    if sumValue < 2 then (
        return new MutableHashTable;
    );

    -- Generate grid
    Grid = {};
    for i from 0 to groupSize // gcdList#0 do (
        for j from 0 to groupSize // gcdList#1 do (
            for k from 0 to groupSize // gcdList#2 do (
                if (weightList#0*i + weightList#1*j + weightList#2*k) % groupSize == 0 and (i,j,k) != (0,0,0) then (
                    Grid = append(Grid, (i,j,k));
                );
            );
        );
    );

    ggenerators = trimGenerators(Grid);  -- Assuming trimGenerators works similarly to the trim_generators in your Python code
    print(Grid);
    print(ggenerators);

    -- Create R_dict
    RDict = new MutableHashTable;

    for R in Grid do (
        if R#0 > 0 and R#1 > 0 and R#2 > 0 then (
            toCheck = {};
            if R#0 > 1 and gcdList#1 >= R#0 and gcdList#2 >= R#0 then (
                toCheck = append(toCheck, 0);
            );
            if R#1 > 1 and gcdList#0 >= R#1 and gcdList#2 >= R#1 then (
                toCheck = append(toCheck, 1);
            );
            if R#2 > 1 and gcdList#0 >= R#2 and gcdList#1 >= R#2 then (
                toCheck = append(toCheck, 2);
            );
            if #toCheck > 0 then (
                print(R);
                T1ofR = getT1(R, ggenerators);  -- Assuming getT1 works similarly to getT1 in your Python code
                if T1ofR > 0 then (
                    RDict#R = T1ofR;
                );
            );
        );
    );

    return RDict;
);
compute3d = (groupSize, weightList) -> (
    -- Assertions
    if groupSize <= 1 then (
        error "AssertionError: groupSize must be greater than 1";
    );
    if #weightList != 3 then (
        error "AssertionError: weightList must contain exactly 3 elements";
    );
    for w in weightList do (
        if w < 0 then (
            error "AssertionError: All weights must be non-negative";
        );
    );

    -- Fast threefold calculation
    Interior = checkInterior(3, groupSize, weightList);  -- Assuming checkInterior works like check_interior in Python

    -- Compute what the boundaries are
    gcdList = apply(0 .. 2, i -> gcd(weightList#i, groupSize));

    Boundaries = {};

    weights = (weightList#1 % gcdList#0, weightList#2 % gcdList#0);
    Boundaries = join(Boundaries, toList apply(compute2d(gcdList#0, weights), entry -> {(-1, entry#0#0, entry#0#1), entry#1}));

    weights = (weightList#0 % gcdList#1, weightList#2 % gcdList#1);
    Boundaries = join(Boundaries, toList apply(compute2d(gcdList#1, weights), entry -> {(entry#0#0, -1, entry#0#1), entry#1}));

    weights = (weightList#0 % gcdList#2, weightList#1 % gcdList#2);
    Boundaries = join(Boundaries, toList apply(compute2d(gcdList#2, weights), entry -> {(entry#0#0, entry#0#1, -1), entry#1}));

    -- Combine Interior and Boundaries
    for l in keys Interior do (
        print(l);
        Boundaries = append(Boundaries, {l,Interior#l});
    );

    return toList Boundaries;
);
print compute3d(12,(1,4,9));
