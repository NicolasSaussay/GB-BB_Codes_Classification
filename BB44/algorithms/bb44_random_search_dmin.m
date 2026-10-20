// =======================================================================
// BB(4,4): random search of bivariate bicycle codes and their minimum distance
//
// Randomly samples pairs of 4-term polynomials A(X,Y), B(X,Y) over
// GF(2)[X,Y]/(X^l - 1, Y^m - 1), builds the CSS code of length N = l*m
// (H_X = [A | B], H_Z = [B^T | A^T]), checks H_X H_Z^T = 0, discards
// codes with k <= dimension_threshold, computes the minimum distance d of the
// others, and prints the results grouped by [[N, k, d]] with an example (A, B).
//
// Usage:
//   Edit l, m, dimension_threshold and num_trials at the top of the file, then:
//   magma bb44_random_search_dmin.m
//
// Output:
//   For each k, then each d: "[[N, k, d]] (number of trials)" with an example A and B.
// =======================================================================

l := 5;
m := 5;
dimension_threshold := 2; // if the dimension is below this, d_min is not computed

num_trials := 5000; // Number of valid draws (k > dimension_threshold)

red := func< T, l, m | {[t[1] mod l, t[2] mod m] : t in T} >;

// Generates 4 terms of the form X^a1 + X^a2 + Y^a3 + Y^a4
function SampleTerms2_2(l, m)
    xA := SetToSequence(RandomSubset({0..l-1}, 2));
    yA := SetToSequence(RandomSubset({0..m-1}, 2));
    return [[xA[1], 0], [xA[2], 0], [0, yA[1]], [0, yA[2]]];
end function;

// Generates 4 terms of the form X^a1 + Y^a2 + Y^a3 + Y^a4
function SampleTerms1_3(l, m)
    yA := SetToSequence(RandomSubset({0..m-1}, 3));
    xA := Random([0..l-1]);
    return [[xA, 0], [0, yA[1]], [0, yA[2]], [0, yA[3]]];
end function;

// Generates 4 terms of the form X^a1 + X^a2 + X^a3 + Y^a4
function SampleTerms3_1(l, m)
    xA := SetToSequence(RandomSubset({0..l-1}, 3));
    yA := Random([0..m-1]);
    return [[xA[1], 0], [xA[2], 0], [xA[3], 0], [0, yA]];
end function;

// Generates 4 terms of the form X^a1 + X^a2 + X^a3 + X^a4 : [[a1,0], [a2,0], [a3,0], [a4,0]]
function SampleTerms4_0(l, m)
    error if l lt 4, "l must be at least 4 to draw 4 distinct exponents in X.";
    xA := SetToSequence(RandomSubset({0..l-1}, 4));
    return [[xA[1], 0], [xA[2], 0], [xA[3], 0], [xA[4], 0]];
end function;

// Generates 4 terms of the form Y^a1 + Y^a2 + Y^a3 + Y^a4 : [[0,a1], [0,a2], [0,a3], [0,a4]]
function SampleTerms0_4(l, m)
    error if m lt 4, "m must be at least 4 to draw 4 distinct exponents in Y.";
    yA := SetToSequence(RandomSubset({0..m-1}, 4));
    return [[0, yA[1]], [0, yA[2]], [0, yA[3]], [0, yA[4]]];
end function;

// Generic formatting of a list of 4 monomials
function FormatPoly(terms)
    s := "";
    for i in [1..#terms] do
        t := terms[i];
        mono := "";
        if t[1] ne 0 and t[2] ne 0 then
            mono := Sprintf("X^%o*Y^%o", t[1], t[2]);
        elif t[1] ne 0 then
            mono := Sprintf("X^%o", t[1]);
        elif t[2] ne 0 then
            mono := Sprintf("Y^%o", t[2]);
        else
            mono := "1";
        end if;
        if i gt 1 then
            s cat:= " + " cat mono;
        else
            s := mono;
        end if;
    end for;
    return s;
end function;

function PolyMatrix(l, m, terms)
    n := l*m;
    M := ZeroMatrix(GF(2), n, n);
    for a in [0..l-1], b in [0..m-1] do
        r := a*m + b + 1;
        for t in terms do
            c := ((a + t[1]) mod l)*m + ((b + t[2]) mod m) + 1;
            M[r, c] := M[r, c] + 1;
        end for;
    end for;
    return M;
end function;

BuildBBCode := function(l, m, Aterms, Bterms)
    error if #Aterms ne 4 or #Bterms ne 4,
        "Aterms and Bterms must contain exactly 4 monomials.";

    n := l*m;
    A := PolyMatrix(l, m, Aterms);
    B := PolyMatrix(l, m, Bterms);
    HX := HorizontalJoin(A, B);
    HZ := HorizontalJoin(Transpose(B), Transpose(A));
    error if not IsZero(HX * Transpose(HZ)), "H_X * H_Z^T non-zero: non-CSS code";

    return HX, HZ, A, B, n;
end function;

function DistanceBB(HX, HZ)
    N := Ncols(HX);
    k := N - 2*Rank(HX);
    if k eq 0 then
        return <N, 0, 0>;
    end if;

    C := Dual(LinearCode(HX));   // Ker(H_X)
    D := LinearCode(HZ);         // Im(H_Z^T)
    d := MinimumWeight(C);
    while d le N do
        nD := #Words(D, d : Method := "Zimmermann");
        nC := #Words(C, d : Method := "Zimmermann", NumWords := nD + 1);
        if nC ne nD then
            return <N, k, d>;
        end if;
        d := d + 1;
    end while;
    error "Inconsistency: k > 0 but no logical operator found";
end function;


///////////////////////////
///////////////////////////


results_by_k := AssociativeArray();

// Random sampling among the 5 forms:
// 1: 2+2, 2: 1+3, 3: 3+1, 4: 4+0 (pure X), 5: 0+4 (pure Y)
function SampleRandomTerms(l, m)
    choice := Random([1, 2, 3, 4, 5]);
    if choice eq 1 then
        return SampleTerms2_2(l, m);
    elif choice eq 2 then
        return SampleTerms1_3(l, m);
    elif choice eq 3 then
        return SampleTerms3_1(l, m);
    elif choice eq 4 then
        return SampleTerms4_0(l, m);
    else
        return SampleTerms0_4(l, m);
    end if;
end function;

valid_count := 0;
while valid_count lt num_trials do
    Aterms := SampleRandomTerms(l, m);
    Bterms := SampleRandomTerms(l, m);

    // Filtering without interruption: exactly 4 distinct monomials mod (l,m)
    if #red(Aterms, l, m) ne 4 or #red(Bterms, l, m) ne 4 then
        continue;
    end if;

    HX, HZ, A, B, n := BuildBBCode(l, m, Aterms, Bterms);

    // Quick calculation of k before the distance search step (MinimumWeight)
    N := Ncols(HX);
    k := N - 2*Rank(HX);

    // Dimension filter: the draw is discarded if it is below the threshold
    if k le dimension_threshold then
        continue;
    end if;

    valid_count := valid_count + 1;

    P := DistanceBB(HX, HZ);
    d := P[3];

    if not IsDefined(results_by_k, k) then
        results_by_k[k] := AssociativeArray();
    end if;

    if not IsDefined(results_by_k[k], d) then
        results_by_k[k][d] := [* *];
    end if;

    Append(~results_by_k[k][d], <Aterms, Bterms, N>);
end while;

// Sorting and display: grouped by k, then ascending d
k_values := Sort([k : k in Keys(results_by_k)]);

for k in k_values do
    printf "====================================================\n";
    printf "                  DIMENSION k = %o                  \n", k;
    printf "====================================================\n\n";

    d_values := Sort([d : d in Keys(results_by_k[k])]);

    for d in d_values do
        entries := results_by_k[k][d];
        N := entries[1][3];
        e := entries[1]; // A representative

        polyA_str := FormatPoly(e[1]);
        polyB_str := FormatPoly(e[2]);

        printf "--- Parameters = [[%o, %o, %o]] (%o trial(s)) ---\n", N, k, d, #entries;
        printf "  Example : A = %o  |  B = %o\n\n", polyA_str, polyB_str;
    end for;
end for;
