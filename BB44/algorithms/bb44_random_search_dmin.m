l := 5;
m := 5;
dimension_threshold := 2; // if the dimension is below this, d_min is not computed

num_trials := 5000; // Number of valid draws (k > dimension_threshold)

// Percentage of polynomials drawn with general monomials x^a*y^b
// (0 = previous behaviour: pure powers of x or y only; 100 = general monomials only)
general_percent := 50;

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

// NEW: generates 4 distinct general monomials X^a*Y^b, drawn uniformly
// among the l*m elements of Z_l x Z_m (pure powers of X or Y may also occur)
function SampleTermsGeneral(l, m)
    error if l*m lt 4, "l*m must be at least 4 to draw 4 distinct monomials.";
    grid := {[a, b] : a in [0..l-1], b in [0..m-1]};
    return SetToSequence(RandomSubset(grid, 4));
end function;

// NEW: true if at least one monomial involves both X and Y
IsMixed := func< terms | exists{t : t in terms | t[1] ne 0 and t[2] ne 0} >;

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

BuildIBMCode := function(l, m, Aterms, Bterms)
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

function DistanceIBM(HX, HZ)
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

// Random sampling:
// - with probability general_percent %: 4 general monomials X^a*Y^b (SampleTermsGeneral)
// - otherwise, one of the 5 previous forms:
//   1: 2+2, 2: 1+3, 3: 3+1, 4: 4+0 (pure X), 5: 0+4 (pure Y)
function SampleRandomTerms(l, m, general_percent)
    if Random(1, 100) le general_percent then
        return SampleTermsGeneral(l, m);
    end if;
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
total_draws := 0;
while valid_count lt num_trials do
    total_draws := total_draws + 1;
    Aterms := SampleRandomTerms(l, m, general_percent);
    Bterms := SampleRandomTerms(l, m, general_percent);

    // Filtering without interruption: exactly 4 distinct monomials mod (l,m)
    if #red(Aterms, l, m) ne 4 or #red(Bterms, l, m) ne 4 then
        continue;
    end if;

    HX, HZ, A, B, n := BuildIBMCode(l, m, Aterms, Bterms);

    // Quick calculation of k before the distance search step (MinimumWeight)
    N := Ncols(HX);
    k := N - 2*Rank(HX);

    // Dimension filter: the draw is discarded if it is below the threshold
    if k le dimension_threshold then
        continue;
    end if;

    valid_count := valid_count + 1;

    P := DistanceIBM(HX, HZ);
    d := P[3];

    if not IsDefined(results_by_k, k) then
        results_by_k[k] := AssociativeArray();
    end if;

    if not IsDefined(results_by_k[k], d) then
        results_by_k[k][d] := [* *];
    end if;

    // 4th component: true if A or B contains a mixed monomial X^a*Y^b (a, b != 0)
    Append(~results_by_k[k][d], <Aterms, Bterms, N, IsMixed(Aterms) or IsMixed(Bterms)>);
end while;

printf "Grid (l, m) = (%o, %o), general_percent = %o, valid draws = %o, total draws = %o\n\n",
    l, m, general_percent, valid_count, total_draws;

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

        mixed_entries := [e : e in entries | e[4]];
        pure_entries  := [e : e in entries | not e[4]];

        printf "--- Parameters = [[%o, %o, %o]] (%o trial(s): %o with mixed monomials, %o pure) ---\n",
            N, k, d, #entries, #mixed_entries, #pure_entries;

        // One representative of each kind when available
        if #pure_entries gt 0 then
            e := pure_entries[1];
            printf "  Example (pure)  : A = %o  |  B = %o\n", FormatPoly(e[1]), FormatPoly(e[2]);
        end if;
        if #mixed_entries gt 0 then
            e := mixed_entries[1];
            printf "  Example (mixed) : A = %o  |  B = %o\n", FormatPoly(e[1]), FormatPoly(e[2]);
        end if;
        printf "\n";
    end for;
end for;
