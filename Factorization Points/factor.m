/*
    This code is related to Section 10 of the paper.

    We use the five elliptic curves  whose twists are factors of the Prym associated
    to the double cover of D_{u,c}.  We fix c=1 and e=1 (so trivial quadratic twists)
*/
ZZ := Integers();
F<u> := FunctionField(Rationals());
r := -u^2*(u^2-9)^2/(u^2-1)^2;
ainvs := [[0,r-18,0,81-2*r,r],
          [0,-1-l,0,l,0] where l:=4*u/((u-1)*(u+3)),
          [0,-1-l,0,l,0] where l:=16*u^2/(u^2+3)^2,
          [0,-1-l,0,l,0] where l:=(u-3)^2*(u+1)^2/(u^2+3)^2,
          [0,-1-l,0,l,0] where l:=(u-1)^2*(u+3)^2/(u^2+3)^2];

// The code below puts the curve equations and j-invariants in a nicer form, removing redundant factors
minvs := [Coefficients(MinimalModel(EllipticCurve(a))):a in ainvs];
assert &and[&and[Denominator(a) eq 1:a in b]:b in minvs];
R<u>:=PolynomialRing(Rationals());
minvs := [[R|a:a in b]:b in minvs];
fden := func<f|LCM([ZZ|Denominator(c):c in Coefficients(f)])>;
fnum := func<f|LCM([ZZ|Numerator(c):c in Coefficients(f)])>;
w := [1,2,3,4,6];
for i:=1 to #minvs do
    a := minvs[i];
    den := LCM([fden(f):f in a]);
    P := PrimeDivisors(den);
    for p in P do
        e := Max([Ceiling(Valuation(fden(a[i]),p)/w[i]):i in [1..5]]);
        a := [p^(e*w[i])*a[i]:i in [1..5]];
    end for;
    assert LCM([fden(f):f in a]) eq 1;
    num := LCM([fnum(f):f in a]);
    P := PrimeDivisors(den);
    for p in P do
        e := Min([Floor(Valuation(fnum(a[i]),p)/w[i]):i in [1..5]]);
        a := [p^(-e*w[i])*a[i]:i in [1..5]];
    end for;
    assert LCM([fden(f):f in a]) eq 1;
    assert IsIsomorphic(EllipticCurve([F!b:b in a]),EllipticCurve([F!b:b in minvs[i]]));
    minvs[i] := a;
end for;
R<u>:=PolynomialRing(Integers());
minvs := [[R|a:a in b]:b in minvs];
jinvs := [jInvariant(EllipticCurve([F|b:b in a])):a in minvs];
jnums := [R|Numerator(j):j in jinvs];
jdens := [R|Denominator(j):j in jinvs];
jinvsmodp := func<x|[ZZ|d ne 0 select ZZ!(Evaluate(jnums[i],x)/d) else #Parent(x) where d:=Evaluate(jdens[i],x):i in [1..#jinvs]]>;

// The code below sanity checks the formulas used in the low level C code to compute jinvs
u2:=u^2; t1:=u^2-1; t2:=u^2-9; v:=u^2+3;
d2:=u*t1*t2; d1:=d2*t1^2; w:= d2*v^2; d3:=w*u;
d4:=w*(u-3)*(u+1); d5:=w*(u+3)*(u-1);
d1:=d1^2; d2:=d2^2; d3:=d3^2; d4:=d4^2; d5:=d5^2;
assert [d1,d2,d3,d4,d5] eq jdens;
f1:=u^3-15*u^2+75*u+3; f2:=u^4-4*u^3+214*u^2-36*u+81;
f3:=u^8-4*u^7+20*u^6+52*u^5-26*u^4-156*u^3+180*u^2+108*u+81;
n1:=v*Evaluate(f1,u^2); n2:=2*v^2; n3:=Evaluate(f2,u^2);
n4:=2*f3; n5:=2*Evaluate(f3,-u);
n1:=n1^3; n2:=2*n2^3; n3:=n3^3; n4:=2*n4^3; n5:=2*n5^3;
assert [n1,n2,n3,n4,n5] eq jnums;

// For faster testing we create a lookup table of absolute values of traces indexed by j-invariant
// For bad u we just set the trace to -1. Note that it doesn't matter which twist of j=0,1728 we use.
jtraces := func<J,x,shifts|&cat[[J[Integers()!j+1]:j in js] where js:=jinvsmodp(x+s): s in shifts]>;
jmap := func<p|[Abs(TraceOfFrobenius(EllipticCurveWithjInvariant(j))):j in GF(p)] cat [-1]>;
jgoodp := func<p,s|#{jtraces(J,x,s):x in GF(p)} eq p where J:=jmap(p)>;

// If our fast test fails we do a more careful test using exact traces and excluding bad u
subs := func<f,x|&+[Parent(x)|c[i]*x^(i-1):i in [1..#c]] where c:=Coefficients(f)>;
trace := func<a,x|Discriminant(R![F|subs(a[5],x),subs(a[4],x),subs(a[2],x),1]) eq 0 select #F
                  else TraceOfFrobenius(EllipticCurve([subs(c,x):c in a])) where R:=PolynomialRing(F) where F:=Parent(x)>;
traces := func<x,shifts|&cat[[ZZ|trace(a,x+s):a in minvs] : s in shifts]>;
goodp := func<p,s|#b eq #Set(b) where b := [r:r in a|Max(r) lt p] where a:=[traces(x,s):x in GF(p)]>;

shifts := {0,-22,-30};  // these shifts appear to work for all p > 3, verified up to 10^6

// this takes about 5 minutes
printf "Checking shifts %o for primes p in [5,10000]", shifts;
time for p in PrimesInInterval(5,10000) do if jgoodp(p,shifts) or goodp(p,shifts) then printf "."; else print p; end if; end for;
print "done";
