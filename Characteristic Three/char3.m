k<a>:=RationalFunctionField(GF(9));
F<x>:=FunctionField(k);
R<y>:=PolynomialRing(F);
K<s>:=ext<F|y^3+y^2+y-(x^6+x^4+x^2+a)*(x^3+x)^2>;
Genus(K);
P:=Zeros(s)[1];
P;
S:=Basis(3*P);
f:=S[1];
g:=S[2];
f;
g;
g^2-(f^3+f^2+f-a^2);
