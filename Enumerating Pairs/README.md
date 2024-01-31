# Enumerating doubly isogenous pairs of <i>D</i><sub>6</sub> curves

This repository provides [magma](http://magma.maths.usyd.edu.au/magma/) programs for computing all pairs of doubly isogenous genus-2 <i>D</i><sub>6</sub> curves over a given finite field; see the full explanation below. It is used in the paper &ldquo;Doubly isogenous curves of genus two with a rational action of <i>D</i><sub>6</sub>&rdquo; by Jeremy Booher, Everett W. Howe, Andrew V. Sutherland, and José Felipe Voloch.

# Files

There is only one file:

- fourfold.magma

# Theory

For full details, see &ldquo;Doubly isogenous curves of genus two with a rational action of <i>D</i><sub>6</sub>&rdquo; by Jeremy Booher, Everett W. Howe, Andrew V. Sutherland, and José Felipe Voloch.

We compute examples of doubly isogenous curves over finite fields, where we only look at curves of genus 2 with all their Weierstrass points rational and with geometric and rational automorphism groups both isomorphic to <i>D</i><sub>6</sub>, the dihedral group of order 12. Each such curve is given by a pair (<i>u, c</i>) of elements of <b>F</b><sub><i>q</i></sub>, as described in the paper. 

Given a prime power q that is not a power of 2 or 3, we look at all such curves <i>C</i> over <b>F</b><sub><i>q</i></sub> and compute an initial signature, which consists of the following.

- The Jacobian of such a <i>C</i> is isogenous to the square of an elliptic curve <i>E</i><sub>0</sub>. The first element of the signature is the trace of <i>E</i><sub>0</sub>.
    
- The Prym variety of the maximal abelian 2-extension of <i>C</i> in which a given Weierstrass point splits is isogenous to a product of 15 elliptic curves over <b>F</b><sub><i>q</i></sub>. The second element of the signature is the sorted list of the traces of these 15 curves.
    
We keep track of all pairs (<i>u</i><sub>1</sub>, <i>c</i><sub>1</sub>), (<i>u</i><sub>2</sub>, <i>c</i><sub>2</sub>) over <b>F</b><sub><i>q</i></sub> with the same signature.

For each such pair of curves, we compute extra information as well. For such a curve <i>C</i> there are natural degree-2 maps from <i>C</i> to two elliptic curves <i>E</i><sub>1</sub> and <i>E</i><sub>2</sub>, and there are natural 3-isogenies between <i>E</i><sub>1</sub> and <i>E</i><sub>2</sub>. The additional information we add to the signature is:

- We consider the degree-2 map <i>C</i>&rarr;<i>E</i><sub>1</sub> that takes a given Weierstrass point to the origin, and we produce the triple cover D&rarr;<i>C</i> obtained by fibering this double cover with the 3-isogeny <i>E</i><sub>2</sub>&rarr;<i>E</i><sub>1</sub>. The Prym of this triple cover of <i>C</i> is isogenous to a product of elliptic curves. The third element of the signature is the sorted list of the traces of these two curves, together with the traces obtained from the same construction with the roles of <i>E</i><sub>1</sub> and <i>E</i><sub>2</sub> swapped.
    
- We again consider the degree-2 map <i>C</i>&rarr;<i>E</i><sub>1</sub>, but now we fiber it with the multiplication-by-3 map <i>E</i><sub>1</sub>&rarr;<i>E</i><sub>1</sub>. We compute the Weil polynomial of the Prym. The fourth element of the signature is the sorted list of factors of this polynomial, combined with the factors of the Weil polynomial of the Prym obtained from the same construction with the roles of <i>E</i><sub>1</sub> and <i>E</i><sub>2</sub> swapped.

We also compute the Weil polynomial of the Prym of the pullback of the 
multiplication-by-4 map on the Jacobian of C. This pullback is a cover of the
pullback of the multiplication-by-2 map, and we already know the Weil polynomial
for that factor, so we only compute the Weil polynomial of the quotient.


# Functions


There are two main functions that the user should be aware of:

- matches(q)

For a prime power <i>q</i>, the function matches(q) runs through the outline above, and outputs two sequences [<i>a,b,c,d</i>] and [<i>e,f,g,h,i</i>], together with a list of the pairs of doubly isogenous curves over <b>F</b><sub><i>q</i></sub>, each curve represented by a sequence [<i>u,c</i>]. The elements of the two sequences are defined as follows:

1. <i>a</i> is the number of doubly isogenous pairs where the curves are twists of one another.
2. <i>b</i> is the number of pairs in (1) that are explained by a certain condition.
3. <i>c</i> is the number of pairs in (1) that remain indistinguishable when we throw in the third element of the signatures.
4. <i>d</i> is the number of pairs in (3) that remain indistinguishable when we throw in the third and fourth elements of the signatures.
5. <i>e</i> is the number of doubly isogenous pairs where the curves are *not* twists of one another.
6. <i>f</i> is the number of pairs in (4) that are reductions of the extraordinary curves in characteristic 0.
7. <i>g</i> is the number of pairs in (4) that remain indistinguishable when we throw in the third element of the signatures.
8. <i>h</i> is the number of pairs in (7) that are reductions of the extraordinary curves in characteristic 0.
9. <i>i</i> is the number of pairs in (7) that remain indistinguishable when we throw in the third and fourth elements of the signatures.

For (2), the &ldquo;certain condition&rdquo; is that the elliptic curve with lambda-invariant equal to 4<i>u</i>/((<i>u</i>  &ndash;  1)(<i>u</i>  +  3)) is supersingular.

For (6) and (7), the &ldquo;extraordinary curves&rdquo; in characteristic 0 are the curves (<i>u,c</i>) where either <i>u</i><sup>12</sup>  &ndash; 9<i>u</i><sup>10</sup> &ndash; 53<i>u</i><sup>8</sup> + 266<i>u</i><sup>6</sup> + 1707<i>u</i><sup>4</sup> + 2183<i>u</i><sup>2</sup> +  1 or <i>u</i><sup>12</sup> + 19647<i>u</i><sup>10</sup> + 138267<i>u</i><sup>8</sup> + 193914<i>u</i><sup>6</sup> &ndash; 347733<i>u</i><sup>4</sup> &ndash; 531441<i>u</i><sup>2</sup> +  531441 is equal to 0.

- signature([u,c])

Given a pair [<i>u,c</i>] of elements of <b>F</b><sub><i>q</i></sub> specifying a curve <i>C</i>, we return the signature described above, including the extended elements. Namely, the function returns a tuple of four elements:

1. The trace of the elliptic curves <i>E</i> and <i>E&prime;</i> such that Jac <i>C</i> = <i>E</i> &times; <i>E&prime;</i>.
2. The sorted list of the 15 elliptic curve traces appearing in the Prym of the maximal 2-cover of <i>C</i> described above.
3. The sorted list of the 4 elliptic curve traces appearing in the Pryms of the two special triple covers of <i>C</i> described above.
4. The sorted list of the factors of the Weil polynomials of the Pryms of the two 3-covers of <i>C</i> coming from the multiplication-by-3 maps on the elliptic factors of the Jacobian of <i>C</i>, as described above.


# Usage


Usage is straightforward. After starting Magma and loading the file via

    load "fourfold.magma";

you can run the program matches() by entering, for example,

    matches(23);
    
which produces the output

    [ 1, 0, 1, 0 ]
    [ 0, 0, 0, 0, 0 ]
    [
       [
         [ 2, 1 ],
         [ 2, 5 ]
       ]
    ]

The first two lines are the sequences described above, for the field <b>F</b><sub>23</sub>. The rest of the output is a sequence of one pair [(<i>u</i><sub>1</sub>, <i>c</i><sub>1</sub>), (<i>u</i><sub>2</sub>, <i>c</i><sub>2</sub>)], representing a pair of curves as described above.

If you want to know the signature of one of these curves, you can enter, for example,

    signature([GF(23)|2,1]);
    
which produces more output than we want to reproduce here. However, by entering

    signature([GF(23)|2,1]) eq signature([GF(23)|2,5]);
    
we get the answer

    false

which shows that while the two curves are doubly isogenous, they can be distinguished by the higher-order covers described above and in the paper.
