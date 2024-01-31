This directory conatins a [magma](http://magma.maths.usyd.edu.au/magma/) program that initializes variables containing data relevant to the paper &ldquo;Doubly isogenous curves of genus two with a rational action of <i>D</i><sub>6</sub>&rdquo; by Jeremy Booher, Everett W. Howe, Andrew V. Sutherland, and José Felipe Voloch.

The one program file in this directory, isogenydata.magma, does nothing except assign values to 19 variables named

    data12
    data13
    data14
    .
    .
    .
    data30

Each variable data<sub><i>n</i></sub> contains data for the 1024 primes closest to 2<sup><i>n</i></sup>. (Note that for the values of <i>n</i> in question we do not have to worry about breaking ties in choosing these 1024 primes.) Each variable is a sequence of pairs

    [ [ GF(p) | u1, c1 ], [ GF(p)| u2, c2 ] ]

that specify curves over <b>F</b><sub><i>p</i></sub> denoted in the paper as <i>D</i><sub><i>u</i><sub>1</sub>,<i>c</i><sub>1</sub></sub> and <i>D</i><sub><i>u</i><sub>2</sub>,<i>c</i><sub>2</sub></sub> (see Notation 3.2) that are doubly isogenous to one another. Every pair of doubly isogenous curves over a field <b>F</b><sub><i>p</i></sub> in the appropriate range is included in the sequence.

Using this data and the function 

    signature()

in the file fourfold.magma found elsewhere in this repository, it is easy to recreate the information presented in Tables 2&ndash;5 of the paper.    
