This directory contains a [magma](http://magma.maths.usyd.edu.au/magma/) program file that can be used to verify computational statements made in the proof of
Theorem 7.4 of the paper &ldquo;Doubly isogenous curves of genus two with a rational action of <i>D</i><sub>6</sub>&rdquo; by Jeremy Booher, Everett W. Howe, Andrew V. Sutherland, and José Felipe Voloch.

The one program file, extraordinary.magma, defines five functions

    verification1()
    verification2()
    verification3()
    verification4()
    verification5()

that, when run, each return the Boolean value *true*. In itself this does not prove much, but if you look at the code for the functions you will see that they verify various claims made in the proof. The statements being proved are spelled out in the comments to the code.
