# Computing Isogeny Correspondences

This repository provides [magma](http://magma.maths.usyd.edu.au/magma/) programs for computing simultaenous isogeny correspondences between families of elliptic curves in characteristic zero.  It is used in the paper "Doubly isogenous curves of genus two with a rational action of D6" by Jeremy Booher, Everett W. Howe, Andrew V. Sutherland, and José Felipe Voloch.

# Files

- isogeny_corresondences.m provides setup and the high-level function process_case

- d6_family.m and d4_family.4 provide the data for two sets of families of elliptic curves. 

- searching.m uses modular polynomials to search for low degree isogenies.

- simplify_diffeq.m solves the differential equations directly, finding relations which cause isogeny correspondences.

- simplify_diffeq_specialize.m solves the differential equations after specializing a variable, and can check there are no additional isogeny correspondences given a list of known ones.

- helper_functions.m provides some basic functions for the above.

# Theory

See "Doubly isogenous curves of genus two with a rational action of D6" by Jeremy Booher, Everett W. Howe, Andrew V. Sutherland, and José Felipe Voloch.

# Usage

The process_case function is the easiest to use.  The arguments are lambda_list, a list of rational functions occurring as lambda-invariants for elliptic curves, and known_correspondences, a list of the isogeny correspondences already known about.

	load "isogeny_correspondences.m";
	//lambda_list and known_correspondences for the D6 family considered is automatically loaded
	process_case(lambda_list,known_correspondences);

It will return a sequence of tuples <a,b,c,d> representing potential isogeny correspondences where the curves with lambda-invariants lambda_list[a] and lambda_list[c] are isogenous as are lambda_list[b] and lambda_list[d].  An empty sequence (what this example gives) indicates there are no additional ones, and takes on the order of 12-24 hours to check all of the possible <a,b,c,d>.  It uses the specialization method.  By default there is quite a bit of terminal output to see how the computation is progressing.

The check_case function provides more control.  It allows specifying four rational functions to be the lambda-invariants of the four elliptic curves to be considered, as well as control over the method to be used.  This will compute all isogeny correspondences directly (without using specialization).  It takes about 5 minutes, and computations with two variable polynomials over the rationals are the bottleneck.

	check_case(lambda_list[1],lambda_list[1],lambda_list[2],lambda_list[2] , S, []); 	 

Since this didn't specify the known isogeny correspondences, the output indicates that there a bunch of unexplained solutions to the differential equations like s-t=0.

Specializaing is faster, but requires knowning the isogeny correspondences in which case the computation verifies there are no additional ones.
	
	//found_factors precomputed for this example in d6_family.m
	check_case(lambda_list[1],lambda_list[1],lambda_list[2],lambda_list[2] , S, found_factors: specialize:=true); 

The computation takes 5 seconds and returns true, indicating there are no additional isogeny correspondences.
