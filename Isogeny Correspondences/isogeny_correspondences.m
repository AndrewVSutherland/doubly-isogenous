/*
Computing Simultaneous Isogeny Correspondences

Let K be a field of characteristic zero.  Given families of elliptic curves over K(u)
E_i : y^2 = x(x-1)(x-lambda(u))
for i = 1,2,3,4, where lambda is a rational function of a parameter u, this will
find all simultaneous isogeny correspondences.  These consist of
irreducible polynomials f(y,b) such that over f: Z = V(f) \subset A^2
we have isogenies
(\pi_1 \circ f)^* (E_1) ~ (\pi_2 \circ f)^*(E_3) 
and 
(\pi_1 \circ f)^* (E_2) ~ (\pi_2 \circ f)^*(E_4) 
where \pi_i : A^1_y \times A^1_b \to A^1_u are projection maps sending y or b to u.

1)  We can produce a list of isogenies using low degree modular polynomials (searching.m).
This is not provably complete as we cannot rule out the existence of a high degree isogeny

2)  We can produce a complete list of isogenies by solving a system of differential equations.
The direct method (simplify_diffeq.m) is slower but will compute a complete set of correspondences
An indirect method (simplify_diffeq_specialize.m) is much faster as many computations are done after 
specializing one of the variables.  It can rule out the existence of additional isogeny correspondences
if you already know a complete set (say from method 1).

For a complete discussion, see "Doubly isogenous curves of genus two with a rational action of D_6"
by Jeremy Booher, Everett W. Howe, Andrew V. Sutherland, and José Felipe Voloch
*/



load "simplify_diffeq.m";
load "searching.m";

//family in this paper
load "d6_family.m";

//family in Math. Comp. paper
//load "d4_family.m";

/*Variables provided by family files:
relevant polynomial rings
lambda_list
j_list
known_correspondences
*/

/*If not precomputed, can create list of known correspondences using
search_isogenies(s,t,j_list);
*/


/*
Search for simultaenous correspondences between elliptic curves with lambda-invariants lambda_list
known_correspondences provide a list of known correspondences. 
 This uses the more efficient specialization method, so can only rule out additional correspondences

 returns a list of problematic indices in which we either missed a correspondence
*/
function process_case(lambda_list,known_correspondences)
    S := Parent(known_correspondences[1][1]);

    //convert format so known_relations is tuples of relations and a set of tuples representing isogenies
    known_relations := [];
    for pair in known_correspondences do
        relation := pair[1];
        new_S := {};
        for tuples in pair[2] do
            new_S := new_S join {<tuples[1][2],tuples[1][3]>, <tuples[2][2],tuples[2][3]>};
        end for;
        Append(~known_relations,<relation,new_S>);
    end for;

    problematic_ones:=[];
    N := # lambda_list;

    for i1 in [1..N] do
        for i2 in [1..N] do
            for i3 in [1..N] do
                for i4 in [1..N] do
                    //look at (i1,i2) ~ (i3,i4)

                    //symmetry allows us to skip
                    if i1 gt i3 or (i1 eq i3 and i2 gt i4) then continue; end if;
                    //asking about the same isogeny twice
                    if (i1 eq i3) and (i2 eq i4) then continue; end if;

                    known_factors := {};
                    for pair in known_relations do
                        if <i1,i2> in pair[2] and <i3,i4> in pair[2] then
                            known_factors := known_factors join {fact[1] : fact in Factorization(pair[1])};
                        end if;
                    end for;

                    result := check_case(lambda_list[i1],lambda_list[i2], lambda_list[i3],lambda_list[i4] , S, known_factors: specialize:=true);

                    if not result then 
                        print "WARNING: Problem with Case", i1,i2,i3,i4; 
                        Append(~problematic_ones,<i1,i2,i3,i4>);
                    end if;

                end for;
            end for;
        end for;
    end for;

    if #problematic_ones gt 0 then
        print "There were problems:";
        print problematic_ones;
    end if;

    return problematic_ones;
end function;

