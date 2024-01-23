/* Finding isogeny correspondences using modular polynomials

The find_isogenies does one particular instance.
The second is optimized for looking at all possibilities at once
*/


//find isogenies (j1 isog j2) as well as (j3 isog j4) where the degree of the isogeny is in a given range
function find_isogenies(s,t,j1,j2,j3,j4: start :=1, upper_bound:=5)
	isogenies1 := [];
	isogenies2 := [];

	for N in [start..upper_bound] do 
		start_time := Cputime();
		polyN := ClassicalModularPolynomial(N);
		
		isog1 := Evaluate(polyN,[Evaluate(j1,s),Evaluate(j2,t)]);
        Append(~isogenies1,Numerator(isog1));

		isog2 := Evaluate(polyN,[Evaluate(j3,s),Evaluate(j4,t)]);
        Append(~isogenies2,Numerator(isog2));

		print "Finished ",N;
	end for;

	correspondences := {};

    for i in [1..upper_bound] do
        for j in [1..upper_bound] do
            common := GCD(isogenies1[i],isogenies2[j]);
            if Degree(common) gt 0 then
                Include(~correspondences,<common,i,j>);
                print "Found something: ",i,j;
            end if;
        end for;
    end for;

    return correspondences;
end function;

/*Search for isogeny correspondences between a list of possible j-invariants.  
In the kind of examples we've been looking at, we only find low degree isogenies (<10)
This is lucky as MAGMA hasn't precomputed larger modular polynomials 
(this could be worked around for composite degree, which are the first ones left out)

Return all simultaneous correspondences found
simultaneous[relation] is a set of tuples indicating the indices that have isogenies
*/
function search_isogenies(s,t,list_js : upper_bound := 10, printing := true)
    L := # list_js;

    //indexed by N,i,j: a degree N isogeny between ith and jth j-invariants
    correspondences := AssociativeArray();

    for N in [1..upper_bound] do
        polyN := ClassicalModularPolynomial(N);
        for i in [1..L] do
            for j in [1..L] do
                correspondences[ <N,i,j>] := Numerator(Evaluate(polyN,[Evaluate(list_js[i],s),Evaluate(list_js[j],t)]));
                //could exploit symmetry
                if printing then
                    print "Finished case",N,i,j;
                end if;
            end for;
        end for;

    end for;

    if printing then
        print "Starting gcds";
    end if;

    //associative array: indexed by the relation, value a list of indexes and degrees where there is an isogeny
    simultaneous := AssociativeArray();
    for key in Keys(correspondences) do
        N1,i1,j1 := Explode(key);
        for other_key in Keys(correspondences) do
            N2,i2,j2 := Explode(other_key);

            //symmetric in key /other_key, so order and only compute one of these
            if (N1 gt N2) or (N1 eq N2 and i1 gt i2) or (N1 eq N2 and i1 eq i2 and j1 ge j2) then continue; end if;
            common := GCD(correspondences[key],correspondences[other_key]);

            if Degree(common) gt 0 then
                if common in Keys(simultaneous) then
                    Include(~simultaneous[common],<key,other_key> );
                else
                    simultaneous[common] := {<key,other_key>};
                end if;
            end if;
        end for;
        if printing then
            print "Finished key:",key;
        end if;
    end for;

    return simultaneous;
end function;

