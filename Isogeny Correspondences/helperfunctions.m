/*
Basic helper functions, for removing known factors of polynomials, computing implicit derivatives
reducing rational functions at a place, and saving relations for later use
*/


/*given a relation, remove the known factors in list from it
will remove powers of the given factor
make sure the relation and known factors are typed as a polynomial!
*/
function remove_known(rel,list)
	if rel eq 0 then
		return 0;
	end if;
	
	assert (rel eq Numerator(rel));
	rel := Numerator(rel); //make sure it's typed as polynomial, not rational function
	
	new_rel := rel;

	for term in list do
		poly := Numerator(term);
				
		new_remainder := new_rel;
		divisible := true;
	
		while divisible do
			new_rel := new_remainder;
			assert (Numerator(new_rel) eq new_rel);  //if things were wrong type, could divide in field of fractions forever
			divisible,new_remainder := IsDivisibleBy(new_rel,poly); //if returns false, new_remainder isn't defined
		end while; 
	end for;
	
	return new_rel;
end function;

/*
given known factors (including multiplicities) pull them out of a relation
*/
function pullout(rel,list : printing:=true)
	if rel eq 0 then
		return 0;
	end if;
	
	assert (rel eq Numerator(rel));
	new_rel := Numerator(rel); //make sure it's typed as polynomial, not rational function
	counter := 1;


	for term in list do
		poly := Numerator(term);
		if printing then
			print "Removing: ",counter;
		end if;
		divisible,new_remainder := IsDivisibleBy(new_rel,poly);
		assert(divisible); //if this fails, we were wrong about factor		
		new_rel := new_remainder;
		counter := counter +1;
	end for;
	
	return new_rel;
end function;


//implicit differentiation relation R, variables x (dependent) and y (independent)
function imp_diff(R,x,y)
    return -1 * Derivative(R,x) / Derivative(R,y);
end function;

//chain rule basically
function higher_impl_diff(formula,yp,x,y)
    //Derivative only works on polynomials
    num := Numerator(formula);
    den := Denominator(formula);
    return (Derivative(num,x) + Derivative(num,y) * yp) / den - num / den^2 * (Derivative(den,x) + Derivative(den,y) * yp );
end function;

//reduce a rational function with coefficients in a number field C  at a place P of C. 
//t is variable in poly ring over residue field (or extension thereof)
function reduce_rat(f,t,P)
	num := Coefficients(Numerator(f));
	den := Coefficients(Denominator(f));
	num_red := &+[Evaluate(num[i],P) * t^(i-1) : i in [1..#num]];
	den_red := &+[Evaluate(den[i],P) * t^(i-1) : i in [1..#den]];
	return num_red / den_red;
end function;



ARCHIVE_PATH := "precomputed_";

//WARNING: Magma sometimes has trouble saving objects and causes strange errors
//this might be caused by relations being too big
//don't save things (keep the default name:="none") if this is an issue
procedure save_relations_v4(L,f1,f2,g1,g2,name,r1, rel1 , common)
	I := Open(ARCHIVE_PATH cat name cat ".dat","w");

	print Parent(r1);

	WriteObject(I,"v4");
	WriteObject(I,L);
	WriteObject(I,[f1,f2,g1,g2]);
	WriteObject(I,r1);
	WriteObject(I,rel1);
	WriteObject(I,common);
	
	delete I; //closes file
	return;
end procedure;

function load_relations_v4(name)
	I := Open(ARCHIVE_PATH cat name cat ".dat","r");
	version :=ReadObject(I);
	assert (version eq "v4");
	L:=ReadObject(I);
	rat_funcs:= ReadObject(I);
	//note that serialization loses track of types, so explicitly making sure everything in L
	rel1:= L!ReadObject(I);
	r1:= L!ReadObject(I);
	common := L!ReadObject(I);

	delete I; //closes file
	return L,rat_funcs,rel1,r1,common;
end function;

/*a nicer format for found isogeny correspondences
known relations are tuples < relation, set>
where set contains tuples of the form <N,i,j> to indicate that there is a degree N isogeny from
the ith curve  to the jth curve when the relation is satisfied.
*/
function process_isogeny_correspondences(known_correspondences : keep_degree:=true)
    known_relations := [];
    for pair in known_correspondences do
        relation := pair[1];
        new_S := {};
        for tuples in pair[2] do
			if keep_degree then
            	new_S := new_S join {<tuples[1][1],tuples[1][2],tuples[1][3]>, <tuples[2][1],tuples[2][2],tuples[2][3]>};
			else
            	new_S := new_S join {<tuples[1][2],tuples[1][3]>, <tuples[2][2],tuples[2][3]>};
			end if;
        end for;
        Append(~known_relations,<relation,new_S>);

    end for;
	return known_relations;
end function;

/*
known_relations are relations which cause simultaneous isogeny correspondences, the output of process_isogeny_correspondences
boring_relations are relations caused by duplication in the parametrization

Return true if there are any interesting relations (not explained by the parametrization)
*/
function check_relations(known_relations,boring_relations)
	exceptions := [];

	for item in known_relations do
		found_relations:= { factor[1] : factor in Factorization(item[1])};
		interesting := found_relations diff boring_relations;
		if # interesting gt 0 then
			Append(~exceptions,item);
		end if;
	end for;

	return exceptions;
end function;

