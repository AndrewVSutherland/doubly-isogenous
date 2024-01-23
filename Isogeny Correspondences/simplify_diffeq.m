/*
solve the differential equations
chi(f_1(y)) = chi(f_2(b)) and chi(g_1(y)) = chi(g_2(b))
By Buium's result, a common solution (i.e. curve on which both hold) is a necessary condition for an isogeny correpondence
between the elliptic curves with lambda-invariants f_1 and f_2 and between g1 and g2.

Solving these differential equations returns a list of potential isogeny correspondences.
Some will have one of the parameters constant, some will occur due to
the parametrization having multiple choices of the parameter give the same curve, and some might be spurious 
(caused by algebraic manipulations introducing additional solutions).

check_case is the main function to call, and has options to use the specialization method as well
This file uses the direct (non-specialization) approach.

common parameters: 
name: store relations in a file for later use if not specializing (none causes things not to be stored)
printing : whether to display progress and timing in terminal
*/

load "helperfunctions.m";
load "simplify_diffeq_specialize.m";

/*this finds potential isogeny correspondences using the method above, and then checks
whether they solve the original differential equation.
This uses the direct method, not specialization.
*/
function find_correspondences(f1,f2,g1,g2  : name:="none", printing := true)
	K := Parent(f1);
	u := K.1;
 	//one variable polynomial ring in variable u

	if printing then
		print "Starting with:";
		print f1,f2,g1,g2;
		print "";

		start_time := Cputime();
	end if;
	f1_p := Derivative(f1);
	f1_pp := Derivative(f1_p);
	f1_ppp := Derivative(f1_pp);

	f2_p := Derivative(f2);
	f2_pp := Derivative(f2_p);
	f2_ppp := Derivative(f2_pp);

	g1_p := Derivative(g1);
	g1_pp := Derivative(g1_p);
	g1_ppp := Derivative(g1_pp);

	g2_p := Derivative(g2);
	g2_pp := Derivative(g2_p);
	g2_ppp := Derivative(g2_pp);

	/* Note: paper uses variables s and t.  Derivatives are with respect to t
	Code uses y and b, derivatives with respect to b
	*/

	//compute derivatives of f1(y) with respect to b, doing chain rule by hand and using y, yp, ypp, yppp as y, y', y'', y'''
	
	C := BaseRing(K); //coefficient field
	L<y,yp,ypp,yppp,b> := PolynomialRing(C,5);
	
	//Warning: we're using the specific ordering of variables
	index_y := 1;
	index_yp := 2;
	index_ypp := 3;
	index_yppp := 4;
	index_b :=5;
	
	f1y := Evaluate(f1,y);
	f1y_p := Evaluate(f1_p,y) * yp;
	f1y_pp := Evaluate(f1_pp,y) * yp^2 + Evaluate(f1_p,y)* ypp;
	f1y_ppp := Evaluate(f1_ppp,y) * yp^3 + 3 * Evaluate(f1_pp,y) * yp * ypp + Evaluate(f1_p,y) * yppp;

	//use b as variable for f2
	f2b := Evaluate(f2,b);
	f2b_p := Evaluate(f2_p,b);
	f2b_pp := Evaluate(f2_pp,b);
	f2b_ppp := Evaluate(f2_ppp,b);	
		
	//same for gs
	g1y := Evaluate(g1,y);
	g1y_p := Evaluate(g1_p,y) * yp;
	g1y_pp := Evaluate(g1_pp,y) * yp^2 + Evaluate(g1_p,y)* ypp;
	g1y_ppp := Evaluate(g1_ppp,y) * yp^3 + 3 * Evaluate(g1_pp,y) * yp * ypp + Evaluate(g1_p,y) * yppp;

	g2b := Evaluate(g2,b);
	g2b_p := Evaluate(g2_p,b);
	g2b_pp := Evaluate(g2_pp,b);
	g2b_ppp := Evaluate(g2_ppp,b);	

	//this is the differential operator from Buium
	//ypp represents second derivative of y, etc.
	chi := (2*yp*yppp-3*ypp^2)/(4*yp^2)+yp^2*(y^2-y+1)/(4*y^2*(y-1)^2);
	
	/*r3 and s3 are the relations  chi(f_1(y)) = chi(f_2(b)) and chi(g_1(y)) = chi(g_2(b))
	The goal is to find simultaneous solutions to r3(b,y,y',y'',y''')=s3(b,y,y',y'',y''')=0.
	we can get away with only computing the numerator of the difference
	*/
	
	r3 := Numerator(Evaluate(chi,[f1y,f1y_p,f1y_pp,f1y_ppp,b]) - Evaluate(chi,[f2b,f2b_p,f2b_pp,f2b_ppp,b]));
	
	s3 := Numerator(Evaluate(chi,[g1y,g1y_p,g1y_pp,g1y_ppp,b]) - Evaluate(chi,[g2b,g2b_p,g2b_pp,g2b_ppp,b]));


	/*first eliminate y''' from the system
	a general calculation shows that  chi(f1(y)) = chi(f2(b)) and chi(g1(y)) = chi(g2(b)) have
	the same y''' and y'' terms (and these derivatives appear linearly), so the difference is first order
	however, as we're computing with the numerators only we can't just subtract 
	we need to use leading coefficients
	
	r2 is denoted F(t,s,s') in the paper
	*/
	
	coeffs_r3 := Coefficients(r3,index_yppp); //get coefficients with respect to yppp
	coeffs_s3 := Coefficients(s3,index_yppp);

	r2 := coeffs_s3[2]*coeffs_r3[1]-coeffs_r3[2]*coeffs_s3[1];

	assert( (Degree(r2,index_ypp) eq 0) and (Degree(r2,index_yppp) eq 0) );
	assert ( r2 eq Numerator(r2));
	r2 := Numerator(r2); //make sure type is polynomial
	
	//Compute y'' for a solution to r2(b,y,y') = 0 by implicit differentiation
	y2 :=-(Derivative(r2,index_b)+Derivative(r2,index_y)*yp)/Derivative(r2,index_yp);	
	
	//Similarly for y''', except only numerator
	n3 := -(Derivative(r2,index_yp)*(Derivative(Derivative(r2,index_b),index_b)+2*Derivative(Derivative(r2,index_b),index_y)*yp+
		Derivative(Derivative(r2,index_b),index_yp)*y2+Derivative(Derivative(r2,index_y),index_y)*yp^2+
		Derivative(Derivative(r2,index_y),index_yp)*yp*y2+Derivative(r2,index_y)*y2)-
		(Derivative(r2,index_b)+Derivative(r2,index_y)*yp)*(Derivative(Derivative(r2,index_b),index_yp)+Derivative(Derivative(r2,index_y),index_yp)*yp+
		Derivative(Derivative(r2,index_yp),index_yp)*y2));

	//now substitute higher derivatives into r3.  Do it in stages, and be aware some computations were numerators only
	
	//remove y'' term from r3 by substituting expression for y''.  We substitute numerator, and clear the 
	//denominator which is Derivative(r2,index_yp)

	r3_coeffypp := Coefficients(r3,index_ypp); //coefficients with respect to ypp
	mp := # r3_coeffypp;
	t2 := &+ [(r3_coeffypp[1+j]*(-(Derivative(r2,index_b)+Derivative(r2,index_y)*yp))^(j)*Derivative(r2,index_yp)^(mp-1-j)) : j in [0..mp-1]];
	
	//remove y''' term from t2 by substituting, and remember only the numerator
	//we substitute n3, and clear the denominator
	//which is Derivative(r2,index_yp)^2 - apply quotient rule to y2.
	t2_numerator:= Numerator(t2); //treat t2 as a polynomial, so can get coefficients
	assert(t2_numerator eq t2);
	t2_coeffyppp := Coefficients(t2_numerator,index_yppp);
	mp := # t2_coeffyppp ;
	s2 := &+ [  t2_coeffyppp[1+j] * n3^j * Derivative(r2,index_yp) ^(2 * (mp-j-1)) : j in [0..mp-1]];
	//s2 represents the relation G(t,s,s') in the paper
	
	assert(s2 eq Numerator(s2));
	s2:= Numerator(s2);
		
	//don't want to deal with constant functions, so strip out yp=0
	r2_mod := remove_known(r2,[yp]);
	s2_mod := remove_known(s2,[yp]);
		
	assert(Numerator(s2_mod) eq s2_mod);
	assert(Numerator(r2_mod) eq r2_mod);
	
	//change type to polynomial
	s2_mod := Numerator(s2_mod);
	r2_mod := Numerator(r2_mod);
	
	if printing then
		print "Finished Basic Manipulations",Cputime(start_time);
	end if;

	/*Compute the resultant, which gives a condition on b,y 
	when there is common solution yp to s2_mod and r2_mod
	This is r1 in the code and R in the paper
	this computation is still quite fast
	*/

	r1 := Resultant(r2_mod,s2_mod,index_yp);
	
	if printing then
		print "Computed Resultant",Cputime(start_time);
	end if;


	//substitute formula for y' obtained from r1 into r2
	//yp := -Derivative(r1,index_b) / Derivative(r1,index_y);

	//or rather, compute the numerator for that relation which 
	//is what we really care about


	yp_num := -Derivative(r1,index_b);
	yp_den := Derivative(r1,index_y);
	
	r2_coeffyp := Coefficients(r2_mod,index_yp);
	mp := # r2_coeffyp;

	precomputed_powers_num:= [L!1];
	precomputed_powers_den:= [L!1];
	power_num := L!1;
	power_den := L!1;
	counter := 1;
	while counter lt mp do
		power_num := power_num * yp_num;
		power_den := power_den * yp_den;
		Append(~precomputed_powers_num, power_num);
		Append(~precomputed_powers_den, power_den);

		if printing then
			print "Powers: ",counter, " out of ",mp-1;
		end if;
		counter := counter+1;
	end while;

	//computing the powers (even just the second) is now slowest step

	rel1 := &+ [r2_coeffyp[j+1] * precomputed_powers_num[j+1] * precomputed_powers_den[mp-j] : j in [0..mp-1]];
	
	assert(rel1 eq Numerator(rel1));
	rel1 := Numerator(rel1);
	
	if printing then
		print "Computed second relation",Cputime(start_time);
	end if;

    common := GCD(rel1,r1);
	//factors are the potential isogeny correspondeces

	//might be good to store for later use, as the above computations can take a while
	if name ne "none" then
		save_relations_v4(L,f1,f2,g1,g2,name, rel1 ,r1, common);
	end if;

	if printing then
		print "Time to Compute Potential Isogeny Correspondences: ",Cputime(start_time);
	end if;

	return L,rel1,r1,common;
end function;


//verify that the function y implicity determined by relation 
//solves chi(f1(y)) = chi(f(b)) and chi(g1(y)) = chi(g2(b))
function verify_solution(relation,L,S, f1 ,f2,g1,g2 : printing:=false)
	y := L.1;
	yp := L.2;
	ypp := L.3;
	yppp := L.4;
	b := L.5;

    f1_p := Derivative(f1);
	f1_pp := Derivative(f1_p);
	f1_ppp := Derivative(f1_pp);

	f2_p := Derivative(f2);
	f2_pp := Derivative(f2_p);
	f2_ppp := Derivative(f2_pp);

	g1_p := Derivative(g1);
	g1_pp := Derivative(g1_p);
	g1_ppp := Derivative(g1_pp);

	g2_p := Derivative(g2);
	g2_pp := Derivative(g2_p);
	g2_ppp := Derivative(g2_pp);


	//compute derivatives of f1(y) with respect to b, doing chain rule by hand and using y, yp, ypp, yppp as y, y', y'', y'''
	
	f1y := Evaluate(f1,y);
	f1y_p := Evaluate(f1_p,y) * yp;
	f1y_pp := Evaluate(f1_pp,y) * yp^2 + Evaluate(f1_p,y)* ypp;
	f1y_ppp := Evaluate(f1_ppp,y) * yp^3 + 3 * Evaluate(f1_pp,y) * yp * ypp + Evaluate(f1_p,y) * yppp;

	//use b as variable for f2
	f2b := Evaluate(f2,b);
	f2b_p := Evaluate(f2_p,b);
	f2b_pp := Evaluate(f2_pp,b);
	f2b_ppp := Evaluate(f2_ppp,b);	
		
	//same for gs
	g1y := Evaluate(g1,y);
	g1y_p := Evaluate(g1_p,y) * yp;
	g1y_pp := Evaluate(g1_pp,y) * yp^2 + Evaluate(g1_p,y)* ypp;
	g1y_ppp := Evaluate(g1_ppp,y) * yp^3 + 3 * Evaluate(g1_pp,y) * yp * ypp + Evaluate(g1_p,y) * yppp;

	g2b := Evaluate(g2,b);
	g2b_p := Evaluate(g2_p,b);
	g2b_pp := Evaluate(g2_pp,b);
	g2b_ppp := Evaluate(g2_ppp,b);	

	//this is the differential operator from Buium
	//ypp represents second derivative of y, etc.
	chi := (2*yp*yppp-3*ypp^2)/(4*yp^2)+yp^2*(y^2-y+1)/(4*y^2*(y-1)^2);

    yp_formula := imp_diff(relation,S.2,S.1);

    ypp_formula := higher_impl_diff(yp_formula,yp_formula,S.2,S.1);
    yppp_formula := higher_impl_diff(ypp_formula,yp_formula,S.2,S.1);

	if printing then
    	print "Computed derivatives using relation";
	end if;

    //compute 
    chi_f1y := Evaluate(chi,[f1y,f1y_p,f1y_pp,f1y_ppp,b]);
	chif1_S := Evaluate(chi_f1y,[S.1,yp_formula,ypp_formula,yppp_formula,S.2]);

	chi_f2 := Evaluate(chi,[f2b,f2b_p,f2b_pp,f2b_ppp,b]);
	chif2_S := Evaluate(chi_f2,[S.1,0,0,0,S.2]);

    chi_g1y := Evaluate(chi,[g1y,g1y_p,g1y_pp,g1y_ppp,b]);
	chig1_S := Evaluate(chi_g1y,[S.1,yp_formula,ypp_formula,yppp_formula,S.2]);

	chi_g2 := Evaluate(chi,[g2b,g2b_p,g2b_pp,g2b_ppp,b]);
	chig2_S := Evaluate(chi_g2,[S.1,0,0,0,S.2]);

	if printing then
		print "Substituted into Differential Equations";
	end if;

	//check if differential equation satisfied (on the curve defined by the relation)
    difference1_S := Numerator(chif1_S) * Denominator(chif2_S) - Numerator(chif2_S) * Denominator(chif1_S);
	difference2_S := Numerator(chig1_S) * Denominator(chig2_S) - Numerator(chig2_S) * Denominator(chig1_S);
	relation_S := Evaluate(relation,[S.1,S.2]); //make sure in same set of variables

	div1 := IsDivisibleBy( difference1_S , relation_S);
	div2 := IsDivisibleBy( difference2_S , relation_S);
	if printing then
		print "First Diffeq:", div1, "Second Diffeq:",div2;
	end if;
	return div1 and div2;
end function;

/* Given potential isogeny correspondences, check if they are actually interesting.  Boring ones include:
1)  Ones we already knew about 
2)  Ones cause by parametrization giving the same curve multiple times
3)  Ones which involve only a single variable
The first (resp. second)) can be found through searching for isogenies of bounded degree (resp. degree 1) using modular
polynomials, as in "searching.m";
It is also possible that
4)  The algebraic manipulations involved in solving the differential equations introduce a spurious solution
To recognize 4, we can try two things.  Either use implicit differentiation to compute derivatives and check if the original differential equation is satisfies
or just specializaing the variables and work over a finite field, and verify if the elliptic curves are isogenous as claimed.
The first method seems fast enough in practice

Returns true if there is anything we didn't know about
(Also returns some auxiliary data as additional values)
*/
function check_correspondences( f1,f2,g1,g2 ,S, common, known_correspondences : name:="none", printing := true)
	s := S.1;
	t := S.2;
	L := Parent(common);
	anything_unknown:=false;
	unexplained := [];

	factors := Factorization(Evaluate(common,[s,0,0,0,t]));

	for factor in factors do
		relation := factor[1]; //can ignore multiplicity

		//is it already accounted for?
		if relation in known_correspondences then continue; end if;

		//is it only in one variable (i.e. constant)
		if Degree(relation,s) eq 0 or Degree(relation,t) eq 0 then continue; end if;


		if  verify_solution(relation,L,S, f1 ,f2,g1,g2) then 
			//this is an actual solution we didn't know about
			anything_unknown := true;
			Append(~unexplained,relation);
			if printing then
				print "Couldn't explain:",relation;
			end if;
		end if;

	end for;

	return anything_unknown,unexplained;
end function;

//return true if everything is accounted for isogeny correspondences
function verify_case(f1,f2,g1,g2 , S ,found_factors   : name:="none", printing := true)
	L,rel1,r1,common := find_correspondences(f1,f2,g1,g2  : name:=name , printing := printing);

	missed_something,suspect_relations := check_correspondences(f1,f2,g1,g2 ,S, common , found_factors : name:=name  , printing := printing);

	return not missed_something;
end function;

//continue using precomputed relations
function continue_verify_case(S,found_factors,name)
	 L,rat_funcs,rel1,r1,common := load_relations_v4(name);
	 f1,f2,g1,g2 := Explode(rat_funcs);

	missed_something,suspect_relations := check_correspondences(f1,f2,g1,g2 ,S, common,found_factors : name:=name );

	return not missed_something;
end function;

//high level function with options to compute exactly or after specialization
function check_case(f1,f2,g1,g2 , S, known_correspondences : printing := true, name:="none", specialize:=false)
	if specialize then
		return verify_case_specialize(f1,f2,g1,g2 , known_correspondences : printing:=printing);
	else
		return verify_case(f1,f2,g1,g2 , S ,known_correspondences : name:=name, printing:=printing);
	end if;
end function;