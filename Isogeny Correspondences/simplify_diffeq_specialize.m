/*
See the documentation for simplify_diffeq.m for what this is doing.
We specialize one of the variables to speed up computations, at the expense of only being able to check
if there are no further isogeny correspondences beyond a list we know about and removed

See the paper for more details
*/

/*
Verify a solution after specialization, using point counting over a finite field to check for isogeny
which is 1 if variable remaining is y, 5 if b
other_val is what we specialized the other variable at
*/
function verify_solution_specialize(f1,f2,g1,g2,relation_L,which,other_val)
	C:= BaseRing(Parent(f1));
    F<t> := PolynomialRing(C);
    if which eq 1 then
        relation := Evaluate(relation_L,[t,0,0,0,other_val]);//other variable shouldn't show up so other_val is superflyous
    else
        relation := Evaluate(relation_L,[other_val,0,0,0,t]);
    end if;

	if Degree(relation) eq 0 then
		return false;
	end if;

    // find a prime to reduce mod p, where we have a point over a small extension of F_p
    p:=10^4;
    deg :=0;
    min_deg:=10^4;
    min_p := 0;
    tries := 0;

    while deg ne 1 and tries lt 1000 do
        p:=NextPrime(p);
		P := Decomposition(C,p)[1][1];
		k:=ResidueClassField(P);
        R<z> := PolynomialRing(k);
        relation_modp := Numerator(reduce_rat(relation,z,P)); //doesn't have denominator

        factorization:= Factorization(relation_modp);	
        
        deg := Maximum( [Degree(fact[1]): fact in factorization]);
        
        if deg lt min_deg then
            min_deg := deg;
            min_p := p;
        end if;

        tries := tries + 1;
    end while;


    if min_deg gt 8 then
        print "Warning: Didn't find a low degree factor";
    end if;

    p:=min_p;
	P := Decomposition(C,p)[1][1];
	base:=ResidueClassField(P);
    k := ext<base | min_deg>;

    R<x> := PolynomialRing(k);
    relation_modp := Numerator(reduce_rat(relation,x,P)); //doesn't have denominator
    factorization:= Factorization(relation_modp); 

    //this will have a linear factor
    for factor in factorization do
        if Degree(factor[1]) eq 1 then
            root := -Coefficients(factor[1])[1];
            break;
        end if;
    end for;

    assert (Evaluate(relation_modp,root) eq 0);

    if which eq 1 then //specialized b
        specialize_y := root;
        specialize_b := k!other_val; 
    else
        specialize_b := root;
        specialize_y := k!other_val; 
    end if;

	//type errors if not working over rationals
    //assert( Evaluate(relation_L,[specialize_y,0,0,0,specialize_b]) eq 0);

	f1_bar := reduce_rat(f1,x,P);
	f2_bar := reduce_rat(f2,x,P);
	g1_bar := reduce_rat(g1,x,P);
	g2_bar := reduce_rat(g2,x,P);

    l1 := Evaluate(Numerator(f1_bar),specialize_y) / Evaluate(Denominator(f1_bar),specialize_y);
    l2 := Evaluate(Numerator(f2_bar),specialize_b) / Evaluate(Denominator(f2_bar),specialize_b);
    l3 := Evaluate(Numerator(g1_bar),specialize_y) / Evaluate(Denominator(g1_bar),specialize_y);
    l4 := Evaluate(Numerator(g2_bar),specialize_b) / Evaluate(Denominator(g2_bar),specialize_b);

    curves := [EllipticCurve( x * (x-1)*(x- l ) ): l in [l1,l2,l3,l4]];

    return IsIsogenous(curves[1],curves[2]) and IsIsogenous(curves[3],curves[4]);
end function;



/*
solve the differential equations
chi(f_1(y)) = chi(f_2(b)) and chi(g_1(y)) = chi(g_2(b))
By Buium's result, a common solution (i.e. curve on which both hold) is a necessary condition for an isogeny correpondence

This returns true if there are no unexplained ones.
We first need to remove known ones, as otherwise can't distinguish after specialization

printing : whether to display progress in terminal
*/

function verify_case_specialize(f1,f2,g1,g2 , known_correspondences  : printing := true)
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
	This is r1 in the code and H in the paper
	this computation is still quite fast
	*/

	r1 := Resultant(r2_mod,s2_mod,index_yp);


	
	if printing then
		print "Computed Resultant",Cputime(start_time);
	end if;

    //it is actually feasible to factor r1 despite it being rather big, which
    //helps us pull out factors

    factorization := Factorization(r1);

    if printing then
        print "Factored Resultant",Cputime(start_time);
    end if;
    

    known_factors_L := [ Evaluate(factor,[y,b]) : factor in known_correspondences];

    r1_mod := L!1;

    for factor in factorization do
        if Degree(factor[1],b) eq 0 or Degree(factor[1],y) eq 0 or factor[1] in known_factors_L then
            //don't use
            if printing then
                print "Reject",factor[1];
            end if;
        else 
            r1_mod := r1_mod * factor[1]^(factor[2]);
        end if;
    end for;

    //r1_mod is r1 with all all one-variable factors and known factors removed

    //select values to specialize at that don't reduce degree.
    degy:=Degree(r1_mod,y);
    degb:=Degree(r1_mod,b);
	if printing then
		print "Degree with Respect to y: ",degy;
		print "Degree with Respect to b: ",degb;
	end if;

    val_b := 2; 
    while Degree(Evaluate(r1_mod,[y,0,0,0,val_b]),y) ne degy do
        val_b := val_b+1;
    end while;

    val_y := 2; 
    while Degree(Evaluate(r1_mod,[val_y,0,0,0,b]),b) ne degb do
        val_y := val_y+1;
    end while;

    if printing then
        print "Specializing at y=",val_y," and b=",val_b;
    end if;

	//substitute formula for y' obtained from r1 into r2
	//yp := -Derivative(r1,index_b) / Derivative(r1,index_y);

	//or rather, compute the numerator for that relation which 
	//after specializing b

	yp_num := -Derivative(r1,index_b);
	yp_den := Derivative(r1,index_y);


    //substitute for b

    yp_num_b := Evaluate(yp_num,[y,0,0,0,val_b]);
    yp_den_b := Evaluate(yp_den,[y,0,0,0,val_b]);

	r2_coeffyp := Coefficients(r2_mod,index_yp);
	mp := # r2_coeffyp;

	precomputed_powers_num:= [L!1];
	precomputed_powers_den:= [L!1];
	power_num := L!1;
	power_den := L!1;
	counter := 1;
	while counter lt mp do
		power_num := power_num * yp_num_b;
		power_den := power_den * yp_den_b;
		Append(~precomputed_powers_num, power_num);
		Append(~precomputed_powers_den, power_den);

		counter := counter+1;
	end while;

	rel1_b := &+ [Evaluate(r2_coeffyp[j+1],[y,0,0,0,val_b]) * precomputed_powers_num[j+1] * precomputed_powers_den[mp-j] : j in [0..mp-1]];
	
	assert(rel1_b eq Numerator(rel1_b));
	rel1_b := Numerator(rel1_b);
	
    common_b := GCD(rel1_b,Evaluate(r1_mod,[y,0,0,0,val_b]));
	//factors are the potential isogeny correspondeces, specialized



   //substitute for y
    yp_num_y := Evaluate(yp_num,[val_y,0,0,0,b]);
    yp_den_y := Evaluate(yp_den,[val_y,0,0,0,b]);

	r2_coeffyp := Coefficients(r2_mod,index_yp);
	mp := # r2_coeffyp;

	precomputed_powers_num:= [L!1];
	precomputed_powers_den:= [L!1];
	power_num := L!1;
	power_den := L!1;
	counter := 1;
	while counter lt mp do
		power_num := power_num * yp_num_y;
		power_den := power_den * yp_den_y;
		Append(~precomputed_powers_num, power_num);
		Append(~precomputed_powers_den, power_den);

		counter := counter+1;
	end while;
	rel1_y := &+ [Evaluate(r2_coeffyp[j+1],[val_y,0,0,0,b]) * precomputed_powers_num[j+1] * precomputed_powers_den[mp-j] : j in [0..mp-1]];
	
	assert(rel1_y eq Numerator(rel1_y));
	rel1_y := Numerator(rel1_y);

    common_y := GCD(rel1_y,Evaluate(r1_mod,[val_y,0,0,0,b]));
	
	if printing then
		print "Time to Compute Potential Isogeny Corresopndences: ",Cputime(start_time);
	end if;

    //now verify if these relations actually produce isogenies

    check1:=verify_solution_specialize(f1,f2,g1,g2,common_b,1,val_b);
    check2:=verify_solution_specialize(f1,f2,g1,g2,common_y,5,val_y);

    if printing then
        print check1,check2;
        print "Time to Check for Actual Isogeny Corresopndences: ",Cputime(start_time);
    end if;

	//return false if it seems like a real isogeny, true if it is ruled out
	return not (check1 and check2);
end function;