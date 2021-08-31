DATA_SECTION

	init_int nobs
	init_int npop
	init_vector pop_mother(1,nobs)
	init_vector pop_father(1,nobs)
	init_vector trait(1,nobs)
	
PARAMETER_SECTION

	init_number mu
	init_bounded_number sigma_a_pop(0.000001,1000000,1)
	init_bounded_number sigma_R(0.000001,1000000,1)
	
	vector pred_trait(1,nobs)

	random_effects_vector u_a_pop(1,npop)
	
	objective_function_value f
	
PROCEDURE_SECTION
	
	f = 0;

	// Prior part for random effects
	f -= - 0.5*norm2(u_a_pop);
	
	// Predicted phenotypes
	for (int i = 1; i <= nobs; i++) {
		pred_trait[i] = mu 
		                + sigma_a_pop*u_a_pop(pop_mother(i)) 
		                + sigma_a_pop*u_a_pop(pop_father(i));
	}
	
	f -= -nobs*log(sigma_R) - 0.5*norm2((pred_trait-trait)/sigma_R);
	

	
