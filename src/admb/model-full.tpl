DATA_SECTION

	init_int nobs
	init_int npop
	init_int nlin
	init_int ninter_pop
	init_int ninter_lin
	
	init_vector pop_mother(1,nobs)
	init_vector pop_father(1,nobs)
	init_vector lin_mother(1,nobs)
	init_vector lin_father(1,nobs)
	init_vector pop_inter(1,nobs)
	init_vector lin_inter(1,nobs)
	
	// For simplicity, model coefficients are provided along with the data
	// It could be recalculated here, but faster/easier/more flexible in R
	init_vector model_a_pop(1,nobs)
	init_vector model_a_lin(1,nobs)
	init_vector model_d_pop(1,nobs)
	init_vector model_d_lin(1,nobs)
	init_vector model_aa_pop(1,nobs)
	init_vector model_aa_lin(1,nobs)
	init_vector model_ad_pop(1,nobs)
	init_vector model_ad_lin(1,nobs)
	init_vector model_dd_pop(1,nobs)
	init_vector model_dd_lin(1,nobs)
	
	init_vector trait(1,nobs)
	
PARAMETER_SECTION

	init_number mu(1)
	init_number mean_d_pop(0)
	init_number mean_d_lin(0)
	init_number mean_e_pop(0)
	init_number mean_e_lin(0)
	init_bounded_number sigma_a_pop(0.000001,1000000,1)
	init_bounded_number sigma_a_lin(0.000001,1000000,1)
	init_bounded_number sigma_d_pop(0.000001,1000000,0)
	init_bounded_number sigma_d_lin(0.000001,1000000,0)
	init_bounded_number sigma_e_pop(0.000001,1000000,0)
	init_bounded_number sigma_e_lin(0.000001,1000000,0)
	init_bounded_number sigma_R(0.000001,1000000,1)
	
	vector pred_trait(1,nobs)

	random_effects_vector u_a_pop(1,npop,1)
	random_effects_vector u_a_lin(1,nlin,1)
	random_effects_vector u_d_pop(1,ninter_pop,0)
	random_effects_vector u_d_lin(1,ninter_lin,0)
	random_effects_vector u_e_pop(1,ninter_pop,0)
	random_effects_vector u_e_lin(1,ninter_lin,0)
	
	objective_function_value f
	
PROCEDURE_SECTION
	
	f = 0;

	// Prior part for random effects
	f -= - 0.5*norm2(u_a_pop);
	f -= - 0.5*norm2(u_a_lin);
	f -= - 0.5*norm2(u_d_pop);
	f -= - 0.5*norm2(u_d_lin);
	f -= - 0.5*norm2(u_e_pop);
	f -= - 0.5*norm2(u_e_lin);
	
	// Predicted phenotypes
	for (int i = 1; i <= nobs; i++) {
		pred_trait[i] = mu
		                + model_a_pop(i) *  sigma_a_pop*u_a_pop(pop_mother(i)) 
		                + model_a_pop(i) *  sigma_a_pop*u_a_pop(pop_father(i))
		                + model_a_lin(i) *  sigma_a_lin*u_a_lin(lin_mother(i))
		                + model_a_lin(i) *  sigma_a_lin*u_a_lin(lin_father(i))
		                + model_d_pop(i) * (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))
		                + model_d_lin(i) * (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i)))
		                + model_aa_pop(i)*  sigma_a_pop*u_a_pop(pop_mother(i))
		                                 *  sigma_a_pop*u_a_pop(pop_father(i)) 
		                                 * (mean_e_pop + sigma_e_pop*u_e_pop(pop_inter(i)))
		                + model_aa_lin(i)*  sigma_a_lin*u_a_lin(lin_mother(i))
		                                 *  sigma_a_lin*u_a_lin(lin_father(i)) 
		                                 * (mean_e_lin + sigma_e_lin*u_e_lin(lin_inter(i)))
		                + model_ad_pop(i)*  sigma_a_pop*(u_a_pop(pop_mother(i))+u_a_pop(pop_father(i)))
		                                 * (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))
		                                 * (mean_e_pop + sigma_e_pop*u_e_pop(pop_inter(i)))
		                + model_ad_lin(i)*  sigma_a_lin*(u_a_lin(lin_mother(i))+u_a_lin(lin_father(i)))
		                                 * (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i)))
		                                 * (mean_e_lin + sigma_e_lin*u_e_lin(lin_inter(i)))
		                + model_dd_pop(i)* (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))
		                                 * (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))
		                                 * (mean_e_pop + sigma_e_pop*u_e_pop(pop_inter(i)))
		                + model_dd_lin(i)* (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i))) 
		                                 * (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i)))
		                                 * (mean_e_lin + sigma_e_lin*u_e_lin(lin_inter(i)));
	}
	
	f -= -nobs*log(sigma_R) - 0.5*norm2((pred_trait-trait)/sigma_R);

TOP_OF_MAIN_SECTION

	arrmblsize = 400000000L;
	// gradient_structure::set_ARRAY_MEMBLOCK_SIZE(200000L);
	gradient_structure::set_MAX_NVAR_OFFSET(100000);
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000000L);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000L);
	
