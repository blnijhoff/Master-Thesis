function solve_and_summarize_model_with_M_scenarios(
	sp_dir::String,
	M::Integer,
	time_limit_sec::Float64,
	risk_measure_k::Integer = 4
)
	risk_measure_alpha_values = [(i - 1)/risk_measure_k for i in 1:risk_measure_k]
	sp, scenario_elements = read_SPMS(sp_dir)
	M_scenarios = sample_from_scenario_dimensions(scenario_elements, M)
	model = initialize_model_with_scenarios(
		sp,
		M_scenarios,
		risk_measure_alpha_values,
		time_limit_sec
	)
	# println(model)
	optimize!(model)
	# println(value(model[:x]))
	# println(value(model[:y]))
	# println(value(model[:z]))
	# println(value(model[:u]))
	assert_is_solved_and_feasible(model)
	solution_summary(model)
end