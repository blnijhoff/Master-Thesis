abstract type abstract_scenario_generation_method end
struct sampling_benchmark_method <: abstract_scenario_generation_method end
struct method_1 <: abstract_scenario_generation_method
	time_limit_sec::Float64
	num_of_trivial_samples::Integer
	size_of_trivial_samples::Integer
	max_clustering_iters::Integer
end
struct method_2 <: abstract_scenario_generation_method
	time_limit_sec::Float64
	size_of_trivial_sample::Integer
end
struct method_3 <: abstract_scenario_generation_method
	time_limit_sec::Float64
	size_of_trivial_sample::Integer
end
struct method_4 <: abstract_scenario_generation_method
	time_limit_sec::Float64
	num_of_trivial_samples::Integer
	size_of_trivial_samples::Integer
	p::Float64
end

function invoke_method(_::abstract_scenario_generation_method,
	sp::sp_instance,
	risk_measure_alpha_values::Vector{Float64},
	original_scenarios::Vector{scenario_vec},
	number_to_generate::Integer
)::Vector{scenario_vec}
	error("Not implemented")
end

function invoke_method(_::sampling_benchmark_method,
	sp::sp_instance,
	risk_measure_alpha_values::Vector{Float64},
	original_scenarios::Vector{scenario_vec},
	number_to_generate::Integer
)::Vector{scenario_vec}
	return sample_m_scenarios_from_M(
		original_scenarios,
		number_to_generate
	)
end

function invoke_method(m::method_1,
	sp::sp_instance,
	risk_measure_alpha_values::Vector{Float64},
	original_scenarios::Vector{scenario_vec},
	number_to_generate::Integer
)::Vector{scenario_vec}
	# Use k-medoids clustering with k = m.
	# Create M x M distance matrix between each scenario
	# d(s_1, s_2) = || [g_x(s_1)]_x in X - [g_x(s_2)]_x in X ||
	# Sample m trivially from M (set S) multiple times
	class_of_trivial_samples = [
		sample_m_scenarios_from_M(original_scenarios, m.size_of_trivial_samples)
		for _ in 1:m.num_of_trivial_samples
	]
	# Solve model to give x_hat ∈ argmin ρ_S [fₓ] for each S in sample sets
	candidate_sols = [
		solve_model_with_scenarios_safe(
			sp,
			S,
			risk_measure_alpha_values,
			m.time_limit_sec
		)[1]
		for S in class_of_trivial_samples
	]
	# Calculate f(x, ξ) for every ξ in M scenarios for every x in candidate_sols
	function_values_for_each_scenario_and_candidate_sol::Vector{
		Vector{Tuple{Float64, scenario_vec}}
	} = [
		collect_function_values(sp, m.time_limit_sec, original_scenarios, x) for
		x in candidate_sols
	]
	if length(function_values_for_each_scenario_and_candidate_sol) == 0
		error()
	end

	# Determine g_x = ζ_x * f_x on M scenarios for every x in candidate_sols
	g_dicts::Vector{Dict{Vector{Float64}, Float64}} = []
	for function_values_x in function_values_for_each_scenario_and_candidate_sol
		ζ_dict = ζ(risk_measure_alpha_values, function_values_x)
		g_dict = Dict{Vector{Float64}, Float64}()
		for (f, scenario) in function_values_x
			g_dict[scenario.vec] = ζ_dict[scenario.vec] * f
		end
		push!(g_dicts, g_dict)
	end

	M::Integer = length(original_scenarios)

	dist::Matrix{Float64} = zeros(Float64, M, M) # d(i,i) = 0 and d(i,j) = d(j, i)

	# Iterate over i < j
	for i in 1:(M-1)
		for j in (i+1):M
			s_i = original_scenarios[i].vec
			s_j = original_scenarios[j].vec
			vector_i = [g_dict[s_i] for g_dict in g_dicts] # [g_x(s_i)]_x in X
			vector_j = [g_dict[s_j] for g_dict in g_dicts] # [g_x(s_j)]_x in X
			distance = norm(vector_i - vector_j) # 2-norm of vector_i - vector_j
			dist[i, j] = distance
			dist[j, i] = distance
		end
	end
	# Set all Inf values to the maximum distance without inf
	max_value_without_inf = maximum(filter(x -> isfinite(x), dist))
	for idx in eachindex(dist)
		if !isfinite(dist[idx])
			dist[idx] = max_value_without_inf
		end
	end

	# Cluster with k-medoids
	clustering_result = kmedoids(dist, number_to_generate; maxiter = m.max_clustering_iters)
	# Find the generated scenarios
	generated_scenarios_dict = Dict{Vector{Float64}, Float64}()
	for i in 1:M
		medoid_index = clustering_result.medoids[clustering_result.assignments[i]]
		medoid_scenario = original_scenarios[medoid_index].vec
		generated_scenarios_dict[medoid_scenario] =
			get(generated_scenarios_dict, medoid_scenario, 0) + original_scenarios[i].prob
	end
	return [
		scenario_vec(vec, prob) for (vec, prob) in generated_scenarios_dict
	]
end

# Algorithm 1 with arpon as CVaR method.
function invoke_method(m::method_2,
	sp::sp_instance,
	risk_measure_alpha_values::Vector{Float64},
	original_scenarios::Vector{scenario_vec},
	number_to_generate::Integer
)::Vector{scenario_vec}
	# Sample m trivially from M (set S)
	trivial_sample = sample_m_scenarios_from_M(original_scenarios, m.size_of_trivial_sample)
	# Solve model to give x_hat ∈ argmin ρ_S [fₓ]
	candidate_sol, _ = solve_model_with_scenarios_safe(
		sp,
		trivial_sample,
		risk_measure_alpha_values,
		m.time_limit_sec
	)
	# Determine the qunatiles of a guess on the optimal solution!
	function_values =
		collect_function_values(sp, m.time_limit_sec, original_scenarios, candidate_sol)
	# sort function values
	sort!(function_values, by = x->x[1])
	# Take the cumulative sums of probabilities
	cumulative_probs = cumsum([tup[2].prob for tup in function_values])
	# A scenario is included in CVaR_α if the value is above or at the α-quantile,
	# i.e., when the cdf at that value is ≥ α (or the cumulative_probs equivalently)
	# These are thus given a weight of ∑ μ_j for j s.t. α_j ≤ α = cum_prob
	# First a dict: ξ → w_s
	weighting_dict = Dict{Vector{Float64}, Float64}(
		[
		(
			s.vec,
			count(alpha -> alpha <= cum_prob, risk_measure_alpha_values) /
			length(risk_measure_alpha_values)
		) for
		(cum_prob, (_, s)) in zip(cumulative_probs, function_values)
	]
	)
	weighted_scenarios::Vector{scenario_vec} =
		[scenario_vec(s.vec, s.prob * weighting_dict[s.vec]) for s in original_scenarios]
	# normalize scenarios to have probability 1
	weighted_scenarios = normalize_probs(weighted_scenarios)

	# Sample from weighted and normalized scenario probabilities
	weighted_sample = sample_m_scenarios_from_M(
		weighted_scenarios,
		number_to_generate
	)
	# Divide probabilities by weights and normalize
	return normalize_probs(
		perturb_scenario_probs_with_inverse_dict(weighted_sample, weighting_dict)
	)
end

function invoke_method(m::method_3,
	sp::sp_instance,
	risk_measure_alpha_values::Vector{Float64},
	original_scenarios::Vector{scenario_vec},
	number_to_generate::Integer
)::Vector{scenario_vec}
	return method_3_and_4_inner(
		sp,
		risk_measure_alpha_values,
		original_scenarios,
		number_to_generate,
		m.time_limit_sec,
		1,
		m.size_of_trivial_sample,
		1.0
	)
end

function invoke_method(m::method_4,
	sp::sp_instance,
	risk_measure_alpha_values::Vector{Float64},
	original_scenarios::Vector{scenario_vec},
	number_to_generate::Integer
)::Vector{scenario_vec}
	return method_3_and_4_inner(
		sp,
		risk_measure_alpha_values,
		original_scenarios,
		number_to_generate,
		m.time_limit_sec,
		m.num_of_trivial_samples,
		m.size_of_trivial_samples,
		m.p
	)
end

function method_3_and_4_inner(
	sp::sp_instance,
	risk_measure_alpha_values::Vector{Float64},
	original_scenarios::Vector{scenario_vec},
	number_to_generate::Integer,
	time_limit_sec::Float64,
	trivial_scenario_sets_to_generate::Integer,
	size_of_trivial_scenario_sets::Integer,
	p::Float64
)::Vector{scenario_vec}
	# Sample m trivially from M (set S) multiple times
	class_of_trivial_samples = [
		sample_m_scenarios_from_M(original_scenarios, size_of_trivial_scenario_sets)
		for _ in 1:trivial_scenario_sets_to_generate
	]
	# Solve model to give x_hat ∈ argmin ρ_S [fₓ] for each S in sample sets
	candidate_sols = [
		solve_model_with_scenarios_safe(sp, S, risk_measure_alpha_values, time_limit_sec)[1]
		for S in class_of_trivial_samples
	]
	# Calculate f(x, ξ) for every ξ in M scenarios for every x in candidate_sols
	function_values_for_each_scenario_and_candidate_sol::Vector{
		Vector{Tuple{Float64, scenario_vec}}
	} = [
		collect_function_values(sp, time_limit_sec, original_scenarios, x) for
		x in candidate_sols
	]
	if length(function_values_for_each_scenario_and_candidate_sol) == 0
		error()
	end
	# Determine ζ_x on M scenarios for every x in candidate_sols
	ζ_dicts = [
		ζ(risk_measure_alpha_values, function_values_x)
		for function_values_x in function_values_for_each_scenario_and_candidate_sol
	]

	# Determine perturbation dict: ζ_dicts[1] if only one, otherwise p-moment
	perturb_dict =
		trivial_scenario_sets_to_generate == 1 ?
		ζ_dicts[1] : p_moment_of_zetas(ζ_dicts, p)

	# Perturb scenario probabilities of M scenarios with perturb_dict
	# Normalize scenario probabilities afterwards!
	perturbed_original_scenarios = normalize_probs(
		perturb_scenario_probs_with_dict(
			original_scenarios,
			perturb_dict
		)
	)
	# Generate m scenarios from perturbed probabilities
	m_generated_scenarios =
		sample_m_scenarios_from_M(perturbed_original_scenarios, number_to_generate)
	# perturb m_scenarios with 1 / perturb_dict and normalize!
	return normalize_probs(
		perturb_scenario_probs_with_inverse_dict(
			m_generated_scenarios,
			perturb_dict
		)
	)
end

function collect_function_values(
	sp::sp_instance,
	time_limit_sec::Float64,
	original_scenarios::Vector{scenario_vec},
	x::Vector{Float64}
)::Vector{Tuple{Float64, scenario_vec}}
	return collect(
		zip(
			compute_outcome_values_for_each_scenario_and_fixed_x(
				sp,
				time_limit_sec,
				original_scenarios,
				x
			), original_scenarios
		)
	)
end

# Create a sample of predefined size, not a predefined number of samples!
function sample_m_scenarios_from_M(
	M_scenarios::Vector{scenario_vec},
	m::Integer
)::Vector{scenario_vec}
	sample_counts = Dict{Vector{Float64}, Float64}()
	scenarios = getproperty.(M_scenarios, :vec)
	probs = getproperty.(M_scenarios, :prob)

	# m can be at most the length of M_scenarios!
	m = min(m, length(M_scenarios))

	sampled = 0
	while length(sample_counts) < m
		sampled_scenario = wsample(scenarios, probs)
		sample_counts[sampled_scenario] = get(sample_counts, sampled_scenario, 0) + 1
		sampled += 1
	end

	return [
		scenario_vec(vec, count / sampled) for (vec, count) in sample_counts
	]
end

function sample_from_scenario_dimensions(
	scenario_elements::Vector{Vector{scenario_elem}},
	M::Integer
)::Vector{scenario_vec}
	sample_counts = Dict{Vector{Float64}, Float64}()
	sampled = 0

	while length(sample_counts) < M
		sampled_scenario = [
			wsample(getproperty.(vec, :elem), getproperty.(vec, :prob))
			for vec in scenario_elements
		]
		sample_counts[sampled_scenario] = get(sample_counts, sampled_scenario, 0) + 1
		sampled += 1
	end

	return [
		scenario_vec(vec, count / sampled) for (vec, count) in sample_counts
	]
end

function perturb_scenario_probs_with_dict_and_normalize(
	scenarios::Vector{scenario_vec},
	dict::Dict
)::Vector{scenario_vec}
	return normalize_probs(
		perturb_scenario_probs_with_dict(
			scenarios,
			dict
		)
	)
end

function perturb_scenario_probs_with_inverse_dict_and_normalize(
	scenarios::Vector{scenario_vec},
	dict::Dict
)::Vector{scenario_vec}
	return normalize_probs(
		perturb_scenario_probs_with_inverse_dict(
			scenarios,
			dict
		)
	)
end

function perturb_scenario_probs_with_dict(
	scenarios::Vector{scenario_vec},
	dict::Dict
)::Vector{scenario_vec}
	return perturb_scenario_probs(
		scenarios,
		(v, d) -> d[v],
		dict
	)
end

function perturb_scenario_probs_with_inverse_dict(
	scenarios::Vector{scenario_vec},
	dict::Dict
)::Vector{scenario_vec}
	return perturb_scenario_probs(
		scenarios,
		(v, d) -> 1 / d[v],
		dict
	)
end

function perturb_scenario_probs(
	scenarios::Vector{scenario_vec},
	perturbation_lambda::Function,
	dict_in_lambda::Dict
)::Vector{scenario_vec}
	return [
		scenario_vec(
			scenario.vec,
			scenario.prob * perturbation_lambda(scenario.vec, dict_in_lambda)
		) for
		scenario in scenarios
		if perturbation_lambda(scenario.vec, dict_in_lambda) > 0
	]
end

function normalize_probs(scenarios::Vector{scenario_vec})::Vector{scenario_vec}
	total_probs = sum(getproperty.(scenarios, :prob))
	return [
		scenario_vec(scenario.vec, scenario.prob / total_probs) for scenario in scenarios
	]
end

function w_rho(
	risk_measure_alpha_values::Vector{Float64},
	alpha::Float64
)::Float64
	# Clamp alpha value to [0,1]
	# (can be just above 1 because of usage with floating point division in probabilities)
	alpha = clamp(alpha, 0, 1)

	k = length(risk_measure_alpha_values)

	return sum(
		(alpha - alpha_i) / (k * (1 - alpha_i)) for
		alpha_i in risk_measure_alpha_values if alpha_i <= alpha
	)
end

function ζ(
	risk_measure_alpha_values::Vector{Float64},
	function_values::Vector{Tuple{Float64, scenario_vec}} # Tuples (f, (ξ, p))
)::Dict{Vector{Float64}, Float64} # Dict ξ => ζ

	# Sort function_values
	sort!(function_values, by = x->x[1])
	# cumulative sum of probabilities for sorted function_values
	cumulative_probs = cumsum([tup[2].prob for tup in function_values])
	w_vals = [w_rho(risk_measure_alpha_values, p) for p in cumulative_probs]

	ζ_dict = Dict{Vector{Float64}, Float64}()
	# Compute zeta[i] with (cumulative_w_vals[i] - (cumulative_w_vals[i-1] or 0 if i = 1)) / p[i]
	for i in eachindex(function_values)
		ζ_dict[function_values[i][2].vec] =
			(w_vals[i] - (i == 1 ? 0 : w_vals[i-1])) /
			function_values[i][2].prob
	end
	return ζ_dict
end

function p_moment_of_zetas(
	ζ_dicts::Vector{Dict{Vector{Float64}, Float64}}, # List [Dict ξ => ζ_x for different x]
	p::Float64
)::Dict{Vector{Float64}, Float64} # Dict ξ => p-moment of ζ_x for all given x
	@assert(p >= 1)
	@assert(length(ζ_dicts) >= 1)
	@assert(all(length(ζ_dicts[1]) == length(ζ_dicts[i]) for i in 2:length(ζ_dicts)))

	factor = (1 / length(ζ_dicts)) ^ (1 / p)

	return Dict{Vector{Float64}, Float64}([
		(ξ, factor * norm([ζ_dict[ξ] for ζ_dict in ζ_dicts], p)) for
		ξ in keys(ζ_dicts[1])
	])
end
