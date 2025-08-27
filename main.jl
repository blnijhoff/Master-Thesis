using Random,
	StatsPlots,
	LinearAlgebra,
	StatsPlots,
	HiGHS,
	JuMP,
	Clustering,
	StatsBase,
	CSV,
	DataFrames,
	HiGHS,
	LaTeXStrings,
	Dates,
	Glob

include("scenarios.jl")
include("stochastic_program.jl")
include("scenario_generation.jl")
include("tests.jl")
include("extra-plots.jl")

function perform_all_tests()
	stochastic_programs_and_M_values::Vector{Tuple{String, Integer}} = [
		("simple-recourse", 10)
		("LandS", 5000)
		("20term", 500)
		("gbd", 5000)
		("ssn", 250)
		("storm", 500)
	]
	for (sp_dir, M) in stochastic_programs_and_M_values
		start_time = now()
		println("Starting $sp_dir with $M")
		perform_full_test_for_one_sp(sp_dir, M)
		end_time = now()
		open("total_runtimes.txt", "a") do f
			println(f, "$end_time: $sp_dir completed full test in $(end_time - start_time)")
		end
	end
	return nothing
end
# Use environment once
# const env = Gurobi.Env()
"""
# perform_full_test_for_one_sp
Performs `perform_tests` for all possible parameter sets given as below: \\
	- k ∈ {2, 3, 4} \\
	- n₁ = 10 \\
	- n₂ = 1 \\
	- M given. \\
	- m = M / 10 \\
	- time_limit_sec = 600. \\
It stores *all* results. \\
M is first maximized as the product of the length of the scenario dimensions. 
Then m is taken to be 25 % of that.
"""
function perform_full_test_for_one_sp(sp_dir::String, M::Integer)
	start_time = now()
	n_1::Integer = 10
	n_2::Integer = 1

	# Determine actual value of M
	_, scenario_elements = read_SPMS(sp_dir)
	# M can be at most the product of the dimensions of the one-dimension independent scenarios
	# Since this product can explode, do a running product and check if the product exceeds M
	running_product = 1
	for vec in scenario_elements
		if M < running_product
			break
		end
		running_product *= length(vec)
	end
	# If broken out early M < running_product, so minimum = M
	M = min(M, running_product)

	# m is M / 10 (rounded)
	m::Integer = div(M, 10)
	time_limit_sec::Float64 = 600.0

	k_values::Vector{Integer} = [2, 3, 4]

	for k in k_values
		perform_tests(
			sp_dir,
			k,
			M, m, n_1, n_2;
			store_results_values = true,
			store_results_boxplots = true,
			store_results_individual_values = true,
			store_results_individual_boxplots = true,
			time_limit_sec = time_limit_sec
		)
	end

	end_time = now()
	println("Test for $sp_dir completed in $(end_time - start_time).")
end

function perform_tests(
	sp_dir::String,
	risk_measure_k::Integer,
	M::Integer, m::Integer, n_1::Integer, n_2::Integer;
	store_results_values::Bool = true,
	store_results_boxplots::Bool = true,
	store_results_individual_values::Bool = false,
	store_results_individual_boxplots::Bool = false,
	seed::Integer = 40895238092384309,
	time_limit_sec::Float64 = 60.0
)

	risk_measure_alpha_values = [(i - 1)/risk_measure_k for i in 1:risk_measure_k]

	# Size of trivial samples need to be at most the size M!
	num_of_trivial_samples::Integer = 5
	size_of_trivial_samples::Integer = min(M, 10)

	sp, scenario_elements = read_SPMS(sp_dir)



	sg_methods::Vector{Tuple{String, abstract_scenario_generation_method}} = [
		("Benchmark", sampling_benchmark_method()),
		(
			"Method 1",
			method_1(time_limit_sec, num_of_trivial_samples, size_of_trivial_samples, 10)
		),
		(
			"Method 2",
			method_2(time_limit_sec, size_of_trivial_samples)
		),
		(
			"Method 3",
			method_3(time_limit_sec, size_of_trivial_samples)
		),
		(
			"Method 4",
			method_4(time_limit_sec, num_of_trivial_samples, size_of_trivial_samples, 2.0)
		)
	]

	G_values = nothing
	G_SP_values = nothing
	if store_results_values || store_results_boxplots
		G_values = Vector{Vector{Float64}}(undef, length(sg_methods))
		G_SP_values = Vector{Vector{Float64}}(undef, length(sg_methods))
	end

	for (i, (sg_display_name, sg_method)) in enumerate(sg_methods)
		println("Performing $sg_display_name")
		G, G_SP = perform_test(
			sp,
			scenario_elements,
			risk_measure_alpha_values,
			sg_display_name,
			sg_method,
			M, m, n_1, n_2,
			store_results_individual_boxplots,
			store_results_individual_values,
			seed,
			time_limit_sec
		)
		if store_results_values || store_results_boxplots
			G_values[i] = G
			G_SP_values[i] = G_SP
		end
	end

	if store_results_values || store_results_boxplots
		risk_measure_string = reduce((x, y)->"$x-$y", risk_measure_alpha_values)
		file_name = "_$(sp.name)_rho$(risk_measure_string)_M-$(M)_m-$(m)_n1-$(n_1)_n2-$(n_2)"
		if store_results_values
			headers = [tuple[1] for tuple in sg_methods]
			CSV.write(
				"output/G$file_name.csv",
				DataFrame([col => G_values[i] for (i, col) in enumerate(headers)])
			)
			CSV.write(
				"output/GSP$file_name.csv",
				DataFrame([col => G_SP_values[i] for (i, col) in enumerate(headers)])
			)
		end

		if store_results_boxplots
			make_box_plots(
				G_values,
				[name for (name, _) in sg_methods],
				"Scenario Generation Method",
				L"G",
				"Boxplots of " * string(L"G\;") *
				"for different scenario generation methods",
				save_file = true,
				file_name = "output/G$file_name.png"
			)

			make_box_plots(
				G_SP_values,
				[name for (name, _) in sg_methods],
				"Scenario Generation Method",
				L"G_{\textrm{SP}}",
				"Boxplots of " * string(L"G_{\mathrm{SP}}\;") *
				"for different scenario generation methods",
				save_file = true,
				file_name = "output/GSP$file_name.png"
			)
		end
	end
end

function perform_test(
	sp::sp_instance,
	scenario_elements::Vector{Vector{scenario_elem}},
	risk_measure_alpha_values::Vector{Float64},
	sg_display_name::String,
	sg_method::abstract_scenario_generation_method, # sg interface
	M::Integer,
	m::Integer,
	n_1::Integer,
	n_2::Integer,
	store_box_plot::Bool,
	store_values::Bool,
	seed::Integer,
	time_limit_sec::Float64
)::NTuple{2, Vector{Float64}}

	# Reset seed of rng
	Random.seed!(seed)

	# True optimality gaps (as vector)
	G::Vector{Float64} = []
	# SP solution value gaps
	G_SP::Vector{Float64} = []

	# Sample all M scenarios first, 
	# so that they coincide for each test with the same seed
	M_scenario_samples::Vector{Vector{scenario_vec}} =
		[sample_from_scenario_dimensions(scenario_elements, M) for _ in 1:n_1]

	for i in 1:n_1
		# Take M samples from ξ, with each dimension independently 
		# (assuming the scenarios are the cartesian product of the individual supports!)
		M_scenarios::Vector{scenario_vec} = M_scenario_samples[i]

		m_sols::Vector{Vector{Float64}} = []
		m_objs::Vector{Float64} = []


		for j in 1:n_2
			println("(i, j) = ($i, $j)")
			# Generate m scenarios from the M scenarios using scenario_generation_method
			m_generated_scenarios = invoke_method(sg_method,
				sp,
				risk_measure_alpha_values,
				M_scenarios,
				m
			)

			# Solve m scenario SP, giving solution x_m and optimal objective v_m
			x_m, v_m = solve_model_with_scenarios_safe(
				sp,
				m_generated_scenarios,
				risk_measure_alpha_values,
				time_limit_sec
			)

			push!(m_sols, x_m)
			push!(m_objs, v_m)
		end

		_, v_M = solve_model_with_scenarios_safe(
			sp,
			M_scenarios,
			risk_measure_alpha_values,
			time_limit_sec
		)

		for j in 1:n_2
			objective_M, model_solved = solve_model_with_scenarios_and_fixed_x(
				sp,
				M_scenarios,
				risk_measure_alpha_values,
				time_limit_sec,
				m_sols[j]
			)
			# If model did not solve, the given x_m is not feasible for the second-stage
			# Thus we do not save the objective
			if model_solved
				push!(G, (objective_M - v_M) / v_M)
			end
			push!(G_SP, (v_M - m_objs[j]) / v_M)

		end
	end

	println("G:")
	display(G)
	println("G_SP:")
	display(G_SP)

	if store_box_plot || store_values
		risk_measure_string = reduce((x, y)->"$x-$y", risk_measure_alpha_values)

		file_name = "_$(sg_display_name)_$(sp.name)_rho$(risk_measure_string)_M-$(M)_m-$(m)_n1-$(n_1)_n2-$(n_2)"

		if store_box_plot
			make_box_plot(
				G,
				L"G",
				"Boxplot of " * string(L"G\;") * "for $sg_display_name",
				save_file = store_box_plot,
				file_name = "output-individual/G$file_name.png",
				x_name = "$(sg_display_name)"
			)

			make_box_plot(
				G_SP,
				L"G_{\textrm{SP}}",
				"Boxplot of " * string(L"G_{\mathrm{SP}}\;") * "for $sg_display_name",
				save_file = store_box_plot,
				file_name = "output-individual/GSP$file_name.png",
				x_name = "$(sg_display_name)"
			)
		end

		if store_values
			open("output-individual/G$file_name.txt", "w") do io
				for val in G
					println(io, val)
				end
			end

			open("output-individual/GSP$file_name.txt", "w") do io
				for val in G_SP
					println(io, val)
				end
			end
		end
	end

	return (G, G_SP)

	# TODO: Store G and G_SP to make a collection of box plots
end

function make_box_plot(
	values::Vector{Float64},
	y_label,
	plot_name;
	save_file::Bool = false,
	file_name::String = "",
	x_name::String = "",
	x_label::String = "Scenario Generation Method"
)
	make_box_plots([values], [x_name], x_label, y_label, plot_name; save_file, file_name)
end

function make_box_plots(
	vectors_of_values::Vector{Vector{Float64}},
	x_names::Vector{String},
	x_label,
	y_label,
	plot_name;
	save_file::Bool = false,
	file_name::String = ""
)
	@assert(length(x_names) == length(vectors_of_values))
	gr(size = (length(x_names) * 200, 600))

	# Size of groups needs to be the same as size of values!
	groups = Vector{String}()
	values = Vector{Float64}()
	for (name, vector) in zip(x_names, vectors_of_values)
		append!(groups, fill(name, length(vector)))
		append!(values, vector)
	end

	# Plot
	plot = boxplot(
		groups,
		values,
		legend = false,
		title = plot_name,
		xlabel = x_label,
		ylabel = y_label,
		# fillalpha = 0.7,
		linewidth = 2,
		markerstrokewidth = 2,
		fillcolor = "#a6cee3",
		markercolor = "#FFFFFF",
		tickfontfamily = "cmunbx",
		tickfontsize = 11,
		titlefontfamily = "cmunbx",
		# titlefontsize = 20,
		gridlinewidth = 0.5,
		gridalpha = 0.5,
		guidefontfamily = "cmunrm",
		guidefontsize = 14,
		# gridstyle = :,
		background_color = :transparent,
		grid = :all,
		margin = (10.0, :mm),
		dpi = 400
	)
	# Displays the plot if it is not saved
	if save_file
		savefig(plot, file_name)
	else
		display(plot)
	end
end


# perform_all_tests()