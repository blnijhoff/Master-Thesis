@enum constraint_sense begin
	equal
	greater_or_equal
	less_or_equal
end

"""
Stochastic program of the form:
min c_1'x + ρ[c_2'y]
s.t. A_1x ? b_1
	A_2x + A_3y(ξ) ? b_2 + ξ, ∀ξ
	[x ≤ dₗ] (optional)
	[x ≥ dᵤ] (optional)
	[y(ξ) ≤ d¹ₗ, ∀ξ] (optional)
	[y(ξ) ≥ d¹ᵤ, ∀ξ] (optional)

This is a 2 stage SP recourse model.  
Each row of b_1 and b_2 can have a different `sense`, i.e.,
equality, greater than or equal, or less than or equal.
The support of ξ is not stored in this struct, as it differs every time
"""
struct sp_instance
	name::String
	first_stage_cols::Integer
	second_stage_cols::Integer
	first_stage_rows::Integer
	second_stage_rows::Integer # Also the dimension of ξ
	c_1::Vector{Float64}
	c_2::Vector{Float64}
	A_1::Matrix{Float64}
	A_2::Matrix{Float64}
	A_3::Matrix{Float64}
	senses_1::Vector{constraint_sense}
	senses_2::Vector{constraint_sense}
	b_1::Vector{Float64}
	b_2::Vector{Float64}
	d_l::Union{Vector{Float64}, Nothing}
	d_u::Union{Vector{Float64}, Nothing}
	d1_l::Union{Vector{Float64}, Nothing}
	d1_u::Union{Vector{Float64}, Nothing}
end

struct row_info
	sense::constraint_sense
	name::String
end

"""
directory_name must be the same as the different files in the directory!
This assumes it is a 2-stage SP!
"""
function read_SPMS(
	directory_name::String
)::Tuple{sp_instance, Vector{Vector{scenario_elem}}}
	cor_file_name = joinpath("SPMS-instances", directory_name, "$directory_name.cor")
	tim_file_name = joinpath("SPMS-instances", directory_name, "$directory_name.tim")
	sto_file_name = joinpath("SPMS-instances", directory_name, "$directory_name.sto")
	@assert(isfile(cor_file_name) && isfile(tim_file_name) && isfile(sto_file_name))

	first_stage_col, first_stage_row, second_stage_col, second_stage_row =
		ntuple(_ -> "", 4)


	open(tim_file_name, "r") do io
		status::String = ""
		stages_seen::Integer = 0
		for line in eachline(io)
			if status == ""
				if startswith(line, "TIME")
					status = "TIME"
				end
				continue
			end
			if status == "TIME"
				if startswith(line, "PERIODS")
					status = "PERIODS"
				end
				continue
			end
			if status == "PERIODS"
				if !startswith(line, "ENDATA")
					split_line = split(line)
					if stages_seen == 0
						first_stage_col = split_line[1]
						first_stage_row = split_line[2]
					elseif stages_seen == 1
						second_stage_col = split_line[1]
						second_stage_row = split_line[2]
					end
					stages_seen += 1
					continue
				end
				break
			end
		end
	end

	name::String = ""
	c_1, c_2, b_1, b_2 = ntuple(_ -> Vector{Float64}(), 4)
	senses_1, senses_2 = ntuple(_ -> Vector{constraint_sense}(), 2)
	A_1, A_2, A_3 = ntuple(_ -> Matrix{Float64}(undef, 0, 0), 3)
	d_l, d_u, d1_l, d1_u = ntuple(_ -> nothing, 4)

	first_stage_row_names::Vector{String} = []
	first_stage_rows::Vector{row_info} = []
	second_stage_row_names::Vector{String} = []
	second_stage_rows::Vector{row_info} = []
	first_stage_col_names::Vector{String} = []
	second_stage_col_names::Vector{String} = []

	open(cor_file_name, "r") do io
		status::String = ""
		row_status::String = ""
		column_status::String = ""
		objective_name = ""

		current_column_first_stage_rows = []
		current_column_second_stage_rows = []
		for line in eachline(io)
			if status == ""
				if startswith(line, "NAME")
					status = "NAME"
					name = split(line)[2]
				end
				continue
			end

			if status == "NAME"
				if startswith(line, "ROWS")
					status = "ROWS"
				end
				continue
			end

			if status == "ROWS"
				if startswith(line, "COLUMNS")
					status = "COLUMNS"
					b_1 = zeros(Float64, length(first_stage_rows))
					b_2 = zeros(Float64, length(second_stage_rows))
					A_1 = Matrix{Float64}(undef, length(first_stage_rows), 0)
					A_2 = Matrix{Float64}(undef, length(second_stage_rows), 0)
					A_3 = Matrix{Float64}(undef, length(second_stage_rows), 0)
					current_column_first_stage_rows =
						zeros(Float64, length(first_stage_rows))
					current_column_second_stage_rows =
						zeros(Float64, length(second_stage_rows))
					continue
				end

				split_line = split(line)

				if split_line[1] == "N"
					objective_name = split_line[2]
					continue
				end

				if split_line[2] == first_stage_row
					row_status = "FIRST"
				elseif split_line[2] == second_stage_row
					row_status = "SECOND"
				end

				sense = nothing
				if split_line[1] == "E"
					sense = equal
				elseif split_line[1] == "G"
					sense = greater_or_equal
				elseif split_line[1] == "L"
					sense = less_or_equal
				else
					error("unknown sense $(split_line[1])!")
				end

				if row_status == "FIRST"
					push!(senses_1, sense)
					push!(first_stage_row_names, split_line[2])
					push!(first_stage_rows, row_info(sense, split_line[2]))
				elseif row_status == "SECOND"
					push!(senses_2, sense)
					push!(second_stage_row_names, split_line[2])
					push!(second_stage_rows, row_info(sense, split_line[2]))
				else
					error("Reached Here without row status! with line $line")
				end

				continue
			end

			if status == "COLUMNS"
				if startswith(line, "RHS")
					status = "RHS"
					if length(second_stage_col_names) == 0 &&
					   length(first_stage_col_names) > 0
						A_1 = hcat(A_1, current_column_first_stage_rows)
						A_2 = hcat(A_2, current_column_second_stage_rows)
					elseif length(second_stage_col_names) > 0
						A_3 = hcat(A_3, current_column_second_stage_rows)
					end

					if length(c_1) != length(first_stage_col_names)
						c_1 = vcat(c_1, zeros(length(first_stage_col_names) - length(c_1)))
					end
					if length(c_2) != length(second_stage_col_names)
						c_2 = vcat(c_2, zeros(length(second_stage_col_names) - length(c_2)))
					end
					continue
				end

				col_name, row_name, coefficient = split(line)[1:3]
				coefficient = parse(Float64, coefficient)

				if col_name == first_stage_col
					column_status = "FIRST"
				elseif col_name == second_stage_col
					column_status = "SECOND"
				end

				if column_status == "FIRST"
					if !(col_name in first_stage_col_names)
						if length(first_stage_col_names) > 0
							A_1 = hcat(A_1, current_column_first_stage_rows)
							A_2 = hcat(A_2, current_column_second_stage_rows)
							current_column_first_stage_rows =
								zeros(Float64, length(first_stage_rows))
							current_column_second_stage_rows =
								zeros(Float64, length(second_stage_rows))
						end
						push!(first_stage_col_names, col_name)
					end
				end

				if column_status == "SECOND"
					if !(col_name in second_stage_col_names)
						if length(second_stage_col_names) == 0 &&
						   length(first_stage_col_names) > 0
							A_1 = hcat(A_1, current_column_first_stage_rows)
							A_2 = hcat(A_2, current_column_second_stage_rows)
						elseif length(second_stage_col_names) > 0
							A_3 = hcat(A_3, current_column_second_stage_rows)
						end
						current_column_first_stage_rows =
							zeros(Float64, length(first_stage_rows))
						current_column_second_stage_rows =
							zeros(Float64, length(second_stage_rows))
						push!(second_stage_col_names, col_name)
					end
				end

				if row_name == objective_name
					if column_status == "FIRST"
						if length(c_1) != length(first_stage_col_names) - 1
							c_1 = vcat(
								c_1,
								zeros(length(first_stage_col_names) - 1 - length(c_1))
							)
						end
						push!(c_1, coefficient)
					elseif column_status == "SECOND"
						if length(c_2) != length(second_stage_col_names) - 1
							c_2 = vcat(
								c_2,
								zeros(length(second_stage_col_names) - 1 - length(c_2))
							)
						end
						push!(c_2, coefficient)
					else
						error("Reached here without column status!")
					end
					continue
				end

				if row_name in first_stage_row_names
					if column_status == "FIRST"
						for i in eachindex(first_stage_rows)
							if first_stage_rows[i].name == row_name
								current_column_first_stage_rows[i] = coefficient
								break
							end
						end
					elseif column_status == "SECOND"
						error("Second stage column cannot have first stage row!")
					else
						error("Reached here without column status!")
					end
					continue
				end

				if row_name in second_stage_row_names
					for i in eachindex(second_stage_rows)
						if second_stage_rows[i].name == row_name
							current_column_second_stage_rows[i] = coefficient
							break
						end
					end
				end

				continue
			end

			if status == "RHS"
				if startswith(line, "BOUNDS")
					status = "BOUNDS"
					continue
				elseif startswith(line, "ENDATA")
					break
				end

				row_name, coefficient = split(line)[2:3]
				coefficient = parse(Float64, coefficient)

				if row_name in first_stage_row_names
					for i in eachindex(first_stage_rows)
						if first_stage_rows[i].name == row_name
							b_1[i] = coefficient
							break
						end
					end
				elseif row_name in second_stage_row_names
					for i in eachindex(second_stage_rows)
						if second_stage_rows[i].name == row_name
							b_2[i] = coefficient
							break
						end
					end
				else
					error("invalid row name $(row_name)!")
				end
				continue
			end

			if status == "BOUNDS"
				if startswith(line, "ENDATA")
					break
				end

				# Only works for lower bounds (LO BND)! #TODO

				split_line = split(line)
				if !(split_line[1] == "LO" && split_line[2] == "BND")
					error("Bounds other than lower bound not supported!")
				end

				col_name, coefficient = split_line[3:4]
				coefficient = parse(Float64, coefficient)

				processed::Bool = false
				for i in eachindex(first_stage_col_names)
					if first_stage_col_names[i] == col_name
						if d_l === nothing
							d_l = zeros(length(first_stage_col_names))
						end
						d_l[i] = coefficient
						processed = true
						break
					end
				end

				if !processed
					for i in eachindex(second_stage_col_names)
						if second_stage_col_names[i] == col_name
							if d1_l === nothing
								d1_l = zeros(length(second_stage_col_names))
							end
							d1_l[i] = coefficient
							processed = true
							break
						end
					end
				end

				if !processed
					error("invalid column name $(col_name)!")
				end

				continue
			end
		end
	end

	scenario_elements::Vector{Vector{scenario_elem}} =
		[Vector{scenario_elem}() for _ in eachindex(second_stage_rows)]

	open(sto_file_name, "r") do io
		status::String = ""
		current_row_index::Integer = 0
		current_row_name::String = ""
		for line in eachline(io)

			if status == ""
				if startswith(line, "STOCH")
					status = "STOCH"
				end
				continue
			end

			if status == "STOCH"
				if startswith(line, "INDEP")
					status = "INDEP"
				end
				continue
			end

			if status == "INDEP"
				if startswith(line, "ENDATA")
					# check a probability of one
					@assert(
						abs(1 - sum(s -> s.prob, scenario_elements[current_row_index])) <
						1e-3
					)
					break
				elseif startswith(line, "*")
					continue
				end

				split_line = split(line)
				row_name, value = split_line[2:3]
				prob = last(split_line)

				value = parse(Float64, value)
				prob = parse(Float64, prob)

				if row_name != current_row_name && current_row_index > 0
					# check a probability of one
					@assert(
						abs(1 - sum(s -> s.prob, scenario_elements[current_row_index])) <
						1e-3
					)
				end

				if row_name != current_row_name
					processed::Bool = false
					for i in eachindex(second_stage_rows)
						if second_stage_rows[i].name == row_name
							current_row_index = i
							processed = true
							break
						end
					end
					if !processed
						error("Could not find row $(row_name)!")
					end
					current_row_name = row_name
				end

				push!(scenario_elements[current_row_index], scenario_elem(value, prob))
			end
		end
	end

	# Each row in the second stage must have a dimension of uncertainty. 
	# If none is found, add a trivial scenario of 0 w.p. 1
	for i in eachindex(scenario_elements)
		if length(scenario_elements[i]) == 0
			scenario_elements[i] = [scenario_elem(0.0, 1.0)]
		end
	end

	return (
		sp_instance(
			name,
			length(first_stage_col_names),
			length(second_stage_col_names),
			length(first_stage_row_names),
			length(second_stage_row_names),
			c_1, c_2,
			A_1, A_2, A_3,
			senses_1, senses_2,
			b_1, b_2,
			d_l, d_u, d1_l, d1_u
		),
		scenario_elements
	)
end


function initialize_model_with_scenarios(
	sp::sp_instance,
	scenarios::Vector{scenario_vec},
	risk_measure_alpha_values::Vector{Float64},
	time_limit_sec::Float64
)
	# Need correct alpha values
	@assert(all(alpha -> 0 <= alpha < 1, risk_measure_alpha_values))

	# Model:
	model = Model(HiGHS.Optimizer)
	set_silent(model)
	set_time_limit_sec(model, time_limit_sec)

	indices_cols1 = 1:sp.first_stage_cols
	indices_cols2 = 1:sp.second_stage_cols
	indices_rows1 = 1:sp.first_stage_rows
	indices_rows2 = 1:sp.second_stage_rows
	num_of_cvars = length(risk_measure_alpha_values)
	indices_risk_measure = 1:num_of_cvars
	indices_scenarios = 1:length(scenarios)


	# Variables
	nothing
	if sp.d_l === nothing && sp.d_u === nothing
		model[:x] = @variable(model, x[i in indices_cols1])
	elseif sp.d_l === nothing
		model[:x] = @variable(model, x[i in indices_cols1] <= sp.d_u[i])
	elseif sp.d_u === nothing
		model[:x] = @variable(model, x[i in indices_cols1] >= sp.d_l[i])
	else
		model[:x] = @variable(model, x[i in indices_cols1] in [sp.d_l[i], sp.d_u[i]])
	end

	if sp.d1_l === nothing && sp.d1_u === nothing
		@variable(model, y[i in indices_cols2, indices_scenarios])
	elseif sp.d1_l === nothing
		@variable(model, y[i in indices_cols2, indices_scenarios] <= sp.d1_u[i])
	elseif sp.d1_u === nothing
		@variable(model, y[i in indices_cols2, indices_scenarios] >= sp.d1_l[i])
	else
		@variable(
			model,
			y[i in indices_cols2, indices_scenarios] in [sp.d1_l[i], sp.d1_u[i]]
		)
	end

	# CVaR variables (u_i ∈ ℝ and infinite z_i(ξ) ≥ 0)
	@variable(model, u[indices_risk_measure])
	@variable(model, 0 <= z[indices_risk_measure, indices_scenarios])
	# outcome c_1'x + c_2'y variable (for ease of use)
	@variable(model, f[indices_scenarios])
	# Objective
	@objective(
		model,
		Min,
		sum(
			(
				u[i] +
				sum(
					z[i, s] * scenarios[s].prob for s in indices_scenarios
				)/(1 - risk_measure_alpha_values[i])
			) for i in indices_risk_measure
		)/num_of_cvars)

	# Define outcome variable
	@constraint(
		model,
		[s in indices_scenarios],
		f[s] == dot(sp.c_1, x) + dot(sp.c_2, y[:, s])
	)
	# Cvar Constraints (for the outcomes)
	@constraint(
		model,
		[i in indices_risk_measure, s in indices_scenarios],
		z[i, s] + u[i] >= f[s]
	)

	# Constraints (switching for sense)
	for j in indices_rows1
		if sp.senses_1[j] == equal
			@constraint(
				model,
				dot(sp.A_1[j, :], x) == sp.b_1[j]
			)
		elseif sp.senses_1[j] == greater_or_equal
			@constraint(
				model,
				dot(sp.A_1[j, :], x) >= sp.b_1[j]
			)
		elseif sp.senses_1[j] == less_or_equal
			@constraint(
				model,
				dot(sp.A_1[j, :], x) <= sp.b_1[j]
			)
		end
	end

	for j in indices_rows2
		if sp.senses_2[j] == equal
			for s in indices_scenarios
				@constraint(
					model,
					dot(sp.A_2[j, :], x) + dot(sp.A_3[j, :], y[:, s]) ==
					sp.b_2[j] + scenarios[s].vec[j]
				)
			end
		elseif sp.senses_2[j] == greater_or_equal
			for s in indices_scenarios
				@constraint(
					model,
					dot(sp.A_2[j, :], x) + dot(sp.A_3[j, :], y[:, s]) >=
					sp.b_2[j] + scenarios[s].vec[j]
				)
			end
		elseif sp.senses_2[j] == less_or_equal
			for s in indices_scenarios
				@constraint(
					model,
					dot(sp.A_2[j, :], x) + dot(sp.A_3[j, :], y[:, s]) <=
					sp.b_2[j] + scenarios[s].vec[j]
				)
			end
		end
	end
	return model
end

function solve_model(model::Model)::Tuple{Vector{Float64}, Float64, Bool}
	optimize!(model)
	model_solved = is_solved_and_feasible(model)
	x_sol = []
	obj = 0.0
	if model_solved
		x_sol = value.(model[:x])
		obj = objective_value(model)
	end
	return (x_sol, obj, model_solved)
end

# Returns positive infinity if not solved
function solve_model_for_objective(model::Model)::Float64
	optimize!(model)
	if is_solved_and_feasible(model)
		return objective_value(model)
	else
		return Inf
	end
end

function solve_model_with_scenarios_safe(
	sp::sp_instance,
	scenarios::Vector{scenario_vec},
	risk_measure_alpha_values::Vector{Float64},
	time_limit_sec::Float64
)::Tuple{Vector{Float64}, Float64} # Returns solution and objective
	sol, obj, model_solved = solve_model(
		initialize_model_with_scenarios(
			sp,
			scenarios,
			risk_measure_alpha_values,
			time_limit_sec
		)
	)
	if !model_solved
		error("Model was not feasible!")
	end
	return (sol, obj)
end

# function initialize_model_with_scenario(
# 	sp::sp_instance,
# 	risk_measure_alpha_values::Vector{Float64},
# 	time_limit_sec::Float64,
# 	fixed_xi::Vector{Float64}
# )::Model
# 	return initialize_model_with_scenarios(
# 		sp,
# 		[scenario_vec(fixed_xi, 1.0)],
# 		risk_measure_alpha_values,
# 		time_limit_sec
# 	)
# end

function compute_outcome_values_for_each_scenario_and_fixed_x(
	sp::sp_instance,
	time_limit_sec::Float64,
	scenarios::Vector{scenario_vec},
	fixed_x::Vector{Float64}
)::Vector{Float64}
	# Set the f value to +∞ if f(ξ, x) is infeasible
	return [
		solve_model_for_objective(
			initialize_model_with_scenario_and_fixed_x(
				sp,
				scenario.vec,
				fixed_x,
				time_limit_sec
			)
		) for scenario in scenarios
	]
end

# function compute_outcome_values_for_each_scenario_and_fixed_x(
# 	sp::sp_instance,
# 	time_limit_sec::Float64,
# 	scenarios::Vector{scenario_vec},
# 	fixed_x::Vector{Float64}
# )::Union{Vector{Float64}, Nothing}
# 	model = initialize_model_with_scenarios(sp, scenarios, [0.0], time_limit_sec)

# 	xref = model[:x]
# 	lb = Vector{Float64}(undef, length(xref))
# 	has_lower_bounds = false
# 	ub = Vector{Float64}(undef, length(xref))
# 	has_upper_bounds = false
# 	for i in eachindex(fixed_x)
# 		if has_lower_bound(xref[i])
# 			lb[i] = lower_bound(xref[i])
# 			has_lower_bounds = true
# 		end
# 		if has_upper_bound(xref[i])
# 			ub[i] = upper_bound(xref[i])
# 			has_upper_bounds = true
# 		end
# 		fix(xref[i], fixed_x[i]; force = (has_lower_bounds || has_upper_bounds))
# 	end

# 	optimize!(model)
# 	if is_solved_and_feasible(model)
# 		return value.(model[:f])
# 	else
# 		return nothing
# 	end
# end

function evaluate_sp_for_fixed_x(
	model::Model,
	fixed_x::Vector{Float64}
)::Tuple{Model, Float64, Bool}

	xref = model[:x]

	@assert(length(xref) == length(fixed_x))

	lb = Vector{Float64}(undef, length(xref))
	has_lower_bounds = false
	ub = Vector{Float64}(undef, length(xref))
	has_upper_bounds = false
	for i in eachindex(fixed_x)
		if has_lower_bound(xref[i])
			lb[i] = lower_bound(xref[i])
			has_lower_bounds = true
		end
		if has_upper_bound(xref[i])
			ub[i] = upper_bound(xref[i])
			has_upper_bounds = true
		end
		fix(xref[i], fixed_x[i]; force = (has_lower_bounds || has_upper_bounds))
	end

	_, objective, model_solved = solve_model(model)


	for i in eachindex(xref)
		unfix(xref[i])
		if has_lower_bounds
			if lb[i] === UndefInitializer()
				error()
			end
			set_lower_bound(xref[i], lb[i])
		end
		if has_upper_bounds
			if ub[i] === UndefInitializer()
				error()
			end
			set_upper_bound(xref[i], ub[i])
		end
	end

	# GRBupdatemodel(model)

	return (model, objective, model_solved)
end

function initialize_model_with_scenario_and_fixed_x(
	sp::sp_instance,
	scenario::Vector{Float64},
	x::Vector{Float64},
	time_limit_sec::Float64
)
	# Model:
	model = Model(HiGHS.Optimizer)
	set_silent(model)
	set_time_limit_sec(model, time_limit_sec)

	indices_cols2 = 1:sp.second_stage_cols
	indices_rows2 = 1:sp.second_stage_rows

	if sp.d1_l === nothing && sp.d1_u === nothing
		@variable(model, y[i in indices_cols2])
	elseif sp.d1_l === nothing
		@variable(model, y[i in indices_cols2] <= sp.d1_u[i])
	elseif sp.d1_u === nothing
		@variable(model, y[i in indices_cols2] >= sp.d1_l[i])
	else
		@variable(
			model,
			y[i in indices_cols2] in [sp.d1_l[i], sp.d1_u[i]]
		)
	end

	# Objective
	@objective(model, Min, dot(sp.c_1, x) + dot(sp.c_2, y))

	for j in indices_rows2
		if sp.senses_2[j] == equal
			@constraint(
				model,
				dot(sp.A_2[j, :], x) + dot(sp.A_3[j, :], y) ==
				sp.b_2[j] + scenario[j]
			)
		elseif sp.senses_2[j] == greater_or_equal
			@constraint(
				model,
				dot(sp.A_2[j, :], x) + dot(sp.A_3[j, :], y) >=
				sp.b_2[j] + scenario[j]
			)
		elseif sp.senses_2[j] == less_or_equal
			@constraint(
				model,
				dot(sp.A_2[j, :], x) + dot(sp.A_3[j, :], y) <=
				sp.b_2[j] + scenario[j]
			)
		end
	end
	return model
end

function solve_model_with_scenarios_and_fixed_x(
	sp::sp_instance,
	scenarios::Vector{scenario_vec},
	risk_measure_alpha_values::Vector{Float64},
	time_limit_sec::Float64,
	fixed_x::Vector{Float64}
)::Tuple{Union{Nothing, Float64}, Bool}

	# Need correct alpha values
	@assert(all(alpha -> 0 <= alpha < 1, risk_measure_alpha_values))

	# Model:
	model = Model(HiGHS.Optimizer)
	set_silent(model)
	set_time_limit_sec(model, time_limit_sec)

	indices_cols1 = 1:sp.first_stage_cols
	indices_cols2 = 1:sp.second_stage_cols
	indices_rows2 = 1:sp.second_stage_rows
	num_of_cvars = length(risk_measure_alpha_values)
	indices_risk_measure = 1:num_of_cvars
	indices_scenarios = 1:length(scenarios)

	if sp.d1_l === nothing && sp.d1_u === nothing
		@variable(model, y[i in indices_cols2, indices_scenarios])
	elseif sp.d1_l === nothing
		@variable(model, y[i in indices_cols2, indices_scenarios] <= sp.d1_u[i])
	elseif sp.d1_u === nothing
		@variable(model, y[i in indices_cols2, indices_scenarios] >= sp.d1_l[i])
	else
		@variable(
			model,
			y[i in indices_cols2, indices_scenarios] in [sp.d1_l[i], sp.d1_u[i]]
		)
	end

	# CVaR variables (u_i ∈ ℝ and infinite z_i(ξ) ≥ 0)
	@variable(model, u[indices_risk_measure])
	@variable(model, 0 <= z[indices_risk_measure, indices_scenarios])
	# outcome c_1'x + c_2'y variable (for ease of use)
	@variable(model, f[indices_scenarios])
	# Objective
	@objective(
		model,
		Min,
		sum(
			(
				u[i] +
				sum(
					z[i, s] * scenarios[s].prob for s in indices_scenarios
				)/(1 - risk_measure_alpha_values[i])
			) for i in indices_risk_measure
		)/num_of_cvars)

	# Define outcome variable
	@constraint(
		model,
		[s in indices_scenarios],
		f[s] == dot(sp.c_1, fixed_x) + dot(sp.c_2, y[:, s])
	)
	# Cvar Constraints (for the outcomes)
	@constraint(
		model,
		[i in indices_risk_measure, s in indices_scenarios],
		z[i, s] + u[i] >= f[s]
	)

	for j in indices_rows2
		if sp.senses_2[j] == equal
			for s in indices_scenarios
				@constraint(
					model,
					dot(sp.A_2[j, :], fixed_x) + dot(sp.A_3[j, :], y[:, s]) ==
					sp.b_2[j] + scenarios[s].vec[j]
				)
			end
		elseif sp.senses_2[j] == greater_or_equal
			for s in indices_scenarios
				@constraint(
					model,
					dot(sp.A_2[j, :], fixed_x) + dot(sp.A_3[j, :], y[:, s]) >=
					sp.b_2[j] + scenarios[s].vec[j]
				)
			end
		elseif sp.senses_2[j] == less_or_equal
			for s in indices_scenarios
				@constraint(
					model,
					dot(sp.A_2[j, :], fixed_x) + dot(sp.A_3[j, :], y[:, s]) <=
					sp.b_2[j] + scenarios[s].vec[j]
				)
			end
		end
	end

	optimize!(model)
	model_solved = is_solved_and_feasible(model)

	return (model_solved ? objective_value(model) : Nothing, model_solved)
end