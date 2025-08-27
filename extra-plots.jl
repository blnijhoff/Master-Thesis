# Create the plot for a provided csv file
function make_box_plot_from_csv(csv_file::String)
	@assert(endswith(csv_file, ".csv"))


	file_name = splitext(basename(csv_file))[1]
	@assert(startswith(file_name, "G"))

	is_G = split(file_name, "_"; limit = 2)[1] == "G"

	values::Vector{Vector{Float64}} = []
	headers::Vector{String} = []
	# Open file with the one line `vals = parse.(Float64, readlines(file_name))`
	open(csv_file, "r") do io
		lines = eachline(io)
		# Read first line separately
		first_line = first(lines)
		headers = split(first_line, ",")
		number_of_cols = length(headers)
		values = [Vector{Float64}() for _ in 1:number_of_cols]

		for line in lines
			line = split(line, ",")
			@assert(length(line) == number_of_cols)

			for i in 1:number_of_cols
				try
					value = parse(Float64, line[i])
					push!(values[i], value)
				catch e
					if !isa(e, ArgumentError)
						rethrow(e)
					end
				end
			end
		end
	end

	make_box_plots(
		values,
		headers,
		"Scenario Generation Method",
		is_G ? L"G" : L"G_{\textrm{SP}}",
		"Boxplots of " *
		(is_G ? string(L"G\;") : string(L"G_{\mathrm{SP}}\;")) *
		"for different scenario generation methods",
		save_file = true,
		file_name = "extra-output/$file_name.png"
	)
end
