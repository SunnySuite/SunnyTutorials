using Literate, Dates, Git

root = dirname(@__FILE__)
example_dir = joinpath(root, "Sunny.jl", "examples", "longer_examples")
save_dir = joinpath(root, "Tutorials")

# If no arguments given, rebuild all notebooks
println(length(ARGS))
if length(ARGS) == 1
    tutorials = filter(name -> split(name, ".")[end] == "jl", readdir(example_dir))
else
    tutorials = ARGS[2:end]
end

# Build the notebooks
map(tutorials) do tutorial
    Literate.notebook(joinpath(example_dir, tutorial), save_dir; execute=true)
end

# Sync with github
run(`$(git()) pull`) # Make sure up to date
run(git(["add", save_dir*"/*.ipynb"]))
run(`$(git()) commit -am "Auto-generated tutorials $(string(Dates.now()))"`)
run(`$(git()) push`)
