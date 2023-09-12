using Literate, Dates, Git

root = dirname(@__FILE__)
example_dir = joinpath(root,"..", "Sunny.jl", "examples", "longer_examples")
save_dir = joinpath(root, "..", "Tutorials")

# If no arguments given, rebuild all notebooks
tutorials = if length(ARGS) == 0
    filter(name -> split(name, ".")[end] == "jl", readdir(example_dir))
else
    ARGS
end

# Build the notebooks
map(tutorials) do tutorial
    Literate.notebook(joinpath(example_dir, tutorial), save_dir; execute=true)
end

# Sync with github
cd("..")
run(`$(git()) pull`) # Make sure up to date
run(git(["add", save_dir*"/*.ipynb"]))
run(`$(git()) commit -am "Auto-generated tutorials $(string(Dates.now()))"`)
run(`$(git()) push`)
