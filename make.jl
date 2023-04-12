using Literate, Dates, Git

root = dirname(@__FILE__)
example_dir = joinpath(root, "Sunny.jl", "examples", "longer_examples")
save_dir = joinpath(root, "Tutorials")
tutorials = filter(name -> split(name, ".")[end] == "jl", readdir(example_dir))

map(tutorials) do tutorial
    Literate.notebook(joinpath(example_dir, tutorial), save_dir; execute=true)
end

date = string(Dates.now())
run(git(["add", save_dir*"/*.ipynb"]))
run(`$(git()) commit -am "Auto-generated tutorials $date"`)
run(`$(git()) push`)