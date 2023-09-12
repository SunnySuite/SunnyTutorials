# > ![](https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.jpg)
# _This is a [tutorial](https://github.com/SunnySuite/SunnyTutorials/tree/main/tutorials)
#  for the [Sunny](https://github.com/SunnySuite/Sunny.jl/) package,
#  which enables dynamical simulations of ordered and thermally disordered spins with dipole
#  and higher order moments._

## Welcome to a Sunny Tutorial on the $\sqrt{3}\times\sqrt{3}$ Kagome Antiferromagnet
# **Script**: Distorted Kagome Lattice Calculation <br>
# **Inspired by**:$\sqrt{3}\times\sqrt{3}$ Kagome Antiferromagnet SpinW tutorial 
# (Bjorn Fak and Sandor Toth https://spinw.org/tutorials/08tutorial).
# **Authors**: Harry Lane <br>
# **Date**: September 11, 2023  (Sunny 0.5.4) <br>
# **Goal**: This script is to calculate the linear spin wave theory spectrum for the 
# $\sqrt{3} \times \sqrt{3}$ Kagome antiferromagnet and compare with the results from SpinW.  

# ---
# #### Loading Packages 

using Sunny, GLMakie

# #### Defining Custom Functions For This Script

function anneal!(sys,  integrator;  kTschedule, ndwell)
    nspins = prod(size(sys.dipoles));
    ensys  = zeros(length(kTschedule))
    for (i, kT) in enumerate(kTschedule)
        integrator.kT = kT
        for _ in 1:ndwell
            step!(sys, integrator)
        end
        ensys[i] = energy(sys)
    end
    return ensys/nspins
end


# ---
# ### System Definition 
# Set up a [`Crystal`](@ref) with $P\overline{3}$ space group and Cr$^{+}$ ions on each site.

a = b = 6.0 # (Å)
c = 40.0
latvecs = lattice_vectors(a, b, c, 90, 90, 120)
types = ["Cr"]
crystal = Crystal(latvecs, [[1/2,0,0]], 147;types)

# The next step is to add interactions. The command [`print_symmetry_table`](@ref) shows all 
# symmetry-allowed interactions up to a cutoff distance.
print_symmetry_table(crystal,7.0)

# We can now create a [`System`](@ref) with $S=1$ spins and specify [`dipole`](@ref) mode.
sys = System(crystal,(9,9,1),[SpinInfo(1,S=1,g=2)],:dipole)

# We now add nearest-neighbor antiferromagnetic interactions. 
set_exchange!(sys,1.0,Bond(2,3,[0,0,0]))

# ---
# ### Search for the system ground-state

# We now have two options. If we happen to know that the ground state we 
# can use [`define_magnetic_structure`](@ref) to prescribe a magnetic sturcture, or we can 
# take advantage of the tools that Sunny provides to find the magnetic ground state. 
# Let us try simulated annealing followed by gradient descent. 

randomize_spins!(sys)
Δt = 0.01
λ = 0.1
langevin = Langevin(Δt; kT=0.0, λ)
kTs = [100 * 0.9^k for k in 0:200];
anneal!(sys,langevin;kTschedule=kTs,ndwell=10_000);
langevin.kT = 0.0
for _ ∈ 1:10_000
    step!(sys, langevin)
end
for i ∈ 1:20
    minimize_energy!(sys;maxiters=10_000)
end

# Given the frustrated nature of the interactions, it is difficult to thermalize a single 
# $\mathbf{Q}$ magnetic structure. 

print_wrapped_intensities(sys)
plot_spins(sys)

#Instead let's build the known magnetic structure.
# Let's see if we can do better than the (3,3,1) lattice that we would naively choose.

kvec = -[1/3,1/3,0]
n    = [0,0,1]
Slist=[[0,1,0],[0,1,0],[-1,-1,0]]
for nsub in 1:3
    set_spiral_order_on_sublattice!(sys, nsub ;q=kvec,axis=n,S0=Slist[nsub])
end
plot_spins(sys)
print_wrapped_intensities(sys)

# Success! This is telling us that we can describe the structure using 
# latsize=(3,3,1).

# ---
# ### Calculating Spin Wave Theory

# The next step is to build a [`SpinWaveTheory`](@ref) object.

swt = SpinWaveTheory(sys;energy_ϵ=4e-2);

# Before calculating the linear spin wave intensity, we must define a path in reciprocal space. 
points_rlu = [[-1/2 0 0],[0 0 0],[1 1 0]];
density = 100
path, xticks = reciprocal_space_path(crystal,points_rlu, density);

# Defining the scattering formula
γ = 0.15 # width in meV
broadened_formula = intensity_formula(swt, :perp; kernel=lorentzian(γ))

# Calculating the spin-waves
energies = collect(0:0.01:6)  # 0 < ω < 6 (meV).
is = intensities_broadened(swt, [q for q in path], energies, broadened_formula)
fig = Figure()
ax = Axis(fig[1,1]; ylabel="Energy (meV)", xticks, xticklabelrotation=π/6)
heatmap!(ax, 1:size(is, 1), energies, is,colorrange = (0,2.0))
fig