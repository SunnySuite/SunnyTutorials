using Sunny, GLMakie
using LinearAlgebra, Random

function FeI₂(dims; spin_rescaling=1.0, rng = nothing)
    isnothing(rng) && (rng = MersenneTwister())

    ## Crystal
    a = b = 4.05012
    c = 6.75214
    lat_vecs = lattice_vectors(a, b, c, 90, 90, 120)
    basis_vecs = [[0,0,0]]
    crystal = Crystal(lat_vecs, basis_vecs, 164; setting="1")

    ## Interactions
    J1pm   = -0.236
    J1pmpm = -0.161
    J1zpm  = -0.261
    J2pm   = 0.026
    J3pm   = 0.166
    J′0pm  = 0.037
    J′1pm  = 0.013
    J′2apm = 0.068
    D      = 2.165

    J1zz   = -0.236
    J2zz   = 0.113
    J3zz   = 0.211
    J′0zz  = -0.036
    J′1zz  = 0.051
    J′2azz = 0.073
    J′2bzz = 0.0

    J1xx = J1pm + J1pmpm 
    J1yy = J1pm - J1pmpm
    J1yz = J1zpm
    J₁ = [J1xx  0.0   0.0;
        0.0   J1yy  J1yz;
        0.0   J1yz  J1zz]
    J₂ = [J2pm  0.0  0.0;
        0.0   J2pm 0.0;
        0.0   0.0  J2zz]
    J₃ = [J3pm   0.0   0.0;
        0.0    J3pm  0.0;
        0.0    0.0   J3zz]
    J′₀ = [J′0pm  0.0   0.0;
        0.0    J′0pm 0.0;
        0.0    0.0   J′0zz]
    J′₁ = [J′1pm  0.0   0.0;
        0.0    J′1pm 0.0;
        0.0    0.0   J′1zz]
    J′₂ = [J′2apm  0.0   0.0;
        0.0    J′2apm 0.0;
        0.0    0.0   J′2azz]

    Sz = Sunny.gen_spin_ops(3)[3]
    D = -2.165
    J_ani = Sz*Sz  

    interactions = [
        exchange(J₁, Bond(1,1,[1,0,0])),
        exchange(J₂, Bond(1,1,[1,2,0])),
        exchange(J₃, Bond(1,1,[2,0,0])),
        exchange(J′₀, Bond(1,1,[0,0,1])),
        exchange(J′₁, Bond(1,1,[1,0,1])),
        exchange(J′₂, Bond(1,1,[1,2,1])),
        Sunny.SUNAnisotropy(D*J_ani, 1, "")
    ]

    # SiteInfos
    site_infos = [SiteInfo(1; N=3, g=2*I(3), spin_rescaling)]

    return SpinSystem(crystal, interactions, dims, site_infos; rng)
end

begin
    rng = MersenneTwister(112)
    sys = FeI₂((4, 4, 4); rng)
    rand!(sys)
    kT = 0.01

    live_langevin_integration(sys, 100, 0.002, 0.00;
        linewidth=0.3, arrowsize=0.7,
        arrowlength=2.5, α=0.1, framerate=60,
    )
end