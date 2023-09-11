using Sunny, LinearAlgebra, GLMakie

# Load experiment data (so we have histogram_parameters)
histogram_parameters, data = load_nxs("experiment_data_normalized.nxs")

# Crystallography & Chemistry
cryst = Crystal("example_cif.cif"; symprec=1e-4)
sys_chemical = System(subcrystal(cryst,"Cr"), (1,1,1), [SpinInfo(1,S=3/2,g=2)], :SUN)
sys = reshape_supercell(sys_chemical, [1 1 0; 1 -1 0; 0 0 1]) # Neel state

plot_spins(sys)

function get_Z(is)
  sum(is[7:14,7:14,:,4]) # Qx, Qy, and E range of magnetic bragg peak
end

# Multi-sampling magic numbers
msaa4 = [[0.625, 0.625, 0.125]
        ,[0.875, 0.125, 0.375]
        ,[0.375, 0.375, 0.875]
        ,[0.125, 0.875, 0.625]]
energy_multisample = [(n + 0.5)/5 for n = 1:5]

# Known J2
J2 = 0.16
set_exchange!(sys,J2,Bond(1,1,[1,1,0]))

# Unknown J1, D
function forward_problem(J1,D)
  set_exchange!(sys,J1,Bond(1,1,[1,0,0]))
  
  Sz = spin_operators(sys,1)[3]
  set_onsite_coupling!(sys,D*Sz^2,1)

  # Standard calculation:
  randomize_spins!(sys)
  minimize_energy!(sys)
  minimize_energy!(sys)
  minimize_energy!(sys)

  swt = SpinWaveTheory(sys)
  
  formula = intensity_formula(swt,:perp;
    # TODO: instrument-adapted broadening
    kernel = lorentzian(2.)
    ,mode_fast = true
    ,formfactors = [FormFactor("Cr3")]
    )

  
  intensity, counts = Sunny.intensities_bin_multisample(swt
                                                       ,histogram_parameters
                                                       ,msaa4
                                                       ,energy_multisample
                                                       ,formula)
  return intensity ./ counts
end

J1_guess = 11.3
A_guess = 0.03
intensities_guess = forward_problem(J1_guess, A_guess)
loss = loss_function(data, intensities_guess)


function loss_function(experiment_data,simulation_data)
  Z_experiment = get_Z(experiment_data)
  normalized_exp_data = experiment_data ./ Z_experiment

  Z_sunny = get_Z(simulation_data)
  normalized_sim_data = simulation_data ./ Z_sunny

  weights = 1.

  # Compute squared error over every histogram bin
  squared_errors = (normalized_exp_data .- normalized_sim_data).^2
  squared_errors[isnan.(squared_errors)] .= 0 # Filter out missing experiment data
  squared_errors[:,:,:,1:2] .= 0 # Filter out elastic line
  sqrt(sum(weights .* squared_errors))
end

# Plot the data
#using GLMakie


