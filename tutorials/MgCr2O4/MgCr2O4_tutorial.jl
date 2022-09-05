# === Loading Packages ===
using Sunny # The main package
using ProgressMeter, Statistics, Formatting, StaticArrays, LinearAlgebra, DelimitedFiles, OffsetArrays, Random # Some useful peripheral packages
using GLMakie, Plots # Some useful plotting packages
Sunny.offline_viewers(); 

# === Define the crystal structure of a B-site pyrochlore lattice with only B-site atoms === 
lat_vecs = lattice_vectors(8.3342, 8.3342, 8.3342, 90.0, 90.0, 90.0);  
bas_vecs = [  [0.87500000, 0.62500000, 0.37500000],
              [0.62500000, 0.12500000, 0.62500000],
              [0.87500000, 0.87500000, 0.12500000],  
              [0.62500000, 0.87500000, 0.37500000],  
              [0.87500000, 0.12500000, 0.87500000],  
              [0.62500000, 0.62500000, 0.12500000],  
              [0.87500000, 0.37500000, 0.62500000],  
              [0.62500000, 0.37500000, 0.87500000],  
              [0.37500000, 0.62500000, 0.87500000],  
              [0.12500000, 0.12500000, 0.12500000],  
              [0.37500000, 0.87500000, 0.62500000],  
              [0.12500000, 0.87500000, 0.87500000],  
              [0.37500000, 0.12500000, 0.37500000],  
              [0.12500000, 0.62500000, 0.62500000],  
              [0.37500000, 0.37500000, 0.12500000],  
              [0.12500000, 0.37500000, 0.37500000]];
bas_typs = ["B","B","B","B", "B","B","B","B", "B","B","B","B", "B","B","B","B"]; 
xtal_pyro   = Crystal(lat_vecs, bas_vecs; types=bas_typs); 

# === Return crystalographic information gathered by Sunny about the crystal structure === 
display(xtal_pyro); 
view_crystal(xtal_pyro, 3.2);
print_bond_table(xtal_pyro, 5.9); 
print_allowed_anisotropy(xtal_pyro, 1);

# === Define the crystal structure of MgCr2O4 by copying the info from a .cif file INCLUDING space group and setting === 
lat_vecs     = lattice_vectors(8.3342, 8.3342, 8.3342, 90.0, 90.0, 90.0);  
bas_vecs     = [ [0.12500, 0.12500, 0.12500],
                 [0.50000, 0.50000, 0.50000],
                [0.26070, 0.26070, 0.26070]]; 
bas_typs     = ["Mg","Cr","O"];
lat_spg      = 227; 
lat_set      = "2"; 
xtal_mgcro_2 = Crystal(lat_vecs, bas_vecs, lat_spg; types=bas_typs, setting=lat_set); 
xtal_mgcro   = subcrystal(xtal_mgcro_2,"Cr");

# === Return crystalographic information gathered by Sunny about the crystal structure === 
display(xtal_mgcro); 
view_crystal(xtal_mgcro,5.9)
print_bond_table(xtal_mgcro, 5.9); 
print_allowed_anisotropy(xtal_mgcro, 1);

# === Define Local Hilbert Space for Cr ===
lhs_Cr  = [SiteInfo(1; N=0, g=2.0, spin_rescaling=3/2)]; 
lhs_B   = [SiteInfo(1; N=0, g=2.0, spin_rescaling=3/2)]; 

# === Define Values of Exchange Interactions ===
val_J1      = 3.27; # value of J1 in meV from Bai's PRL paper
val_J_pyro  = [1.00,0.000,0.000,0.000]*val_J1; # pure nearest neighbor pyrochlore
val_J_mgcro = [1.00,0.0815,0.1050,0.085]*val_J1; # further neighbor pyrochlore relevant for MgCr2O4
#val_J_mgcro = [1.00,0.000,0.025,0.025]*val_J1; # this is a funny setting!

# === Build list of exchange interactions for our system ===
exint_pyro  = [heisenberg(val_J_pyro[1],  Bond(1, 3, [0, 0, 0]),"J1"),
               heisenberg(val_J_pyro[2],  Bond(1, 2, [0, 0, 0]),"J2"), 
               heisenberg(val_J_pyro[3],  Bond(2, 6, [0, 0, 0]),"J3a"), 
               heisenberg(val_J_pyro[4],  Bond(1, 5, [0, 0, 0]),"J3b")];
exint_mgcro = [heisenberg(val_J_mgcro[1], Bond(1, 2, [0, 0, 0]),"J1"),
               heisenberg(val_J_mgcro[2], Bond(1, 7, [0, 0, 0]),"J2"),
               heisenberg(val_J_mgcro[3], Bond(1, 3, [0, 0, 0]),"J3a"),   # Watch out for which is J3a see above plots
               heisenberg(val_J_mgcro[4], Bond(1, 3, [1, 0, 0]),"J3b")];   # Watch out for which is J3b see above plots           

## === Define Super Cell Size ===
scd  = (10,10,10); #Super Cell Dimension Let's go small for now

## === Construct the Spin Systems ===
sys_pyro  = SpinSystem(xtal_pyro,  exint_pyro,  scd, lhs_B);
sys_mgcro = SpinSystem(xtal_mgcro, exint_mgcro, scd, lhs_Cr);

## === Randomize the Spins ===
rand!(sys_pyro);
rand!(sys_mgcro);

## === Construct Langevin Sampler ===
nLA = 10;  # Number of Langevin time steps performed each time the Sampler is invoked
α   = 0.2; # Langevin damping, usually 0.05 or 0.1 is good.
Δt  = 0.01; # Time steps in Langevin
kT  = val_J1*20; # Initializing spin system at some finite temperature corresponding to 10 times J1 (to be well paramagnetic)
sam_LA_pyro  = LangevinSampler(sys_pyro, kT, α, Δt, nLA);
sam_LA_mgcro = LangevinSampler(sys_mgcro, kT, α, Δt, nLA);

## === Optional: Construct Metropolis Sampler ===
nMC = 5;  # Number of Langevin time steps performed each time the Sampler is invoked
kT  = val_J1*20; # Initializing spin system at some finite temperature corresponding to 10 times J1 (to be well paramagnetic)
sam_MC_pyro  = MetropolisSampler(sys_pyro, kT, nMC);
sam_MC_mgcro = MetropolisSampler(sys_mgcro, kT, nMC);

## === Thermalize System to the temperature said below using Langevin===
kT     = 1.8; # Target temperature in meV
nTherm = 1000; # Number of times the Sampler will run, here nLA*nTherm
@time begin
    set_temp!(sam_LA_pyro,kT); 
    set_temp!(sam_LA_mgcro,kT); 
    set_temp!(sam_MC_pyro,kT); 
    set_temp!(sam_MC_mgcro,kT); 
    prog = Progress(Int64(round(nTherm/100)); dt=0.10, desc="Thermalizing: ", color=:blue)
    for j in 1:nTherm/100
        thermalize!(sam_LA_pyro,Int64(round(nTherm/100)));
        thermalize!(sam_MC_pyro,Int64(round(nTherm/100)));
        thermalize!(sam_LA_mgcro,Int64(round(nTherm/100)));
        thermalize!(sam_MC_mgcro,Int64(round(nTherm/100)));
        next!(prog);  
    end
end

## === Plot the resulting spin system for the Pyrochlore ===
plot_spins(sys_pyro,arrowlength=0.5, linewidth=0.2, arrowsize=0.5)

## === Plot the resulting spin system for the MgCr2O4 ===
plot_spins(sys_mgcro,arrowlength=0.5, linewidth=0.2, arrowsize=0.5)

##  === Calculate SQ ===
nsam        = 10; # Number of samples that are averaged over  (usually 10 is good)
decor_ratio = 10; # Number of time the Sampler is called to decorelate samples between sampling
bz_size     = (8,8,8); # Size of the resulting extended Brillouin zone after FFT
@time begin
    sq_pyro  = StructureFactor(sys_pyro; bz_size=bz_size, dipole_factor=true)
    sq_mgcro = StructureFactor(sys_mgcro; bz_size=bz_size, dipole_factor=true)
    prog     = Progress(nsam ; dt=1.00, desc="Sampling: ", color=:green)
    for j in 1:nsam
        thermalize!(sam_LA_pyro,decor_ratio); 
        thermalize!(sam_LA_mgcro,decor_ratio); 
        Sunny.update!(sq_pyro,sys_pyro)
        Sunny.update!(sq_mgcro,sys_mgcro);
        next!(prog);        
    end
end

#  === Create a function to plot S(Q) ===
function PlotSQ(input_sq, input_sys; Slice, Imax)
    Q1   = range(input_sq.sfactor.offsets[1],length=size(input_sq.sfactor)[1])/size(input_sys.lattice)[1];
    Q2   = range(input_sq.sfactor.offsets[2],length=size(input_sq.sfactor)[2])/size(input_sys.lattice)[2];
    Q3   = range(input_sq.sfactor.offsets[3],length=size(input_sq.sfactor)[3])/size(input_sys.lattice)[3];
    midQ = Int64(size(Q1)[1]/2);
    Int = input_sq.sfactor[:,:,Slice,0]/prod(size(input_sys.lattice));
    return display(Plots.heatmap(Q1,Q2,Int,clim=(0,Imax)));   
end

#  === Plot the Results ===
PlotSQ(sq_pyro, sys_pyro; Slice=0, Imax=100)
PlotSQ(sq_mgcro, sys_mgcro; Slice=0, Imax=100)


#  === Calculate SQW ===
Nω   = 100;  # Number of Frequencies Calculated
ωmax = val_J1*15;
Δt   = 0.001;
nsam = 2;  # Number of samples that are averaged over  (usually 10 is good)
nLa  = 20;
bz_size = (8,8,8);
@time begin
sqw_pyro = dynamic_structure_factor(sys_pyro, sam_LA_pyro; 
    dynΔt=Δt, meas_rate=convert(Int64,round(2*pi/(Δt*ωmax))), dyn_meas=Nω,
    bz_size=bz_size, 
    thermalize=nLa,nsamples=nsam,  
    verbose=true, reduce_basis=true, dipole_factor=true,
)
end
@time begin
sqw_mgcro = dynamic_structure_factor(sys_mgcro, sam_LA_mgcro; 
    dynΔt=Δt, meas_rate=convert(Int64,round(2*pi/(Δt*ωmax))), dyn_meas=Nω,
    bz_size=bz_size, 
    thermalize=nLa,nsamples=nsam,  
    verbose=true, reduce_basis=true, dipole_factor=true,
)
end

#  === Create a function to plot S(Q,W) ===
function PlotSQW(input_sqw, input_sys; Slice1, Slice2, Imax)
    Q1   = range(input_sqw.sfactor.offsets[1],length=size(input_sqw.sfactor)[1])/size(input_sys.lattice)[1];
    Q2   = range(input_sqw.sfactor.offsets[2],length=size(input_sqw.sfactor)[2])/size(input_sys.lattice)[2];
    Q3   = range(input_sqw.sfactor.offsets[3],length=size(input_sqw.sfactor)[3])/size(input_sys.lattice)[3];
    EN   = range(input_sqw.sfactor.offsets[4],length=size(input_sqw.sfactor)[4]);
    midQ = Int64(size(Q1)[1]/2);
    midE = Int64(size(EN)[1]/2);
    Int  = input_sqw.sfactor[:,:,:,:]/prod(size(input_sys.lattice));
    return display(Plots.heatmap(Q1,EN[1:midE],Int[:,Slice1,Slice2,1:midE]',clim=(0,Imax)));   
end


#  === Plot the Results ===
PlotSQW(sqw_pyro, sys_pyro; Slice1=-5, Slice2=0, Imax=20);
PlotSQW(sqw_mgcro, sys_mgcro; Slice1=0, Slice2=0, Imax=200);
