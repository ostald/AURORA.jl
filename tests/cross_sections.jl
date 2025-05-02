

using Pkg
AURORA_dir = "/mnt/data/oliver/AURORA.jl/"
Pkg.activate(AURORA_dir)
using AURORA
using CairoMakie
set_theme!(Theme(colormap = :hawaii))



## Setting parameters
altitude_lims = [60, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-10:0;            # (°) angle-limits for the electron beams
E_max = 1e5;                   # (eV) upper limit to the energy grid
E_min = 1;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith


h_atmosphere_models = altitude_lims[1]-9:1:altitude_lims[2]+200
msis_file = find_nrlmsis_file(
    year=2006, month=12, day=12, hour=19, minute=30, lat=70, lon=19, height=h_atmosphere_models
    );
iri_file = find_iri_file(
    year=2006, month=12, day=12, hour=19, minute=30, lat=70, lon=19, height=h_atmosphere_models
    );


## Define input parameters
input_type = "constant_onset"
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = altitude_lims[2];          # (km) altitude of the source
#E_min = E_max - 100;         # (eV) bottom energy of the FAB
Beams = 1:length(θ_lims)-1;   # beam numbers for the precipitation, starting with field aligned down
t0 = 0;                     # (s) time of start for smooth transition
t1 = 0;                     # (s) time of end for smooth transition


INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);

h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);


cs_N2 = σ_neutrals.σ_N2
cs_O2 = σ_neutrals.σ_O2
cs_O  = σ_neutrals.σ_O


##
fig = Figure();
ax = Axis(fig[1, 1],
            xscale=log10,
            yscale=log10, 
            xlabel="Energy [eV]",
            ylabel="Cross section [m²]");
for i in axes(cs_N2, 1)
    lines!(E, cs_N2[i, :])
end
display(fig)

##

fig = Figure();
ax = Axis(fig[1, 1],
            xscale=log10,
            yscale=log10, 
            xlabel="Energy [eV]",
            ylabel="Cross section [m²]");
for i in axes(cs_O2, 1)
    lines!(E, cs_O2[i, :])
end
display(fig)
##

fig = Figure();
ax = Axis(fig[1, 1],
            xscale=log10,
            yscale=log10, 
            xlabel="Energy [eV]",
            ylabel="Cross section [m²]");
for i in axes(cs_O, 1)
    lines!(E, cs_O[i, :])
end
display(fig)
##