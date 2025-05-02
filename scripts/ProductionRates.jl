using Pkg
AURORA_dir = "/mnt/data/oliver/AURORA.jl/"
Pkg.activate(AURORA_dir)
using AURORA
using MAT
using CairoMakie
set_theme!(Theme(colormap = :hawaii))


root_savedir = "ionprod_v6"  # name of the root folder
e_grid = 10 .^(range(1,stop=5,length=400))


#include("scripts/Control_script_SteadyState.jl")

## Run the analysis
directory_to_process = joinpath(AURORA_dir, "data", root_savedir)
subdirs = readdir(directory_to_process, join=true)
"""
for dir in subdirs
    make_volume_excitation_file(dir)
end
"""

dE = diff(e_grid)                                   
E = e_grid[1:end-1] + dE/2                          # already centre energy bins

h_atm = []
QN2i      = []
QO1D      = []
QO1S      = []
QO2i      = []
QOi       = []
Q_from_O  = []
Q_from_O2 = []
Q_from_N2 = []

for i in 1:length(e_grid)-1
    E_max = e_grid[i+1]
    E_min = e_grid[i]
   
    name_savedir = string(E_min)*'-'*string(E_max) * "eV"
    f = matopen(joinpath(directory_to_process, name_savedir, "Qzt_all_L.mat"))
        push!(h_atm, read(f, "h_atm"))
        push!(QN2i      , read(f, "QN2i"     ))
        push!(QO1D      , read(f, "QO1D"     ))
        push!(QO1S      , read(f, "QO1S"     ))
        push!(QO2i      , read(f, "QO2i"     ))
        push!(QOi       , read(f, "QOi"      ))
        push!(Q_from_O  , read(f, "Q_from_O" ))
        push!(Q_from_O2 , read(f, "Q_from_O2"))
        push!(Q_from_N2 , read(f, "Q_from_N2"))
    close(f)
end

h_atm = h_atm[1]
QN2i      = reduce(hcat, QN2i     )
QO1D      = reduce(hcat, QO1D     )
QO1S      = reduce(hcat, QO1S     )
QO2i      = reduce(hcat, QO2i     )
QOi       = reduce(hcat, QOi      )

Q_from_O  = cat(Q_from_O...,  dims = 3)
Q_from_O2 = cat(Q_from_O2..., dims = 3)
Q_from_N2 = cat(Q_from_N2..., dims = 3)



fig = Figure()
ax = Axis(fig[1, 1], 
            xlabel="Energy [keV]", 
            ylabel="Height [km]")
hm = heatmap!(E/1e3,
                h_atm/1e3,
                QN2i', 
                colorscale = log10, 
                colorrange = (1e5, 1e11))
cb = Colorbar(fig[1, 2], 
                hm, 
                label = "Ionization rate [m⁻³ s⁻¹]")
display(fig)


fig = Figure()
ax = Axis(fig[1, 1], 
            xlabel="Energy [keV]", 
            ylabel="Height [km]", 
            xscale=log10)
hm = heatmap!(E,
                h_atm/1e3,
                QN2i', 
                colorscale = log10, 
                colorrange = (1e5, 1e11))
cb = Colorbar(fig[1, 2], 
                hm, 
                label = "Ionization rate [m⁻³ s⁻¹]")
display(fig)




fig = Figure();
ax = Axis(fig[1, 1], 
            xscale = log10, 
            xlabel="Ionization rate [m⁻³ s⁻¹]", 
            ylabel="Height [km]")
xlims!(1e1, 1e77)
for i in 1:length(E)
    lines!(QN2i[:, i], h_atm/1e3)
end
display(fig)
save("i_profile.pdf", fig, pdf_version="1.4")


fig = Figure();
ax = Axis(fig[1, 1], 
            xscale = log10, 
            xlabel = "Ionization rate [m⁻³ s⁻¹]", 
            ylabel="Height [km]")
xlims!(1e1, 1e77)
for i in 1:10:length(E)
    lines!(QN2i[:, i], h_atm/1e3)
end
display(fig)



fig = Figure();
ax = Axis3(fig[1, 1], 
            limits=((50, 400), nothing, (2, 15)), 
            xlabel = "Height [km]", 
            ylabel="log₁₀(E) [eV]", 
            zlabel= "log₁₀(qₑ) [m⁻³ s⁻¹]")
hm = surface!(h_atm/1e3, 
                log10.(E), 
                log10.(QN2i))
display(fig)



fig = Figure();
ax = Axis3(fig[1, 1], 
            limits=((50, 600), nothing, nothing), 
            xlabel = "Height [km]", 
            ylabel="log₁₀(E) [eV]", 
            zlabel= "log₁₀(qₑ) [m⁻³ s⁻¹]"
            ) #, zscale=log10 )
hm = surface!(h_atm/1e3, 
                log10.(E), 
                log10.(QN2i.+1e-16)
                )#, colorscale=log10, colorrange=(1e5, 1e12))
display(fig)



QN2i_ = copy(QN2i)
QN2i_[log10.(QN2i_[:]) .< 2 ] .= NaN
QN2i_[log10.(QN2i_[:]) .> 15] .= NaN
fig = Figure();
ax = Axis3(fig[1, 1],
            limits=(nothing, nothing, (2, 15)), 
            azimuth=0.1*pi, 
            xlabel = "Height [km]", 
            ylabel="log₁₀(E) [eV]", 
            zlabel= "log₁₀(qₑ) [m⁻³ s⁻¹]")
hm = surface!(h_atm/1e3, 
                log10.(E), 
                log10.(QN2i_))
display(fig)


##
QN2i_ = copy(QN2i)
QN2i_[log10.(QN2i_[:]) .< 2 ] .= 1e2
QN2i_[log10.(QN2i_[:]) .> 13] .= 1e13
fig = Figure();
ax = Axis3(fig[1, 1], 
            limits=(nothing, nothing, (2, 15)), 
            azimuth=1.47*pi, 
            elevation=0.2*pi, 
            xlabel = "Height [km]", 
            ylabel="log₁₀(E) [eV]", 
            zlabel= "log₁₀(qₑ) [m⁻³ s⁻¹]")
hm = surface!(h_atm/1e3, 
                log10.(E), 
                log10.(QN2i_),
                colormap = :hawaii )
display(fig)
##

fig = Figure();
ax = Axis3(fig[1, 1], 
            xlabel = "Height [km]", 
            ylabel="E [eV]", 
            zlabel= "qₑ [m⁻³ s⁻¹]")#, limits=(nothing, nothing, (1e5, 1e12)))
hm = surface!(h_atm[1:200],
                E[1:200], 
                QN2i[1:200, 1:200].+1e-16)
display(fig)


fig = Figure();
ax = Axis(fig[1, 1], 
            xlabel="Ionization rate [m⁻³ s⁻¹]", 
            ylabel="Height [km]")#, xscale = log10)
xlims!(-1e5, 1e5)#1e77)
for i in 1:length(E)
    lines!(QN2i[:, i], h_atm/1e3)
end
display(fig)


for i in 1:length(E)
    if all(QN2i[:, i] .== 0)
        println(i)
    end
end
