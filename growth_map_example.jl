
# load GrowthMaps and other related packages
using GrowthMaps, Unitful, UnitfulRecipes, Dates, Setfield, Statistics
using GeoData, ArchGDAL, NCDatasets, ModelParameters
using Unitful: °C, K, cal, mol, mm
# using HDF5 # for SMAP

# Set SchoolfieldIntrinsicGrowth model parameters including fields for units and bounds for fitting
ψ_EA =12.254cal/K/mol
ΔH_A = Param(11822.57; units=cal/mol, bounds=(1e5, 2e5))
ΔH_H = Param(90678.50; units=cal/mol, bounds=(0.5e5, 1e5))
ΔH_L = Param(-6603.414; units=cal/mol, bounds=(0.5e5, 1e5))
p27 = Param(0.75; bounds=(0.0, 1.0)) # this is higher than 0.17 as
T_ref = K(27.0°C)
S_H = 295.08cal/mol/K
S_L = -19.59cal/mol/K
Thalf_L = Param(ustrip(K, withunits(ΔH_L)/S_L); units=K, bounds=(290, 350))
Thalf_H = Param(ustrip(K, withunits(ΔH_H)/S_H); units=K, bounds=(290, 350))

growthmodel = SchoolfieldIntrinsicGrowth(p27, ΔH_A, ΔH_L, Thalf_L, ΔH_H, Thalf_H, T_ref)

# Link growth model to hypothetical layer of temperature data.
growthresponse = Layer(:tavg, °C, growthmodel)

# Test growth at 20°C
GrowthMaps.rate(stripparams(growthresponse), K(20.0°C))

# plot temperature response from 0 - 40°C
temprangeC = collect(1.0:1.0:40.0)°C
temprange = K.(temprangeC)
dev = (x -> GrowthMaps.rate(stripparams(growthresponse), x)).(temprange)
p = plot(temprangeC, dev; label=false, ylabel="growth rate (1/d)", xlab = "temperature")


# now run simulation using mean monthly temp and create GeoArray
modelkwargs = (
    series=agser,
    tspan=DateTime(2001, 1, 1):Year(1):DateTime(2001, 12, 1),
);

growthrates = mapgrowth(stripparams(growthresponse);
    modelkwargs...)
plot(growthrates[Ti(1)], clim=(0, 0.15))

# incorporate stress
coldthresh = Param(ustrip(K, 13.0°C); units=K, bounds=(-10, 20))
coldmort = Param(-0.2; units=K^-1, bounds=(-1, 0))
coldstress = Layer(:tavg, °C, LowerStress(coldthresh, coldmort))

heatthresh = 40.0°C |> K
heatmort = Param(-0.02; units=K^-1, bounds=(-1, 0))
heatstress = Layer(:tavg, °C, UpperStress(heatthresh, heatmort))

drythresh = Param(100.0; bounds=(0, 1000))
drymort = Param(-0.0001;  bounds=(-1, 0))
drystress = Layer(:prec, LowerStress(drythresh, drymort));

model = stripparams((growthresponse, coldstress, heatstress, ))

growthrates = mapgrowth(model; modelkwargs...)

plot(growthrates[Ti(1)], axis=false)
plot(growthrates[Ti(1)], clim=(0, 0.15), axis=false)
    