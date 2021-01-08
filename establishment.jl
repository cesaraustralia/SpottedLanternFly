# First, to specify the directory in which the data is to be downloaded, modify your startup.jl file located in your Julia install directory (e.g. Julia\Julia 1.5.2\etc\julia\startup.jl) to include the following line:
# test vscode

# ENV["RASTERDATASOURCES_PATH"] = "\MyDataLocation"
using TimeInterpolatedGeoData
using Interpolations, RasterDataSources, Dates, GeoData, ArchGDAL
# using RasterDataSources: Values, SoilMoisture, Upper, Lower 
using Plots
using Pipe: @pipe

# download monthly mean data from WorldClim at the 10 m resoltution
months = 1:12
layers = (:tavg, :prec)
getraster(WorldClim{Climate}, layers, months, "10m")
getraster(WorldClim{Climate}, :wind; month=1, res="10m") 
ser_immut = series(WorldClim{Climate}, layers; res="10m")

dates = (DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1))
ser = [GeoStack(ser_immut[i][:tavg], float(ser_immut[i][:prec]); 
    keys = (:tavg, :prec)) for i = 1:12] 
ser = GeoSeries(ser, (Ti(dates),))
for i in 1:12    
    ser[i][:prec] .= map(x -> x < 0 ? 0.0 : x, ser[i][:prec])
end

map(x -> x < 0 ? 0.0 : float(x), ser1[i][:prec]).data
plot(ser[1][:prec])

plot(ser[DateTime(2001, 4, 1)][:tavg], clim = (0, 40))

agser = GeoData.aggregate(Center(), ser, (10, 10, 1); keys=layers)

@pipe agser[DateTime(2001, 4, 1)][:tavg] |>
    plot(_, )

@pipe agser[DateTime(2001, 4, 1)][:prec] |>
    plot(_, )   

agser[DateTime(2001, 4, 1)][:tavg].data
agser[DateTime(2001, 4, 1)][:prec].data

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






# ##### SMAP
# smapfolder = "K:/SMAP/SMAP_L4_SM_gph_v4"
# smapfolder = "K:/SMAP/SMAP_monthly_midpoint_subset"
# # This does the same as the above
# days = (1, )
# seriesall = SMAPseries(smapfolder)
# series = seriesall[Where(t -> (t >= DateTime(year) && t < DateTime(year)) && dayofmonth(t) in days)]
#
# #' We can plot a layer from a file at some date in the series:
# seriesall[2][:surface_temp] |> plot

####
