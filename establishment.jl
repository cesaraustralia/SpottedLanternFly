
using Interpolations, TimeInterpolatedGeoData, RasterDataSources, Test, Dates, GeoData
using RasterDataSources: Values, SoilMoisture, Upper, Lower

# download montly mean tmax and tmin data from WorldClim at the 10 m resoltution
layers = (:tmin, :tmax)
months = 1:12
download_raster(WorldClim{Climate}, layers; resolution = "10m", month=months)

# make two band GeoSeries with Time = 1:12 for tmin and tmax
ser = series(WorldClim{Climate}; layers=layers)
ser = set(ser, Ti=(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))

# We use a DimensionalData dim instead of a vector of dates because a
# Dim can have a step size with irregular spaced data - here 1 Hour.
dates = Ti(vec([d + h for h in Hour.(0:23), d in index(ser, Ti)]);
    mode=Sampled(Ordered(), Regular(Hour(1)), Intervals())
)

# Create Min Max specification for intorpolation
# times is a tuple specifying the time of tmin and tmax
# interpmode specifies the interporator
tempspec = MinMaxSpec((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))

#  Create an interpolated GeoSeries using the original stacked series, the new dates to interpolate ,
mmseries = meandayminmaxseries(ser, dates, (), (temp=tempspec,))
mmseries[DateTime(2001, 1, 1, 6)][:temp]
mmseries[DateTime(2001, 4, 1, 1)][:temp]

using Plots
mmseries[DateTime(2001, 11, 1, 2)][:temp] |> plot


# load GrowthMaps and other related packages
using GrowthMaps, Unitful, UnitfulRecipes, Dates, Setfield, Statistics
using GeoData, ArchGDAL, NCDatasets, ModelParameters
using Unitful: °C, K, cal, mol
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
growthresponse = Layer(:temp, °C, growthmodel)

# Test growth at ~20C
GrowthMaps.rate(stripparams(growthresponse), K(20.0°C))

temprangeC = collect(1.0:1.0:40.0)°C
temprange = K.(temprangeC)
# GrowthMaps.rate.(Ref(growthresponse), temprange)
dev = (x -> GrowthMaps.rate(stripparams(growthresponse), x)).(temprange)
p = plot(temprangeC, dev; label=false, ylabel="growth rate (1/d)", xlab = "temperature")


growthrates = mapgrowth(stripparams(growthresponse);
    series=mmseries,
    tspan=DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)
)
plot(growthrates[Ti(1)], clim=(0, 0.15))





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
