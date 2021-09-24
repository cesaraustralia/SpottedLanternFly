# ENV["RASTERDATASOURCES_PATH"] = "\MyDataLocation"
using Interpolations, TimeInterpolatedGeoData, RasterDataSources, Test, Dates, GeoData
using Plots, Unitful, GrowthMaps
# load GrowthMaps and other related packages
using GrowthMaps, Unitful, UnitfulRecipes, Dates, Setfield, Statistics
using GeoData, ModelParameters
using Unitful: °C, K, cal, mol
# using HDF5 # for SMAP

# download montly mean tmax and tmin data from WorldClim at the 10 m resoltution
months = 1:12
layerkeys = (:tmin, :tmax)
ser_immut = GeoSeries(WorldClim{Climate}, layerkeys; month=1:12, res="10m")

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
pl = plot(temprangeC, dev; label=false, ylabel="growth rate (1/d)", xlab = "temperature")

# make two band GeoSeries with Time = 1:12 for tmin and tmax
ser = view(GeoSeries(WorldClim{Climate}, layerkeys; month=1:12, res="10m"), Band(1))
ser = DimensionalData.set(ser, :month => Ti(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))

# We use a DimensionalData dim instead of a vector of dates because a
# Dim can have a step size with irregular spaced data - here 1 Hour.
dates = [d + h for d in index(ser, Ti) for h in Hour.(0:23) ]

# Create Min Max interpolator 
# times is a tuple specifying the time of tmin and tmax
tempinterpolator = MinMaxInterpolator((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))

# create an interpolated GeoSeries using the original stacked series, the new dates to interpolate ,
mmseries = meanday_minmaxseries(ser, dates; step=Hour(1), mm_interpolators=(temp=tempinterpolator,))

mmseries[At(DateTime(2001, 11, 1, 2))][:temp] |> plot

# update growth response to link with the interpolated key
growthresponse = Layer(:temp, °C, growthmodel)

# now run simulation using interpolated series and create GeoArray for each
growthrates = mapgrowth(stripparams(growthresponse);
    series=mmseries,
    tspan=DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)
)
pyplot()
plot(growthrates[Ti(1:3:12)], clim=(0, 0.15))
