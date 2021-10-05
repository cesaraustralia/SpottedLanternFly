ENV["RASTERDATASOURCES_PATH"] = "C:\\RasterData"
# load GrowthMaps and other related packages
using GeoData, ArchGDAL, NCDatasets, ModelParameters
using GrowthMaps, Unitful, UnitfulRecipes, Dates, Setfield, Statistics
using StaticArrays

using RasterDataSources, Dates, GeoData
using Plots
using Pipe: @pipe
using DimensionalData: set
using Interpolations, TimeInterpolatedGeoData
using DataFrames, CSV

using GeoData: Between

include( "functions.jl")

agfact = 10 # aggregation factor (set higher for faster run time)

# download monthly mean data from WorldClim at the 10 m resoltution
months = 1:12
# make two band GeoSeries with Time = 1:12 for tmin and tmax
ser = map(GeoSeries(WorldClim{Climate}, layerkeys; month=1:12, res="10m")) do x
    view(x, Band(1))
end

ser_tavg = GeoSeries(WorldClim{Climate}, :tavg, month=1:12, res="10m")
ser_tavg = DimensionalData.set(ser_tavg, :month => Ti(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))
ser_tavg = GeoData.aggregate(Center(), ser_tavg, (agfact, agfact, 1); keys=(:tavg,))


# make two band GeoSeries with Time = 1:12 for tmin and tmax
ser_tmntmx = map(GeoSeries(WorldClim{Climate}, (:tmin, :tmax); month=1:12, res="10m")) do x
    view(x, Band(1))
end
ser_tmntmx = DimensionalData.set(ser_tmntmx, :month => Ti(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))
ser_tmntmx = GeoData.aggregate(Center(), ser_tmntmx, 10; 
    keys=(:tmin, :tmax))

# We use a DimensionalData dim instead of a vector of dates because a
# Dim can have a step size with irregular spaced data - here 1 Hour.
dates = Ti(vec([d + h for h in Hour.(0:23), d in index(ser_tmntmx, Ti)]);
    mode=Sampled(Ordered(), Regular(Hour(1)), Intervals())
)


# We use a DimensionalData dim instead of a vector of dates because a
# Dim can have a step size with irregular spaced data - here 1 Hour.
dates = [d + h for d in index(ser_tmntmx, Ti) for h in Hour.(0:23) ]


# Create Min Max specification for intorpolation
# times is a tuple specifying the time of tmin and tmax
# interpmode specifies the interporator

#######new
# Create Min Max interpolator 
# times is a tuple specifying the time of tmin and tmax
tempinterpolator = MinMaxInterpolator((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))

# create an interpolated GeoSeries using the original stacked series, the new dates to interpolate ,
mmseries = meanday_minmaxseries(ser_tmntmx, dates; step=Hour(1), mm_interpolators=(temp=tempinterpolator,))

mmseries[At(DateTime(2001, 11, 1, 2))][:temp] |> plot

# update growth response to link with the interpolated key
growthresponse = Layer(:temp, °C, growthmodel)

# now run simulation using interpolated series and create GeoArray for each
growthrates = mapgrowth(stripparams(growthresponse);
    series=mmseries,
    tspan=DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)
)
plot(growthrates[Ti(1:3:12), Band(1)], clim=(0, 0.15))


######new end

# test interpolation @rafaqz
ser_mm = mmseries
t = 1:(24*12)
tT = [ser_mm[i][:temp][186, 72, 1] for i = t]
pl = plot(t,  tT)
plot!(pl, (0:11)*24, [ser_tmntmx[i][(:tmin)][186, 72, 1] for i = 1:12])
plot!(pl, (0:11)*24, [ser_tmntmx[i][(:tmax)][186, 72, 1] for i = 1:12])

#################### did not get passed here ##############


tempspec = MinMaxSpec((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))
tempfracs = MinMaxFracs(tempspec, Day(1), 0/24)
tempfracs.indices == (tmax=0, tmin=1)
tempfracs.frac == 10/15
tempfracs.interpmode == BSpline(Linear())

tempfracs = MinMaxFracs(tempspec, Day(1), 5/24)
tempfracs.indices == (tmin=1, tmax=1)
tempfracs.frac == 0.0

# Kreitman 2020 parameters (Briere and Linear fitted)
l1 = BriereModel(0.000011, 0.9335, 13.05, 43.81)
l2 = BriereModel(0.000062, 1.9766, 12.43, 35.58)
l3 = BriereModel(0.000073, 6.4742,  9.87, 35.00)
l4 = BriereModel(0.000065, 9.6873, 10.00, 35.00) # typo in orginal paper

# l1 = LinearModel(-0.06069, 0.00529)
# l2 = LinearModel(-0.03526, 0.00379)
# l3 = LinearModel(-0.02121, 0.00250)
# l4 = LinearModel(-0.0103,  0.00164)  # typo in original

# # Park 2015
# l1 = DegreeDayModel(11.13,  270.71)
# l2 = DegreeDayModel(11.13,  491.98)
# l3 = DegreeDayModel(11.13,  618.31)
# l4 = DegreeDayModel(11.13,  907.60)  
# a  = DegreeDayModel(11.13, 1820.65) 

#  Smyers 2020 for linear hatch model
# also similar to Park 2015 thesis
e  = LinearModel(-0.03370, 0.00323) 


multi = MultiRate((e, l1, l2, l3, l4,))

# compare dev rate of L1-L4 at 20C to Park 2009
park2009 = 1 ./ [18.8, 20.9, 20.8, 22.2] 
kreit2020 = GrowthMaps.rate.(Ref(multi), 20)[2:5]
adjustdev = park2009 ./ kreit2020

# adjust development by park factor
l1 = BriereModel(adjustdev[1]*0.000011, 0.9335, 13.05, 43.81)
l2 = BriereModel(adjustdev[2]*0.000062, 1.9766, 12.43, 35.58)
l3 = BriereModel(adjustdev[3]*0.000073, 6.4742,  9.87, 35.00)
l4 = BriereModel(adjustdev[4]*0.000065, 9.6873, 10.00, 35.00) # typo in orginal paper
multi = MultiRate((e, l1, l2, l3, l4,))


temp = collect(-10:0.001:40)
mrates = GrowthMaps.rate.(Ref(multi), temp)
mrates = [mrates[i][j] for i = 1:length(mrates), j = 1:5]
pl = plot(temp, mrates)



stagenames = (:egg,:L1,:L2,:L3,:L4)
CSV.write("dataout/dev_temp_response.csv", 
    DataFrame([temp mrates], [:temp, stagenames...]))


# SLF outputs 
# suitability based on development and mortality  based on Kreitman2020 
# monthly map where color represents stage and opacity represents mortality
kwargs = 
    (tspan=DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1),
    initval=SA[0.0, 0.0, 0.0, 0.0, 0.0]) 

devout_tavg = mapgrowth(Layer(:tavg, multi); 
    series=ser_tavg,
    kwargs...
);

devout_mm = mapgrowth(Layer(:temp, multi); 
    series=ser_mm,
    kwargs...
);


plot(map(x -> x[1], devout_tavg[Ti(8)]))
plot(map(x -> x[1], devout_mm[Ti(8)]))


# Survival model fitted in R using Kreitman2020
pars = CSV.read("../../R/plots/poisson_survival_model_coefficients.csv", DataFrame)
const p = pars.val
ep = CSV.read("../../R/plots/poisson_survival_model_coefficients_egg.csv",
    DataFrame)


he = EggHazModel(ep[1,2],ep[2,2], ep[1,3],ep[2,3]) 
h1 = HazModel(p[1], p[2], p[3], p[4], p[5] )
h2 = h1
h3 = h1
h4 = h1

multihaz = MultiRate((he, h1, h2, h3, h4,))

plot(temp, rate.(Ref(he), temp))
rate.(Ref(h4), 3)

temp = collect(-15:0.01:40)

mrates = GrowthMaps.rate.(Ref(multihaz), temp)
mrates = [mrates[i][j] for i = 1:length(mrates), j = 1:5]
pl = plot(temp, mrates, ylim=(0,0.2))

CSV.write("dataout/hazard_temp_response.csv", 
    DataFrame([temp mrates], [:temp, stagenames...]))

hazout_tavg = mapgrowth(Layer(:tavg, multihaz); 
    series=ser_tavg,
    kwargs...
);

hazout_mm = mapgrowth(Layer(:temp, multihaz); 
    series=ser_mm,
    kwargs...
);


function getstagerate(rates, stage) 
    # add modulo function for stage looping 
    # stage = mod(stage,  length(rates))
    if isnan(stage)
        return NaN
    elseif (stage + 1e-9) > length(rates)
        # return 0
        # assume haz accumulates in adult as nymphs
        return rates[length(rates)]
    elseif (stage + 1e-9) < 0
        return 0
    else
       return  rates[convert(Int, ceil(stage + 1e-9))]
    end
end

plot(map(x -> x[2], hazout_mm[Ti(6)]))

function getstagehaz(devout, hazout, startdate)
    stage = map(x -> x[1]*0.0, devout)
    haz   = map(x -> x[1]*0.0, hazout)
    days = DateTime(2001, 1, 1):Day(1):DateTime(2001, 12, 31)
    start = dayofyear(startdate)
    n = length(days)
    dayser = start == 1 ? days : days[vcat(start:n, 1:(start-1))];
    for d = dayser
        m = month(d)
        stage[Ti(m)] += 
            map(getstagerate, devout[Ti(m)], stage[Ti(m)] )
        haz[Ti(m)] += 
            map(getstagerate, hazout[Ti(m)], stage[Ti(m)] )    
        if day(d + Day(1)) == 1  && 
            sum(map(x -> isnan(x) ? 0 : x, stage[Ti(month(d + Day(1)))])) == 0
            m1 = month(d + Day(1))
            stage[Ti(m1)] = stage[Ti(m)] 
            haz[Ti(m1)]   = haz[Ti(m)] 
        end
    end
    return (stage, haz)
end

# run simulation for north and south hemishere starting from winter
stage_tavg,  haz_tavg  = getstagehaz(devout_tavg, hazout_tavg, DateTime(2001,1,1))
stage_tavgs, haz_tavgs = getstagehaz(devout_tavg, hazout_tavg, DateTime(2001,7,1))
stage_mm,  haz_mm  = getstagehaz(devout_mm, hazout_mm, DateTime(2001,1,1))
stage_mms, haz_mms = getstagehaz(devout_mm, hazout_mm, DateTime(2001,7,1))


# shift months 
stage_tavg[Lat(Between(-90, 0)), Ti(1:6)]=stage_tavgs[Lat(Between(-90, 0)), Ti(7:12)]
stage_tavg[Lat(Between(-90, 0)), Ti(7:12)]=stage_tavgs[Lat(Between(-90, 0)), Ti(1:6)]
haz_tavg[Lat(Between(-90, 0)), Ti(1:6)]   = haz_tavgs[Lat(Between(-90, 0)), Ti(7:12)]
haz_tavg[Lat(Between(-90, 0)), Ti(7:12)]  = haz_tavgs[Lat(Between(-90, 0)), Ti(1:6)]

stage_mm[Lat(Between(-90, 0)), Ti(1:6)]   = stage_mms[Lat(Between(-90, 0)), Ti(7:12)]
stage_mm[Lat(Between(-90, 0)), Ti(7:12)]  = stage_mms[Lat(Between(-90, 0)), Ti(1:6)]
haz_mm[Lat(Between(-90, 0)), Ti(1:6)]   = haz_mms[Lat(Between(-90, 0)), Ti(7:12)]
haz_mm[Lat(Between(-90, 0)), Ti(7:12)]  = haz_mms[Lat(Between(-90, 0)), Ti(1:6)]


for m = 1:12
    ms = string(m, pad=2)
    write("dataout/stage_tavg$ms.grd", GRDarray, stage_tavg[Ti(m)])
    write("dataout/haz_tavg$ms.grd"  , GRDarray, haz_tavg[Ti(m)])
    write("dataout/stage_mm$ms.grd"  , GRDarray, stage_mm[Ti(m)])
    write("dataout/haz_mm$ms.grd"    , GRDarray, haz_mm[Ti(m)])
end

plot(haz_mms[Ti(7)])
plot(map(x -> exp(-(x[1])), hazout_mm[Ti(1)]))

aust = Lon(Between(113, 153)), Lat(Between(-43, -11))
usa = Lon(Between(-125.0, -66.96)), Lat(Between(20.0, 50))

l = @layout grid(1,2)
plot(plot(haz_mms[Ti(7), aust...]), plot(haz_mms[Ti(7), aust...]), layout=l)


## load rainfall 
# ser_immut = series(WorldClim{Climate}, (:tavg,:prec); res="10m")
# dates = (DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1))
# ser = [GeoStack(ser_immut[i][:tavg], float(ser_immut[i][:prec]); 
#     keys = (:tavg, :prec)) for i = 1:12] 
# ser = GeoSeries(ser, (Ti(dates),))
# for i in 1:12    
#     ser[i][:prec] .= map(x -> x < 0 ? 0.0 : x, ser[i][:prec])
# end
# plot(ser[3][:prec])
# plot(ser[DateTime(2001, 4, 1)][:tavg], clim = (0, 40))
