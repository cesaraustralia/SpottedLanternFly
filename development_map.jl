ENV["RASTERDATASOURCES_PATH"] = "C:\\RasterData"
# load GrowthMaps and other related packages
using GeoData, ModelParameters
using GrowthMaps, Unitful, UnitfulRecipes, Dates, Setfield, Statistics
using StaticArrays

using RasterDataSources, Dates, GeoData
using Plots
using Pipe: @pipe
using DimensionalData: set
using Interpolations, TimeInterpolatedGeoData
using DataFrames, CSV

using GeoData: Between, Y 

include( "functions.jl")

agfact = 1 # aggregation factor (set higher for faster run time)

# download monthly mean data from WorldClim at the 10 m resoltution
months = 1:12
# make two band GeoSeries with Time = 1:12 for tmin and tmax
ser_tavg = map(GeoSeries(WorldClim{Climate}, (:tavg,); month=1:12, res="10m")) do x
    view(x, Band(1))
end
ser_tavg = DimensionalData.set(ser_tavg, :month => Ti(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))

ser_tmntmx = map(GeoSeries(WorldClim{Climate}, (:tmin, :tmax); month=1:12, res="10m")) do x
    view(x, Band(1))
end
ser_tmntmx = DimensionalData.set(ser_tmntmx, :month => Ti(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))

# aggregate
if agfact != 1
    ser_tavg   = GeoData.aggregate(Center(), ser_tavg, agfact)
    ser_tmntmx = GeoData.aggregate(Center(), ser_tmntmx, agfact)
end

tempfunction = "tavg"
devadjust = true
# ############################# START MAIN LOOP #################################
# for tempfunction = ("tavg", "mm")
#     for devadjust = (true, false)


if tempfunction == "tavg"
    ser = ser_tavg
end

if tempfunction == "mm"
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
    ser = meanday_minmaxseries(ser_tmntmx, dates; step=Hour(1), mm_interpolators=(temp=tempinterpolator,))

    # test interpolation visually
    coord = (1861/agfact, 720/agfact, 1) #random point with data
    coord = convert(Tuple{Int, Int, Int}, round.(coord))
    t = 1:(24*12)
    tT = [ser[i][:temp][coord...] for i = t]
    pl = plot(t,  tT; ylims=(0, 40))
    # We need to add 1 for the first hour to be index 1
    # Then add the hour offset we using in the interpolator
    tmintimes = (0:11)*24 .+1 .+5
    tmaxtimes = (0:11)*24 .+1 .+14
    plot!(pl, tmintimes, [ser_tmntmx[i][(:tmin)][coord...] for i = 1:12])
    plot!(pl, tmaxtimes, [ser_tmntmx[i][(:tmax)][coord...] for i = 1:12])

end


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
 
adjustdev = devadjust ? park2009 ./ kreit2020 : [1.0;1.0;1.0;1.0]
    
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
CSV.write("dataout/dev_temp_response_devadjust$devadjust.csv", 
    DataFrame([temp mrates], [:temp, stagenames...]))


# SLF outputs 
# suitability based on development and mortality  based on Kreitman2020 
# monthly map where color represents stage and opacity represents mortality
kwargs = 
    (tspan=DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1),
    initval=SA[0.0, 0.0, 0.0, 0.0, 0.0]) 
tempkey = tempfunction == "tavg" ? :tavg : :temp 
devout = mapgrowth(Layer(tempkey, multi); 
        series=ser,
        kwargs...
);

plot(map(x -> x[1], devout[Ti(8)]))
# Survival model fitted in R using Kreitman2020
# folder="C:/Users/james/Dropbox (cesar)/cesar Team Folder/_Customer Project Files/2003 APBSF/Completed projects/2003CR2 SLF/2. Milestones/Modelling/R/plots/"
folder = "../../R/plots/"

pars = CSV.read(folder*"poisson_survival_model_coefficients.csv", DataFrame)
p = pars.val
ep = CSV.read(folder*"poisson_survival_model_coefficients_egg.csv",
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

hazout = mapgrowth(Layer(tempkey, multihaz); 
    series=ser,
    kwargs...
);

plot(map(x -> x[1], hazout[Ti(8)]))

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

plot(map(x -> x[2], hazout[Ti(6)]))

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
stage,   haz   = getstagehaz(devout, hazout, DateTime(2001,1,1))
stage_s, haz_s = getstagehaz(devout, hazout, DateTime(2001,7,1))


# shift months 
stage[Y(Between(-90, 0)), Ti(1:6)]  = stage_s[Y(Between(-90, 0)), Ti(7:12)]
stage[Y(Between(-90, 0)), Ti(7:12)] = stage_s[Y(Between(-90, 0)), Ti(1:6)]
haz[Y(Between(-90, 0)), Ti(1:6)]    = haz_s[Y(Between(-90, 0)), Ti(7:12)]
haz[Y(Between(-90, 0)), Ti(7:12)]   = haz_s[Y(Between(-90, 0)), Ti(1:6)]


for m = 1:12
    ms = string(m, pad=2)
    GeoData.write("dataout/stage_agfact$(agfact)_tempfun$(tempfunction)_devadjust$(devadjust)_mon$ms.grd", stage[Ti(m)])
    write("dataout/haz_agfact$(agfact)_tempfun$(tempfunction)_devadjust$(devadjust)_mon$ms.grd", haz[Ti(m)])
end

plot(haz_s[Ti(7)])
plot(map(x -> exp(-(x[1])), hazout[Ti(1)]))

aust = X(Between(113, 153)), Y(Between(-43, -11))
usa = X(Between(-125.0, -66.96)), Y(Between(20.0, 50))

l = @layout grid(1,2)
plot(plot(haz_s[Ti(7), usa...]), plot(haz_s[Ti(7), aust...]), layout=l)

plot(exp.(-haz[Ti(8), usa...]))

# ################# END MAIN LOOP #######################
# end
# end
