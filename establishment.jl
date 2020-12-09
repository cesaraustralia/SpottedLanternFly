

### TEST WEATHER ####
using Interpolations, TimeInterpolatedGeoData, RasterDataSources, Test, Dates, GeoData
using RasterDataSources: Values, SoilMoisture, Upper, Lower

# Weather time-series
dates = (Date(2001), Date(2002))
layers = (:tmin, :tmax)
ser = series(WorldClim{Weather}; layers=layers, dates=dates, window=(Lat(600:1900), Lon(600:1900)))
ser[Date(2001, 1)][:tmax]
# Select Australia, using regular lat/lon selectors
A = geoarray(WorldClim{Weather}, :prec, Date(2001, 05); mappedcrs=EPSG(4326))
A[Lat(Between(-10, -45)), Lon(Between(110, 160))]

# Weather time-series
tempspec = MinMaxSpec((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))
specs = (temp=tempspec,)

mmseries = minmaxseries(ser, DateTime(2000, 12, 31):Hour(1):DateTime(2001, 12, 31), (), specs)
iseries = interpseries(ser, DateTime(2019, 1):Hour(1):DateTime(2019, 12, 31), (tmax=BSpline(Linear()), tmin=BSpline(Linear()),))

mmseries[DateTime(2019, 1, 5, 6)][:temp]


####################

# WorldClim Weather

# Weather time-series
ser = series(WorldClim{Climate}; layers=(:tmin,), resolution="10m", months=Jan:March)
ser[Jan][:tmin]
# Select Australia, using regular lat/lon selectors
A = geoarray(WorldClim{Climate}, :tmin, "10m", Feb; mappedcrs=EPSG(4326))
A[Lat(Between(-10, -45)), Lon(Between(110, 160))]

using Plots
A[Lat(Between(-10, -45)), Lon(Between(110, 160))] |>
   plot


layers = (:tmin, :tmax)
ser = series(WorldClim{Climate}; layers = layers, resolution="10m", months=Jan:Dec)

dates = (Date(2001), Date(2002))
ser[Jan][:tmax]

# Weather time-series
index(ser, Ti)
ser[Jan][:tmin]
ser[Jan][:tmax]

mmseries = minmaxseries(ser, DateTime(2019, 1):Hour(1):DateTime(2019, 12, 31), (), specs)

iseries = interpseries(ser, DateTime(2019, 1):Hour(1):DateTime(2019, 12, 31), (tmax=BSpline(Linear()), tmin=BSpline(Linear()),))
awapstack = stack(AWAP; date=DateTime(2019, 1, 5))
mmseries[DateTime(2019, 1, 5, 6)][:temp]
