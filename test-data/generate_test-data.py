import numpy as np
import pandas as pd
import json
import cartopy.crs as ccrs

# Generate random geographic coordinates (mostly in Austria)
lat_min = 46.372132
lat_max = 49.020703
lon_min = 9.530768
lon_max = 17.160749

number = 1000

lat = np.random.uniform(lat_min, lat_max, number)
lon = np.random.uniform(lon_min, lon_max, number)


# Lambert projection with Cartopy: https://scitools.org.uk/cartopy/docs/latest/
with open('../projection_constants.json', 'r') as f: const = json.load(f)

projection_lambert = ccrs.LambertConformal(
    standard_parallels=(const['phi_1'], const['phi_2']),
    central_latitude=const['phi_F'], central_longitude=const['lambda_F'],
    false_easting=const['E_F'], false_northing=const['N_F'],
    globe = ccrs.Globe(semimajor_axis = const['a'], inverse_flattening = const['f_inv'])
)
projection_geodetic = ccrs.Geodetic(globe=ccrs.Globe(ellipse='WGS84'))
cartesian_coordinates = projection_lambert.transform_points(projection_geodetic, lon, lat)
x = cartesian_coordinates[:,0]
y = cartesian_coordinates[:,1]


# Save the result as a reference ground-truth
pd.DataFrame({'lat':lat, 'lon':lon, 'x':x, 'y':y}).to_csv('test-data.csv', index=False)
