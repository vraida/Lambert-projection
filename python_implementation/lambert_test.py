import numpy as np
import pandas as pd
import json
from lib.lambert import Lambert

# Load the reference coordinates (converted with Cartopy)
ref_data = pd.read_csv('../test-data/test-data.csv')

lat_true = ref_data['lat'].to_numpy()
lon_true = ref_data['lon'].to_numpy()
x_true = ref_data['x'].to_numpy()
y_true = ref_data['y'].to_numpy()


# Test if our custom implementation of Lambert projection yields same results as Cartopy:
with open('../projection_constants.json', 'r') as f: const = json.load(f)

lamb = Lambert(
    standard_parallels=(const['phi_1'], const['phi_2']),
    central_latitude=const['phi_F'], central_longitude=const['lambda_F'],
    false_easting=const['E_F'], false_northing=const['N_F'],
    semimajor_axis = const['a'], inverse_flattening = const['f_inv']
)

x_test, y_test = lamb.geographic2cartesian(lat_true, lon_true)
print('Forward direction - largest absolute errors: x = {}, y = {}'.format(np.max(np.abs(x_test-x_true)), np.max(np.abs(y_test-y_true))))

lat_test, lon_test = lamb.cartesian2geographic(x_true, y_true)
print('Backward direction - largest absolute errors: lat = {}, lon = {}'.format(np.max(np.abs(lat_test-lat_true)), np.max(np.abs(lon_test-lon_true))))

