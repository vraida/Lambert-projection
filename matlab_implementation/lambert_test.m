clear; clc; close all; addpath inc

% Load the reference coordinates (converted with Cartopy)
ref_data = readtable('../test-data/test-data.csv');

lat_true = ref_data.lat;
lon_true = ref_data.lon;
x_true = ref_data.x;
y_true = ref_data.y;

% Test if our custom implementation of Lambert projection yields same results as Cartopy:
const = jsondecode(fileread('../projection_constants.json'));

lamb = Lambert(...
    [const.phi_1, const.phi_2],  ...  % Latitude of the first and second standard parallel
    const.phi_F, const.lambda_F, ...  % Central latitude and longitude
    const.E_F, const.N_F,        ...  % False easting and northing
    const.a, const.f_inv         ...  % Ellipsoid: Semi-major axis and inverse flattening
);

[x_test, y_test] = lamb.geographic2cartesian(lat_true, lon_true);
fprintf('Forward direction - largest absolute errors: x = %f, y = %f\n', max(abs(x_test-x_true)), max(abs(y_test-y_true)));

[lat_test, lon_test] = lamb.cartesian2geographic(x_true, y_true);
fprintf('Backward direction - largest absolute errors: lat = %f, lon = %f\n', max(abs(lat_test-lat_true)), max(abs(lon_test-lon_true)));