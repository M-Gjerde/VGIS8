import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import os
import cartopy.crs as crs
import cartopy.feature as feat
from scipy.interpolate import griddata
import warnings
from scipy.interpolate.rbf import Rbf

warnings.filterwarnings("ignore")

""" Math functions """


# Normalize an ndarray
def normalize(x):
    mag = np.sqrt(np.sum(x * x))
    return x / mag


# Convert degrees to radians, works with both ndarray and scalar values
def degToRad(a):
    return (a * np.pi) / 180


""" For projection """


# Calculate the miller projection only in Y, works with both ndarray and scalar values
def millerY(latitude):
    y_prime = 1.25 * np.log(np.tan((np.pi / 4.0) + (0.4 * latitude)))
    return y_prime


# Calculate the miller projection in both Y and X coordinates
def millerCombined(latitude, longitude):
    lat = degToRad(latitude)
    lon = degToRad(longitude)

    x_prime = lon
    y_prime = millerY(lat)

    return x_prime, y_prime


# Convert a list of latitude, longitude and altitude values to u, v coordinates
def latLongToMiller(size, window, latitude, longitude):
    height, width = size[0], size[1]

    x_prime, y_prime = millerCombined(latitude, longitude)

    topLat = millerY(degToRad(window[0]['lat']))
    bottomLat = millerY(degToRad(window[1]['lat']))
    leftLong = degToRad(window[0]['lon'])
    rightLong = degToRad(window[1]['lon'])

    normU = width / (rightLong - leftLong)
    normV = height / (topLat - bottomLat)

    i = (x_prime - leftLong) * normU
    j = (topLat - y_prime) * normV

    return i.astype(np.int64), j.astype(np.int64)


# Given a projection of orbited values create an interpolated image
def interpolate(img, size, interpolation='linear'):
    output = np.copy(img)

    # Extract points that should be interpolated
    condition = img > 0
    coords = np.where(img > 0)
    ix = np.where(img <= 0)
    vals = np.extract(condition, img)

    interpolated = griddata(coords, vals, ix, method=interpolation, fill_value=-1)

    for i in range(ix[0].shape[0]):
        output[ix[0][i], ix[1][i]] = interpolated[i]

    # Remove values that could not be interpolated
    temp = np.extract(output != -1, output)
    cols = int(temp.shape[0] / size[1])
    output = temp.reshape((cols, size[1]))

    output[output > 500] = 500
    output[output < 0] = 0

    return output


def interpolateRBF(x, y, z, func='gaussian'):
    radialBasis = Rbf(x, y, z, function=func)

    x_new = np.linspace(np.amin(x), np.amax(x), x.shape[0])
    y_new = np.linspace(np.amin(y), np.amax(y), y.shape[0])

    x_grid, y_grid = np.meshgrid(x_new, y_new)
    z_new = radialBasis(x_grid.ravel(), y_grid.ravel()).reshape(x_grid.shape)

    return x_new, y_new, z_new


""" For reading and modifying the data """


def read_variable_data(ncid, varname):
    data = ncid.variables[varname][:]
    longname = ncid.variables[varname].long_name
    standardName = ncid.variables[varname].name
    units = ncid.variables[varname].units
    return data, standardName, longname, units


def gatherData(path, variable):
    latitude = np.array([])
    longitude = np.array([])
    data = np.array([])

    for foldername in os.listdir(path):
        if foldername.endswith(".SEN6"):
            for filename in os.listdir(path + foldername):
                if filename.endswith(".nc") and 'STD' in filename:
                    # Open the ESA CCI SLA file.
                    ncid = nc.Dataset(path + '/' + foldername + '/' + filename)

                    # get variables and metadata
                    ncid = ncid["data_01"]

                    latData, _, latLongName, latUnits = read_variable_data(ncid, "latitude")
                    lonData, _, lonLongName, lonUnits = read_variable_data(ncid, "longitude")
                    lonData = np.array(lonData)
                    altitudeData, slaShortName, _, slaUnits = read_variable_data(ncid, variable)

                    latitude = np.concatenate((latitude, latData))
                    longitude = np.concatenate((longitude, lonData))
                    data = np.concatenate((data, altitudeData))

                    continue
                else:
                    continue
            continue
        else:
            continue

    return latitude, longitude, data


def applyWindow(window, latitude, longitude, data):
    latCondition = (latitude > window[1]['lat']) & (latitude < window[0]['lat'])
    longCondition = (longitude > window[0]['lon']) & (longitude < window[1]['lon'])
    combinedCondition = latCondition & longCondition

    latitude = np.extract(combinedCondition, latitude)
    longitude = np.extract(combinedCondition, longitude)
    data = np.extract(combinedCondition, data)

    return latitude, longitude, data


def windowToAxes(window):
    return [window[0]['lon'], window[1]['lon'], window[1]['lat'], window[0]['lat']]


def main():
    # Initialize variables
    world = [{'lat': 90.0, 'lon': 0.0},
             {'lat': -90.0, 'lon': 360.0}]

    denmark = [{'lat': 58.06, 'lon': 7.72},
               {'lat': 54.0, 'lon': 12.78}]

    australia = [{'lat': -3.32, 'lon': 107.93},
                 {'lat': -45.19, 'lon': 160.49}]

    axes = windowToAxes(australia)
    size = [1000, 1000]

    print("Loading data... \n")

    # Get path of the satellite
    latitude, longitude, data = gatherData("data/", "mean_sea_surface_sol1")

    # Remove values that are outside the specified window
    latitude, longitude, data = applyWindow(australia, latitude, longitude, data)

    # Binarize data

    # lon2d, lat2d, data2d = interpolateRBF(longitude.astype(np.float32), latitude.astype(np.float32), data, 'thin_plate')

    print("Loading data: succesful.\nCreating an image...")

    plt.figure(figsize=(6, 5))
    ax = plt.axes(projection=crs.Mollweide())
    ax.add_feature(feat.LAND, zorder=100, edgecolor='k')
    ax.set_extent(axes)
    ax.coastlines()
    # ax.contourf(longitude, latitude, data2d, transform=crs.PlateCarree(), cmap='jet')
    plt.scatter(x=longitude, y=latitude,
                color='red',
                s=1,
                transform=crs.PlateCarree())
    ax.gridlines(draw_labels=False)
    plt.show()


if __name__ == '__main__':
    main()
