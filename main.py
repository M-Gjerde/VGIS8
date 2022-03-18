import random

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import math

from matplotlib import cm
from matplotlib.ticker import MaxNLocator


def degToRad(arr):
    return (arr * np.pi) / 180


def latLongToCart(latitude, longitude, r):
    lat = degToRad(latitude)
    lon = degToRad(longitude)

    x = r * np.cos(lat) * np.cos(lon)
    y = r * np.cos(lat) * np.sin(lon)
    z = r * np.sin(lat)

    return x, y, z


def latLonToMercator(lat, lon, lambda0):
    x = lon

    cos_value = np.cos(lat)
    sec_value = np.arccos(cos_value)
    tan_value = np.tan(lat)
    y = np.log(tan_value + sec_value)
    return x, y


def latLonToMillerCylindrical(lat, lon):
    x = lon
    y = (5 / 4) * np.arcsinh(np.tan(4 * degToRad(lat) / 5))
    return x, y


def latLonGenerator(n):
    """
    Sentinel 6 lat lon min/max -66.15 to 66.15
    :param n:
    :return:
    """
    lat = []
    lon = []
    for x in range(n):
        lat.append(random.randint(-66, 66))
        lon.append(random.randint(-180, 180))

    return np.array(lat), np.array(lon)


def read_variable_data(ncid, varname):
    data = np.array(ncid.variables[varname][:])
    longname = ncid.variables[varname].long_name
    standardName = ncid.variables[varname].standard_name
    units = ncid.variables[varname].units
    return data, standardName, longname, units


def surfacePlot(Xs, Ys, Zs):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_trisurf(Xs, Ys, Zs, cmap=cm.jet, linewidth=0)
    fig.colorbar(surf)

    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.zaxis.set_major_locator(MaxNLocator(5))

    fig.tight_layout()
    plt.show()


def scatter2DPlot(Xs, Ys, altitudeData):
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_aspect("auto")
    ax.scatter(Xs, Ys, c=altitudeData)
    plt.show()


def scatterPlot(Xs, Ys, Zs, altitudeData):
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_aspect("auto")
    ax.scatter(Xs, Ys, Zs, c=altitudeData)
    plt.show()


def main():
    # Open the ESA CCI SLA file.
    ncid = nc.Dataset(
        'data/S6A_P4_1A_HR______20220216T084504_20220216T094121_20220311T184648_3377_046_244_122_EUM__OPE_NT_F05.SEN6/measurement.nc')

    # get variables and metadata
    ncid = ncid.groups["data_140"].groups["ku"]

    lat, _, latLongName, latUnits = read_variable_data(ncid, "latitude")
    lon, _, lonLongName, lonUnits = read_variable_data(ncid, "longitude")
    lon = lon[:] - 180
    altitudeData, slaShortName, _, slaUnits = read_variable_data(ncid, "altitude")

    # Get cartesian coordinates
    r = 6371
    # x, y, z = latLongToCart(latData, lonData, r)
    # lat, lon = latLonGenerator(len(altitudeData))
    #x, y, z = latLongToCart(lat, lon, r)
    x, y = latLonToMercator(lat, lon, lat[0])
    # x, y = latLonToMillerCylindrical(lat, lon)
    x = x[::100]
    y = y[::100]
    # z = z[::100]
    altitudeData = altitudeData[::100]
    scatter2DPlot(x, y, altitudeData)
    # scatterPlot(x, y, altitudeData, altitudeData)
    # surfacePlot(x, y, z)


if __name__ == '__main__':
    main()
