import random
import glob

import matplotlib.animation
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import math
from datetime import datetime, timedelta
from matplotlib import cm
from matplotlib.ticker import MaxNLocator

earthRadius = r = 6371


def applyWindow(window, latitude, longitude, data):
    latCondition = (latitude > window[1]['lat']) & (latitude < window[0]['lat'])
    longCondition = (longitude > window[0]['lon']) & (longitude < window[1]['lon'])
    combinedCondition = latCondition & longCondition

    latitude = np.extract(combinedCondition, latitude)
    longitude = np.extract(combinedCondition, longitude)
    data = np.extract(combinedCondition, data)

    return latitude, longitude, data


def degToRad(arr):
    return (arr * np.pi) / 180


# function to get unique values
def unique(list1):
    # initialize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
    return unique_list


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


def read_pass_number(ncid):
    data = []
    for nc in ncid:
        data.append(nc.parent.pass_number)
    return np.array(data)


def read_variable_data(ncid, varname):
    if len(ncid) > 1:
        row = []
        for x in range(len(ncid)):
            row.append(ncid[x].variables[varname][:])
        return row
    else:
        data = np.array(ncid.variables[varname][:])
        longname = ncid.variables[varname].long_name
        standardName = ncid.variables[varname].standard_name
        units = ncid.variables[varname].units
        return data


def surfacePlot(Xs, Ys, Zs):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_trisurf(Xs, Ys, Zs, cmap=cm.jet, linewidth=1)
    fig.colorbar(surf)

    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.zaxis.set_major_locator(MaxNLocator(5))

    """
    # Plot plane at origion
    x = np.linspace(-10, 10, 100)
    y = np.linspace(-10, 10, 100)

    x, y = np.meshgrid(x, y)
    eq = 1 * x + 1 * y

    ax.plot_surface(x, y, eq)
    """

    fig.tight_layout()
    return fig, ax


def scatter2DPlot(Xs, Ys, altitudeData):
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_aspect("auto")
    ax.scatter(Xs, Ys, c=altitudeData)
    plt.show()


def scatterPlot(Xs, Ys, Zs, altitudeData, size=0.8):
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_aspect("auto")

    ax.scatter(Xs, Ys, Zs, c=altitudeData, s=size)
    return fig, ax


def scatterCities(ax):
    r = 6371
    ## Aalborg
    X, Y, Z = latLongToCart(57.0469106, 9.9167348, r)
    ax.scatter(X, Y, 0, c="blue", s=5)
    ## CPH
    ## 55.6713442,12.5237847
    X, Y, Z = latLongToCart(55.6713442, 12.5237847, r)
    ax.scatter(X, Y, 0, c="red", s=5)
    ## Esbjerg
    ## 55.5225969,8.3934395
    X, Y, Z = latLongToCart(55.5225969, 8.3934395, r)
    ax.scatter(X, Y, 0, c="green", s=5)


def loadDataFiles():
    # Open the ESA CCI SLA file.
    path = "data/l2_caspian/*/*STD*.nc"
    files = glob.glob(path)
    ncid = []
    for file in files:
        ncid.append(nc.Dataset(file).groups["data_01"])
    return ncid


def getIndicesOfItem(lst, item):
    return [i for i, x in enumerate(lst) if x == item]


def getLatLonBoundingBox():
    # lon = []
    # altitudeData = []
    # for i in altitudeData_raw:
    # altitudeData.append((i[:] - 1300000.0) * 0.0001)# - 180)

    ### LAT LON SQUARE
    # View online
    # http://bboxfinder.com/#53.173110,6.811520,57.973150,14.919430
    # AALBORG
    # Upper left 58.06,7.72
    # Lower left 54, 12.78
    # http://bboxfinder.com/#54,7.72,58.06,12.78

    # CASPBIAN SEA
    offset = 0
    latMin = 33.66 - offset
    latMax = 53.53 + offset
    lonMin = 39.16 - offset
    lonMax = 61.78 + offset
    # AAlborg
    # offset = 0
    # latMin = 54 - offset
    # latMax = 58.06 + offset
    # lonMin = 7.72 - offset
    # lonMax = 12.78 + offset
    print(str(str(latMin) + "," + str(lonMin) + "," + str(latMax) + "," + str(lonMax).strip()))
    return latMin, latMax, lonMin, lonMax


def main():
    # get variables and metadata
    print("Reading data")
    ncid = loadDataFiles()
    latList = read_variable_data(ncid, "latitude")
    lonList = read_variable_data(ncid, "longitude")
    altitudeDataList = read_variable_data(ncid, "mean_sea_surface_sol1")
    passNumbers = read_pass_number(ncid)
    time = read_variable_data(ncid, "time")

    # Remove values that are outside the specified window
    denmark = [{'lat': 58.06, 'lon': 7.72},
               {'lat': 54.0, 'lon': 12.78}]

    ## Create A list of every repeated pass.
    # Orbits object has indices for which every data point repeats its pass over earth ground tracks.
    uniquePassNumbers = unique(passNumbers)
    orbits = []
    for orbitPass in uniquePassNumbers:
        indices = getIndicesOfItem(passNumbers, orbitPass)
        orbits.append(indices)

    latMin, latMax, lonMin, lonMax = getLatLonBoundingBox()

    x = [[]]
    y = [[]]
    z = [[]]
    h = [[]]

    for i in range(len(orbits[0])):
        x.append([])
        y.append([])
        z.append([])
        h.append([])

    for orbit in orbits:
        for index, groundPass in enumerate(orbit):
            latArr = latList[groundPass]
            lonArr = lonList[groundPass]
            altArr = altitudeDataList[groundPass]

            for i in range(len(latArr)):
                if latMin < latArr[i] < latMax:
                    if lonMin < lonArr[i] < lonMax:
                        lat = latArr[i]
                        lon = lonArr[i]
                        altitude = altArr[i]
                        x_, y_, z_ = latLongToCart(lat, lon, earthRadius)
                        # x_, y_ = latLonToMillerCylindrical(lat, lon)
                        x[index].append(x_)
                        y[index].append(y_)
                        z[index].append(z_)
                        h[index].append(altitude)

    # scatter2DPlot(x, y, h)
    # surfacePlot(x[0], y[0], h[0])
    # print("Plotting data")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    title = ax.set_title('3D Test')

    xPlot, yPlot, hPlot = [], [], []
    graph, = ax.plot(xPlot, yPlot, hPlot, linestyle="", marker="o")

    n = 2
    print("Filtering some points")
    for i in range(len(x)):
        x[i] = x[i][::n]
        y[i] = y[i][::n]
        z[i] = z[i][::n]
        h[i] = h[i][::n]

    # plot = [ax.plot_trisurf(x[0], y[0], h[0], cmap=cm.jet, linewidth=1)]

    plt.xlim(min(x[0]), max(x[0]))
    plt.ylim(min(y[0]), max(y[0]))
    ax.set_zlim(20, -45)

    def animate(i):
        ax.clear()
        ax.plot_trisurf(np.array(x[i]), np.array(y[i]), np.array(h[i]), cmap=cm.jet, linewidth=1)
        # graph.set_data(x[i], y[i])
        # graph.set_color(h[i])
        # graph.set_3d_properties(h[i])
        print("Plotting i:", i)
        return title, graph,

    ani = matplotlib.animation.FuncAnimation(fig, animate,
                                             frames=len(orbits[0]) - 1, interval=1000, repeat=True)

    plt.show()

    ## AAlborg coordinates
    # lat, lon: 57.0469106,9.9167348


# surfacePlot(x, y, c)


if __name__ == '__main__':
    main()
