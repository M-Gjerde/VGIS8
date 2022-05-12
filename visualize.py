import glob
import random

import cartopy.crs as crs
import cartopy.feature as feat
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from matplotlib import animation
from scipy.interpolate import Rbf
from numpy import inf
from datetime import datetime, timedelta
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


def interpolateRBF(x, y, z, func='linear'):
    radialBasis = Rbf(x, y, z, function=func)

    x_new = np.linspace(np.amin(x), np.amax(x), x.shape[0])
    y_new = np.linspace(np.amin(y), np.amax(y), y.shape[0])

    x_grid, y_grid = np.meshgrid(x_new, y_new)
    z_new = radialBasis(x_grid.ravel(), y_grid.ravel()).reshape(x_grid.shape)

    return x_new, y_new, z_new


def loadDataFiles():
    # Open the ESA CCI SLA file.
    path = "data/data_NorthAtlantic/fol_1/*/S6A*STD*.nc"
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
    offset = 0

    # CASPBIAN SEA
    # latMin = 33.66 - offset
    # latMax = 53.53 + offset
    # lonMin = 39.16 - offset
    # lonMax = 61.78 + offset
    # AAlborg
    # latMin = 54 - offset
    # latMax = 58.06 + offset
    # lonMin = 7.72 - offset
    # lonMax = 12.78 + offset

    # Atlantic/North sea
    latMin = np.mod(0, 360)
    latMax = np.mod(65, 360)
    lonMin = 255
    lonMax = 360

    print(str(str(latMin) + "," + str(lonMin) + "," + str(latMax) + "," + str(lonMax).strip()))
    return latMin, latMax, lonMin, lonMax


def windowToAxes(window):
    return [window[0]['lon'], window[1]['lon'], window[1]['lat'], window[0]['lat']]


def main():
    # get variables and metadata
    print("Reading data")
    ncid = loadDataFiles()
    latList = read_variable_data(ncid, "latitude")
    lonList = read_variable_data(ncid, "longitude")
    altitudeDataList = read_variable_data(ncid, "mean_sea_surface_sol1")
    passNumbers = read_pass_number(ncid)
    time = read_variable_data(ncid, "time")

    timestamps = []
    base_date = datetime(2000, 1, 1)
    for ti in range(len(time)):
        for t in time[ti]:
            d = base_date + timedelta(seconds=t)
            d.strftime("%Y-%m-%d %H:%M:%S")
            timestamps.append(d)
    timestamps = np.array(timestamps)
    # Remove values that are outside the specified window

    ## Create A list of every repeated pass.
    # Orbits object has indices for which every data point repeats its pass over earth ground tracks.
    uniquePassNumbers = unique(passNumbers)
    orbits = []
    for orbitPass in uniquePassNumbers:
        indices = getIndicesOfItem(passNumbers, orbitPass)
        orbits.append(indices)

    print("Num passes: {}, Identical orbits: {}".format(len(orbits), len(orbits[0])))

    latMin, latMax, lonMin, lonMax = getLatLonBoundingBox()

    window = [{'lat': latMax, 'lon': lonMax}, {'lat': latMin, 'lon': lonMin}]

    # floats (left, right, bottom, top), optional
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
                        x[index].append(lon)
                        y[index].append(lat)
                        z[index].append(altitude)
                        h[index].append(altitude)

    n = 1
    print("Filtering some points")
    for i in range(len(x)):
        x[i] = x[i][::n]
        y[i] = y[i][::n]
        z[i] = z[i][::n]
        h[i] = h[i][::n]

        x[i] = np.nan_to_num(np.array(x[i]))
        x[i][x[i] == -inf] = 0
        x[i][x[i] == inf] = 0

        y[i] = np.nan_to_num(np.array(y[i]))
        y[i][y[i] == -inf] = 0
        y[i][y[i] == inf] = 0

        h[i] = np.nan_to_num(np.array(h[i]))
        h[i][h[i] == -inf] = 0
        h[i][h[i] == inf] = 0


    # _, _, imgInter = interpolateRBF(x[0], y[0], h[0])

    # plt.imshow(imgInter)
    # plt.show()
    # input("jhg")

    fig = plt.figure(figsize=(12, 7))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_global()



    def animate(i):
        # _, _, imgInter = interpolateRBF(x[i], y[i], h[i])
        ax.clear()
        # ax.imshow(imgInter, extent=window, transform=crs.PlateCarree(central_longitude=15))
        ax.add_feature(feat.LAND, zorder=100, edgecolor='k', alpha=0.3)
        ax.coastlines()
        plt.scatter(x=x[i], y=y[i], c=h[i], s=2, transform=crs.PlateCarree(central_longitude=15))
        ax.gridlines(draw_labels=False)
        fig.savefig("map_{}".format(i), orientation="landscape")
        print("Plotting i:", i)

    ani = animation.FuncAnimation(fig, animate,
                                  frames=len(orbits[0]), interval=1000, repeat=True)
    plt.show()

    ## AAlborg coordinates
    # lat, lon: 57.0469106,9.9167348


# surfacePlot(x, y, c)


if __name__ == '__main__':
    main()
