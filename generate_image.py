import glob
import random
from datetime import datetime, timedelta

import cartopy.crs as crs
import cartopy.feature as feat
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from matplotlib import animation
from scipy.interpolate import Rbf

earthRadius = r = 6371
selectedDay = 27
selectedMonth = 11


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
    row = []
    for x in range(len(ncid)):
        row.append(ncid[x].variables[varname][:])
    return row


def read_parent_variable_data(ncid, varname):
    row = []
    for x in range(len(ncid)):
        row.append(ncid[x].parent.group[varname][:])
    return row


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
        ncid.append(nc.Dataset(file).groups["data_20"].groups["ku"])
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
    meanSeaSurfaceList = read_variable_data(ncid, "mean_sea_surface_sol1")
    oceanAltimeterRange = read_variable_data(ncid, "range_ocean")
    time = read_variable_data(ncid, "time")

    timestamps = []
    latArr = []
    lonArr = []
    meanSeaSurfaceArr = []
    oceanAltimeterRangeArr = []

    base_date = datetime(2000, 1, 1)

    # Put all our data into single arrays
    print("Putting all data from different orbits into 1D arrays")
    for ti in range(len(time)):
        latArrTmp = latList[ti]
        lonArrTmp = lonList[ti]
        altArrTmp = meanSeaSurfaceList[ti]
        timeArrTmp = time[ti]
        oceanAltimeterRangeArrTmp = oceanAltimeterRange[ti]
        for x in range(len(latArrTmp)):
            latArr.append(latArrTmp[x])
            lonArr.append(lonArrTmp[x])
            meanSeaSurfaceArr.append(altArrTmp[x])
            oceanAltimeterRangeArr.append(oceanAltimeterRangeArrTmp[x])
            t = timeArrTmp[x]
            d = base_date + timedelta(seconds=t)
            d.strftime("%Y-%m-%d %H:%M:%S")
            timestamps.append(d)

    latArr = np.array(latArr)
    lonArr = np.array(lonArr)
    meanSeaSurfaceArr = np.array(meanSeaSurfaceArr)
    oceanAltimeterRangeArr = np.array(oceanAltimeterRangeArr)
    timestamps = np.array(timestamps)

    latMin, latMax, lonMin, lonMax = getLatLonBoundingBox()
    window = [{'lat': latMax, 'lon': 15}, {'lat': latMin, 'lon': -90}]

    latitudes = []
    longitudes = []
    meanSeaSurfaces = []
    times = []
    oceanRange = []
    oceanRangeSeaSurfaceChange = []

    print("Sorting points, exluding points outside the bounding box")
    for i in range(len(latArr)):
        if latMin < latArr[i] < latMax:
            if lonMin < lonArr[i] < lonMax:
                latitudes.append(latArr[i])
                longitudes.append(lonArr[i])
                times.append(timestamps[i])
                oceanRange.append(oceanAltimeterRangeArr[i])

                if np.isnan(oceanAltimeterRangeArr[i]) or np.isinf(oceanAltimeterRangeArr[i]) or \
                        np.isnan(meanSeaSurfaceArr[i]) or np.isinf(meanSeaSurfaceArr[i]):
                    oceanRangeSeaSurfaceChange.append(0)
                    meanSeaSurfaces.append(1)
                    print("Nan value")
                else:
                    oceanRangeSeaSurfaceChange.append(meanSeaSurfaceArr[i] - oceanAltimeterRangeArr[i])
                    meanSeaSurfaces.append(meanSeaSurfaceArr[i])

    latitudes = np.array(latitudes)
    longitudes = np.array(longitudes)
    meanSeaSurfaces = np.array(meanSeaSurfaces)
    oceanRange = np.array(oceanRange)
    oceanRangeSeaSurfaceChange = np.array(oceanRangeSeaSurfaceChange)

    times = np.array(times)

    fig = plt.figure(figsize=(12, 7))
    ax = plt.axes(projection=crs.PlateCarree())

    def animate(i):
        global selectedDay, selectedMonth

        if selectedMonth % 2 != 0 and selectedDay % 31 == 0:
            selectedMonth += 1
            selectedDay = 1

        if selectedMonth % 2 == 0 and selectedDay % 32 == 0:
            selectedMonth = 1
            selectedDay = 1

        print("mm:dd, ", selectedMonth, selectedDay)
        x = []
        y = []
        h = []
        for j in range(len(latitudes)):
            day = times[j].day
            month = times[j].month
            if day == selectedDay and month == selectedMonth:
                y.append(latitudes[j])
                x.append(longitudes[j])
                h.append(oceanRangeSeaSurfaceChange[j])
        x = np.array(x)
        y = np.array(y)
        h = np.array(h)

        # Use only every nth point
        n = 1000
        x = x[::n]
        y = y[::n]
        h = h[::n]

        print("plotting {}, points: {}".format(i, len(x)))

        _, _, imgInter = interpolateRBF(x, y, h)
        imgInter = np.rot90(imgInter, 2)

        ax.clear()
        ax.set_extent(windowToAxes(window))
        ax.add_feature(feat.LAND, zorder=100, edgecolor='k', alpha=0.3)
        ax.coastlines()
        plt.scatter(x=x, y=y, c="r", s=2, transform=crs.PlateCarree(central_longitude=15))
        ax.imshow(imgInter, extent=windowToAxes(window), transform=crs.PlateCarree(central_longitude=0))
        ax.gridlines(draw_labels=False)
        fig.savefig("output_images/high_res_map_{}".format(i), orientation="landscape")
        selectedDay += 1

    ani = animation.FuncAnimation(fig, animate, frames=100, interval=1000, repeat=True)
    plt.show()


if __name__ == '__main__':
    main()
