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
selectedDay = [9, 10, 11, 12, 13, 14]
selectedMonth = 11
minMaxList = []


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
        row.append(ncid[x].variables[varname])
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


def loadDataFiles(dataType):
    #  the ESA CCI SLA file.
    if "Jason" in dataType:
        path = "data/JSON-3/*/*.nc"
        files = glob.glob(path)
        ncid = []
        for i, file in enumerate(files):
            if i == 25:
                break
            ncid.append(nc.Dataset(file))
        return ncid

    else:
        path = "data/data_CaspianSea/*/S6A*STD*.nc"
        files = glob.glob(path)
        ncid = []
        for file in files:
            if "fol_5" in file:
                break
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
    # latMin = np.mod(0, 360)
    # latMax = np.mod(65, 360)
    # lonMin = 255
    # lonMax = 360

    # Caspian lake
    latMin = np.mod(36, 360)
    latMax = np.mod(48, 360)
    lonMin = 46
    lonMax = 55

    print(str(str(latMin) + "," + str(lonMin) + "," + str(latMax) + "," + str(lonMax).strip()))
    return latMin, latMax, lonMin, lonMax


def windowToAxes(window):
    return [window[0]['lon'], window[1]['lon'], window[1]['lat'], window[0]['lat']]


def main():
    # get variables and metadata
    print("Reading data")
    ncid = loadDataFiles("sentinel")
    latList = read_variable_data(ncid, "lat")
    lonList = read_variable_data(ncid, "lon")
    #seaSurfaceQuality = read_variable_data(ncid, "mean_sea_surface_sol1_qual")
    #seaSurfaceAccuracy = read_variable_data(ncid, "mean_sea_surface_sol1_acc")

    meanSeaSurfaceList = read_variable_data(ncid, "mean_sea_surface")
    #anomalyList = read_variable_data(ncid, "ssha")
    time = read_variable_data(ncid, "time")

    del ncid
    print(len(latList))
    print(len(latList[0]))

    latArr = [0] * (len(time) * len(latList[0] * 2))
    lonArr = [0] * (len(time) * len(latList[0] * 2))
    meanSeaSurfaceArr = [0] * (len(time) * len(latList[0] * 2))
    oceanAltimeterRangeArr = [0] * (len(time) * len(latList[0] * 2))
    anomalies = [0] * (len(time) * len(latList[0] * 2))
    base_date = datetime(2000, 1, 1)
    timestamps = [base_date] * (len(time) * len(latList[0] * 2))

    latMin, latMax, lonMin, lonMax = getLatLonBoundingBox()
    window = [{'lat': latMax, 'lon': lonMax}, {'lat': latMin, 'lon': lonMin}]

    # Put all our data into single arrays
    print("Filtering data from different orbits and inserting into 1D arrays")
    index = 0
    for ti in range(len(time)):
        latArrTmp = latList[ti]
        lonArrTmp = lonList[ti]
        altArrTmp = meanSeaSurfaceList[ti]
        timeArrTmp = time[ti]
        #anomalyListTmp = anomalyList[ti]
        #qualityListTmp = seaSurfaceQuality[ti]
        #seaSurfaceAccuracyTmp = seaSurfaceAccuracy[ti]
        for x in range(len(latArrTmp)):
            #if qualityListTmp[x] == 1 or seaSurfaceAccuracyTmp[x] > 0.05:
            #    continue
            if latMin < latArrTmp[x] < latMax:
                if lonMin < lonArrTmp[x] < lonMax:
                    latArr[index] = latArrTmp[x]
                    lonArr[index] = lonArrTmp[x]
                    meanSeaSurfaceArr[index] = altArrTmp[x]
                    #anomalies[index] = anomalyListTmp[x]
                    t = timeArrTmp[x]
                    d = base_date + timedelta(seconds=t)
                    d.strftime("%Y-%m-%d %H:%M:%S")
                    timestamps[index] = d
                    index += 1

    print("Sorting out Infs and Nans")
    latitudes = np.array(latArr)
    longitudes = np.array(lonArr)
    meanSeaSurfaceArr = np.nan_to_num(np.array(meanSeaSurfaceArr), nan=0.0, posinf=0.0, neginf=0.0)
    oceanAltimeterRangeArr = np.array(oceanAltimeterRangeArr)
    times = np.array(timestamps)
    anomalies = np.nan_to_num(np.array(anomalies), nan=0.0, posinf=0.0, neginf=0.0)

    fig = plt.figure(figsize=(12, 7))
    ax = plt.axes(projection=crs.PlateCarree())
    print("Animating data. There are in total {} points within our region".format(len(latitudes)))

    def animate(i):
        global selectedDay, selectedMonth, minMaxList
        print(minMaxList)

        """
                if selectedMonth % 2 != 0 and selectedDay % 31 == 0:
            selectedMonth += 1
            selectedDay = 1

        if selectedMonth % 2 == 0 and selectedDay % 32 == 0:
            selectedMonth = 1
            selectedDay = 1
        """
        if selectedDay[len(selectedDay) -1 ] == 31:
            return

        print("mm:dd, ", selectedMonth, selectedDay)
        x = []
        y = []
        h = []
        for j in range(len(latitudes)):
            day = times[j].day
            month = times[j].month
            if 0 < selectedDay.count(day) and month == selectedMonth:
                y.append(latitudes[j])
                x.append(longitudes[j])
                h.append(meanSeaSurfaceArr[j])
        x = np.array(x)
        y = np.array(y)
        h = np.array(h)

        # Use only every nth point

        n = 50

        x = x[::n]
        y = y[::n]
        h = h[::n]
        print("plotting {}, points: {}".format(i, len(x)))
        if len(x) < 2:
            selectedDay += 1
            return

        minMaxList.append((np.amin(h), np.amax(h)))
        print(np.amin(h), np.amax(h))
        _, _, imgInter = interpolateRBF(x, y, h)
        imgInter = np.rot90(imgInter, 2)

        ax.clear()
        # ax.set_global()
        ax.set_extent(windowToAxes(window))
        ax.add_feature(feat.LAND, zorder=100, edgecolor='k', alpha=1)
        ax.coastlines()
        plt.scatter(x=x, y=y, c="r", s=2, transform=crs.PlateCarree(central_longitude=0))
        ax.imshow(imgInter, vmin=-43, vmax=-29, cmap="jet", extent=windowToAxes(window),
                  transform=crs.PlateCarree(central_longitude=0))
        ax.gridlines(draw_labels=False)
        fig.savefig("output_images/ssha{}".format(i), orientation="landscape")
        # selectedDay += 1

        # Shift list one day.
        selectedDay[0] = selectedDay[1]
        selectedDay[len(selectedDay) - 1] = selectedDay[len(selectedDay) - 1] + 1


    ani = animation.FuncAnimation(fig, animate, frames=100, interval=1000, repeat=False)
    plt.show()


if __name__ == '__main__':
    main()
