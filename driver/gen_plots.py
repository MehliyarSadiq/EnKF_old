from pylab import *
from numpy import *
from matplotlib.colors import BoundaryNorm

try:
    from matplotlib.toolkits.basemap import Basemap, shiftgrid
except ImportError:
    from mpl_toolkits.basemap import Basemap, shiftgrid


def plot_map(data, rlon=None, rlat=None, use_pcolor=0, **keywords):
    """ display  the data in a map
    keywords: dict which can include the following keywords  minv=0.0, maxv=0.0, dv=0.0,
    show_map=0, map_proj='cyl', show_lat=0, show_lon=0, lat_0=0.0, lon_0=0.0,
    minlat=-90.0, maxlat=90.0, minlon=0.0, maxlon=360.0, use_log=0, level=level
    """

    minv = 0.0
    if ('minv' in keywords):
        minv = keywords['minv']

    maxv = 0.0

    if ('maxv' in keywords):
        maxv = keywords['maxv']

    dv = 0.0
    if ('dv' in keywords):
        dv = keywords['dv']

    if ('dv' in keywords):
        dv = keywords['dv']

    if (maxv > minv):
        rlvl = arange(minv, maxv + dv, dv)
        ncolor = 256
        unorm = BoundaryNorm(rlvl, ncolor)

        # rlvl[0]  =-999.0
        # rlvl[size(rlvl)-1]=999.

    stitle = ""
    add_str = ""
    if ('title' in keywords):
        add_str = keywords['title']
        add_str = add_str.strip()
        stitle = stitle + ' ' + add_str

    if ('unit' in keywords):
        add_str = keywords['unit']
        add_str = add_str.strip()
        stitle = stitle + '(' + add_str + ')'

    cbar_vert = 1
    if ('cbar_vert' in keywords):
        cbar_vert = keywords['cbar_vert']

    if (cbar_vert == 1):
        orientation = 'vertical'
    else:
        orientation = 'horizontal'

    show_lat = 1
    if ('show_lat' in keywords):
        show_lat = keywords['show_lat']

    show_lon = 1
    if ('show_lon' in keywords):
        show_lon = keywords['show_lon']

    nlon, nlat = shape(data)

    vals = array(data)

    map_proj = 'cyl'
    if ('map_proj' in keywords):
        map_proj = keywords['map_proj']

    lat_0 = 0.0
    if ('lat_0' in keywords):
        lat_0 = keywords['lat_0']

    minlat = -90.0
    maxlat = 90.0

    if ('minlat' in keywords):
        minlat = keywords['minlat']

    if ('maxlat' in keywords):
        maxlat = keywords['maxlat']

    if (rlat is None):
        dlat = (maxlat - minlat) / nlat
        rlat = arange(minlat, maxlat, dlat)

    minlon = -180
    maxlon = 180.0

    if ('minlon' in keywords):
        minlon = keywords['minlon']

    if ('maxlon' in keywords):
        maxlon = keywords['maxlon']

    if (rlon is None):
        dlon = (maxlon - minlon) / nlon
        rlon = arange(minlon, maxlon, dlon)

    lon_0 = 0

    do_bdr = 0
    if ('do_bdr' in keywords):
        do_bdr = keywords['do_bdr']

    if ('lon_0' in keywords):
        lon_0 = keywords['lon_0']

        boundinglat = 45

    if ('boundinglat' in keywords):
        boundinglat = keywords['boundinglat']

    if (map_proj == 'npstere' or map_proj == 'spstere'):
        m = Basemap(projection=map_proj, lon_0=lon_0, boundinglat=boundinglat)
    elif (map_proj == 'ortho'):
        m = Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
    else:
        m = Basemap(
            llcrnrlon=minlon,
            llcrnrlat=minlat,
            urcrnrlon=maxlon,
            urcrnrlat=maxlat,
            projection=map_proj,
            lon_0=lon_0,
            lat_0=lat_0)

    if (rlon[-1] < rlon[0] + 360.0):
        rlon = resize(rlon, nlon + 1)
        rlon[-1] = rlon[0] + 360.0
        vals = squeeze(vals)
        vals = resize(vals, [nlon + 1, nlat])

    x, y = m(*meshgrid(rlon, rlat))

    cmap = cm.jet  # cm.Set1
    if ('cmap' in keywords):
        cmap = keywords['cmap']

    if (maxv > minv):
        if (use_pcolor == 1):
            # cs0=m.pcolor(x, y, transpose(vals),  edgecolors='none', norm=unorm,  cmap=cmap)
            cs0 = m.pcolormesh(
                x,
                y,
                transpose(vals),
                norm=unorm,
                edgecolors='none',
                cmap=cmap)

        else:
            cs0 = m.contourf(x, y, transpose(vals), rlvl)
    else:
        if (use_pcolor == 1):
            # cs0=m.pcolor(x, y, transpose(vals), edgecolors='none', cmap=cmap)
            cs0 = m.pcolormesh(x, y, transpose(vals), cmap=cmap)

        else:

            cs0 = m.contourf(x, y, transpose(vals))

    # info(m.drawcoastlines)
    m.drawcoastlines(color='k', linewidth=1.5)
    # m.drawcountries(color=white)
    m.drawmapboundary()
    if (show_lat == 1):
        m.drawparallels(
            arange(
                minlat,
                maxlat +
                30.0,
                30.),
            labels=[
                1,
                0,
                0,
                0])

    if (show_lon == 1):
        m.drawmeridians(arange(minlon, maxlon + 60, 60.), labels=[0, 0, 0, 1])
        title(stitle)

    show_colorbar = 1
    if ('cb' in keywords):
        show_colorbar = keywords['cb']
    if (show_colorbar == 1):
        colorbar(orientation=orientation, extend='both')
        # sel_div=arange(minv, maxv+dv, 2.5*dv)
        # yticks(sel_div, sel_div)

    if (do_bdr == 1):
        lvl = arange(max(vals.flat))
        print(shape(x), shape(y), shape(vals))
        # cs2=m.contour(x[0:-1,:], y[0:-1,:], transpose(vals), lvl, colors='k', linewidth=0.5)

    return m


def plot_track(rlon, rlat=None, m=None, **keywords):
    show_lat = 1
    if ('show_lat' in keywords):
        show_lat = keywords['show_lat']

    show_lon = 1
    if ('show_lon' in keywords):
        show_lon = keywords['show_lon']

    if (m is None):
        map_proj = 'cyl'
        if ('map_proj' in keywords):
            map_proj = keywords['map_proj']

        lat_0 = 0.0
        if ('lat_0' in keywords):
            lat_0 = keywords['lat_0']

        minlat = -90.0
        maxlat = 90.0

        if ('minlat' in keywords):
            minlat = keywords['minlat']

        if ('maxlat' in keywords):
            maxlat = keywords['maxlat']

        minlon = -180
        maxlon = 180.0

        if ('minlon' in keywords):
            minlon = keywords['minlon']

        if ('maxlon' in keywords):
            maxlon = keywords['maxlon']

        lon_0 = 0

        if ('lon_0' in keywords):
            lon_0 = keywords['lon_0']

        boundinglat = 45

        if ('boundinglat' in keywords):
            boundinglat = keywords['boundinglat']

        if (map_proj == 'npstere' or map_proj == 'spstere'):
            m = Basemap(
                projection=map_proj,
                lon_0=lon_0,
                boundinglat=boundinglat)
        elif (map_proj == 'ortho'):
            m = Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        else:
            m = Basemap(
                llcrnrlon=minlon,
                llcrnrlat=minlat,
                urcrnrlon=maxlon,
                urcrnrlat=maxlat,
                projection=map_proj,
                lon_0=lon_0,
                lat_0=lat_0)
        if (show_lat == 1):
            m.drawparallels(
                arange(
                    minlat,
                    maxlat +
                    30.0,
                    30.),
                labels=[
                    1,
                    0,
                    0,
                    0])
        if (show_lon == 1):
            m.drawmeridians(
                arange(
                    minlon,
                    maxlon + 60,
                    60.),
                labels=[
                    0,
                    0,
                    0,
                    1])
    sgn = 'o'
    if ('sgn' in keywords):
        sgn = keywords['sgn']

    if ('color' in keywords):
        print('I am here', size(rlon))

        lcolor = keywords['color']
        print(max(lcolor), min(lcolor))

        if ('minv' in keywords):
            vmin = keywords['minv']
        else:
            vmin = min(lcolor)

        if ('maxv' in keywords):
            vmax = keywords['maxv']
        else:
            vmax = max(lcolor)

        cx = cm.jet
        cx.set_over('r')
        cx.set_under('w')
        x, y = m(rlon, rlat)

        if (len(rlon) < 100):
            scatter(
                x,
                y,
                marker=sgn,
                c=lcolor,
                s=30.0,
                color='w',
                vmin=vmin,
                vmax=vmax,
                cmap=cx)  # , markersize=11)
        else:
            scatter(x, y, marker=sgn, c=lcolor, s=26, edgecolor='w',
                    vmin=vmin, vmax=vmax, cmap=cx)  # , markersize=6)

    cbar_vert = 1
    if ('cbar_vert' in keywords):
        cbar_vert = keywords['cbar_vert']

    if (cbar_vert == 1):
        orientation = 'vertical'
    else:
        orientation = 'horizontal'

    show_colorbar = 0
    if ('cb' in keywords):
        show_colorbar = keywords['cb']
    if (show_colorbar == 1):
        colorbar(orientation=orientation)

    # colorbar(cmap=cx)

    m.drawcoastlines(color='k', linewidth=1.0)
    m.drawmapboundary()
    # m.drawcountries(color='k', linewidth=0.5)

    if ('title' in keywords):
        stitle = keywords['title']
        title(stitle)

    return m


def add_text(rlon, rlat, txt, m=None, **keywords):
    show_lat = 1
    if ('show_lat' in keywords):
        show_lat = keywords['show_lat']
    rlon = array(rlon)
    rlat = array(rlat)

    show_lon = 1
    if ('show_lon' in keywords):
        show_lon = keywords['show_lon']

    if (m is None):
        map_proj = 'cyl'
        if ('map_proj' in keywords):
            map_proj = keywords['map_proj']

        lat_0 = 0.0
        if ('lat_0' in keywords):
            lat_0 = keywords['lat_0']

        minlat = -90.0
        maxlat = 90.0

        if ('minlat' in keywords):
            minlat = keywords['minlat']

        if ('maxlat' in keywords):
            maxlat = keywords['maxlat']

        minlon = -180
        maxlon = 180.0

        if ('minlon' in keywords):
            minlon = keywords['minlon']

        if ('maxlon' in keywords):
            maxlon = keywords['maxlon']

        lon_0 = 0

        if ('lon_0' in keywords):
            lon_0 = keywords['lon_0']

        boundinglat = 45

        if ('boundinglat' in keywords):
            boundinglat = keywords['boundinglat']

        if (map_proj == 'npstere' or map_proj == 'spstere'):
            m = Basemap(
                projection=map_proj,
                lon_0=lon_0,
                boundinglat=boundinglat)
        elif (map_proj == 'ortho'):
            m = Basemap(projection=map_proj, lon_0=lon_0, lat_0=lat_0)
        else:
            m = Basemap(
                llcrnrlon=minlon,
                llcrnrlat=minlat,
                urcrnrlon=maxlon,
                urcrnrlat=maxlat,
                projection=map_proj,
                lon_0=lon_0,
                lat_0=lat_0)
        if (show_lat == 1):
            m.drawparallels(
                arange(
                    minlat,
                    maxlat +
                    30.0,
                    30.),
                labels=[
                    1,
                    0,
                    0,
                    0])
        if (show_lon == 1):
            m.drawmeridians(
                arange(
                    minlon,
                    maxlon + 60,
                    60.),
                labels=[
                    0,
                    0,
                    0,
                    1])
    x, y = m(rlon, rlat)

    for i in range(size(rlon)):
        if (i < 15):
            text(
                x[i],
                y[i],
                txt[i],
                fontsize=10,
                color="w",
                horizontalalignment='center')
        else:
            text(
                x[i],
                y[i],
                txt[i],
                fontsize=10,
                color="k",
                horizontalalignment='center')

    if ('title' in keywords):
        stitle = keywords['title']
        title(stitle)

    return m


if (__name__ == "__main__"):

    rlon = arange(-180, 180, 10)
    rlat = arange(-90, 90, 10)
    lon_m, lat_m = meshgrid(rlon, rlat)
    data = sin(lon_m * pi / 180.0)**2 + cos(lat_m * pi / 180.0)**2
    data = transpose(data)
    subplot(2, 1, 1)
    plot_map(data, rlon, rlat, use_pcolor=1, maxv=1.0, minv=-1.0, dv=0.2)
    show()
