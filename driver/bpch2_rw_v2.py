""" class for reading and writing binary punch files
    history:
      0.5 by lf 2007.10.10

"""
from pylab import *

try:
    from matplotlib.toolkits.basemap import Basemap, shiftgrid
except ImportError:
    from mpl_toolkits.basemap import Basemap, shiftgrid

from numpy import *
import bpch2_rw_py
import time_module as tm
import geo_constant as gc
import gp_axis as gax
import geos_chem_def as gcdf

bpch2_fill_val = -999
bpch2_fill_str_val = "???"
bpch2_fill_flt_val = -999.0


class bpgrid:
    """ the class for the grid used to store GEO-CHEM binary punch data"""

    def __init__(self, ix, jx, lx,
                 ifirst=1, jfirst=1, lfirst=1,
                 halfpolar=0, centre180=1,
                 lonres=None, latres=None):
        """ Initialize the class
        Arguments:
        ix, jx, lx: 		integer -- the sizes of lon, lat, z
        ifirst, jfirst, lfirst: the start point for the data to store
        halfpolar, centre180: 	0 or 1 the type of the lon
        lonres, latres:    	the resolution of the lat and lon
        In the future, we will include the details on the vertical grid as well.
        """
        self.ix = ix
        self.jx = jx
        self.lx = lx
        self.ifirst = ifirst
        self.jfirst = jfirst
        self.lfirst = lfirst
        self.halfpolar = halfpolar
        self.centre180 = centre180
        if (lonres is None):
            self.lonres = 360.0 / (ix)
        else:
            self.lonres = lonres

        if (latres is None):
            self.latres = 180.0 / (jx)
        else:
            self.latres = latres
        self.lats = self.get_lat()
        self.rons = self.get_lon()
        self.zs = self.get_z()

    def get_lat(self):
        lat = arange(0.0, self.jx)
        lat = self.latres * lat - 90.0
        return lat

    def get_lon(self):
        lon = arange(0.0, self.ix)
        if (self.centre180 == 0):
            lon = self.lonres * lon
        else:
            lon2 = self.lonres * lon
            halfid = self.ix / 2
            if (mod(self.ix, 2) > 0):
                lon[0:halfid] = lon2[halfid + 1:] - 360
                lon[halfid:] = lon2[0:halfid + 1]

            else:
                lon[0:halfid] = lon2[halfid:] - 360
                lon[halfid:] = lon2[0:halfid]
        return lon

    def get_z(self):
        level = arange(self.lx)
        return level


class bpch2_data:
    def __init__(self, ntracer, category, unit, bpgrid, data, **keywords):
        """ initialize
        Arguments:
        ntracer:	integer, the order number of the tracer in GEOS-CHEM
        category:	string, the name of the tracer
        unit:		string, the name of the unit
        bpgrid:		class bpgrid, the grid of the data 	def get_lon(self):
        lon=arange(0.0, self.ix)


        the keyword: 	additional information given in form of  tau0=a, tau1=b, modelname='md', reserved='a'
        """
        self.ntracer = ntracer  # the tracer no in model
        self.grid = bpgrid
        self.data = data
        self.category = category
        self.unit = unit

        self.attr = {}
        if (len(keywords) > 0):
            self.attr.update(keywords)

    def write_to_bpch2_file(self, funit):
        """ write data into one open bpch2 file
        Arguments:
        funit:	integer- the unit number of the bpch2 file

        """

        ix = self.grid.ix
        jx = self.grid.jx
        lx = self.grid.lx

        lonres = self.grid.lonres
        latres = self.grid.latres
        ifirst = self.grid.ifirst
        jfirst = self.grid.jfirst
        lfirst = self.grid.lfirst
        halfpolar = self.grid.halfpolar
        centre180 = self.grid.centre180

        ntracer = self.ntracer
        data = self.data
        category = self.category
        unit = self.unit
        tau0 = 0.0
        tau1 = 0.0
        modelname = " "
        reserved = " "
        modelname = " "

        if ('modelname' in self.attr):
            modelname = self.attr['modelname']

        if ('tau0' in self.attr):
            tau0 = self.attr['tau0']

        if ('tau1' in self.attr):
            tau1 = self.attr['tau1']

        if ('reserved' in self.attr):
            reserved = self.attr['reserved']

        reserved = reserved.strip()

        if (len(reserved) == 0):
            reserved = "000000"
        traunit = self.attr['traunit']
        # print unit, traunit

        stat = bpch2_rw_py.write_bpch2_data(
            funit,
            modelname,
            category,
            reserved,
            lonres,
            latres,
            halfpolar,
            centre180,
            ntracer,
            unit,
            tau0,
            tau1,
            ifirst,
            jfirst,
            lfirst,
            data)
        return stat

    def regrid_data(self, new_lon, new_lat):
        """ write data into one open bpch2 file
        Arguments:
        funit:	integer- the unit number of the bpch2 file

        """

        ix = self.grid.ix
        jx = self.grid.jx
        lx = self.grid.lx

        lonres = self.grid.lonres
        latres = self.grid.latres
        ifirst = self.grid.ifirst
        jfirst = self.grid.jfirst
        lfirst = self.grid.lfirst
        halfpolar = self.grid.halfpolar
        centre180 = self.grid.centre180

        ntracer = self.ntracer
        data = self.data
        category = self.category

        lat = self.grid.get_lat()
        lon = self.grid.get_lon()
        ax_lat = gax.gp_axis('lat', lat)
        ax_lon = gax.gp_axis('lat', lon)
        lonp1, lonp2, lonw = ax_lon.getwgt(new_lon)
        latp1, latp2, latw = ax_lat.getwgt(new_lat)
        new_ix = size(new_lon)
        new_jx = size(new_lat)
        new_lx = lx
        # new_data=zeros([new_ix, new_jx, new_lx], float)
        new_data1 = lonw[:, newaxis, newaxis] * data[lonp1, :, :] + \
            (1.0 - lonw[:, newaxis, newaxis]) * data[lonp2, :, :]
        new_data = latw[newaxis, :, newaxis] * new_data1[:, latp1, :] + \
            (1.0 - latw[newaxis, :, newaxis]) * new_data1[:, latp2, :]
        self.data = array(new_data)
        new_grid = bpgrid(new_ix, new_jx, new_lx)

        self.grid = new_grid

    def write_to_netcdf(self, flnm):
        """ write the data into one netcdf """
        ntau = 2
        times = arange(ntau)
        tau0 = 0.0
        tau1 = 0.0
        if ('tau0' in self.attr):
            tau0 = 1.0 * self.attr['tau0']

        if ('tau1' in self.attr):
            tau1 = 1.0 * self.attr['tau1']
        rlat = self.grid.get_lat()
        rlon = self.grid.get_lat()
        rz = self.grid.get_z()
        rtime = arange([tau0, tau1])

    # dimensions as channel number, segment in one array for each tangent point
    dimNames = ['lon', 'lat', 'level', 'time']
    dimTypes = ['f', 'f', 'f', 'f']
    dimVars = [rlon, rlat, rz, rtime]

    nf.netCDF_def_dims(flnm, dimNames, dimTypes, dimVars)
    varType = 'f'
    dimNames = ['lon', 'lat', 'level']
    varName = self.category
    varData = 1.0 * array(self.data)

    nf.netCDF_var_write(flnm, dimNames, varName, varType, varData)

    def set_attr(self, attr_name, value):
        """ set one attribute to the data set """
        self.attr.update({attr_name: value})

    def get_attr(self, attr_names):
        """ get attributes of the data set  """
        if (attr_names is None):
            return self.attr

        at = list()
        ats = list()
        if (isinstance(attr_names, type(at))):
            at = attr_names
        else:
            at.append(attr_names)

        for attr_name in at:
            if (attr_name in self.attr):
                ats.append(self.attr[attr_name])
            else:
                print('no attribute ' + attr_name + 'in data ' + self.category)
                return None

        return ats

    def update_data(self, data):
        if (self.data is not None):
            del self.data
        self.data = data

    def display(self, dirc, **keywords):
        """ display  the data
        keywords: dict which can include the following keywords  minv=0.0, maxv=0.0, dv=0.0,
        show_map=0, map_proj='cyl', show_lat=0, show_lon=0, lat_0=0.0, lon_0=0.0,
        minlat=-90.0, maxlat=90.0, minlon=0.0, maxlon=360.0, use_log=0, level=level
        """

        rlat = self.grid.get_lat()
        rlon = self.grid.get_lon()

        levels = self.grid.get_z()

        minv = 0.0
        if ('minv' in keywords):
            minv = keywords['minv']

        maxv = 0.0

        if ('maxv' in keywords):
            maxv = keywords['maxv']

        dv = 0.0
        if ('dv' in keywords):
            dv = keywords['dv']

        show_map = 0

        if ('dv' in keywords):
            dv = keywords['dv']

        factor = 1.0
        if ('factor' in keywords):
            factor = keywords['factor']

        if (maxv > minv):
            rlvl = arange(minv, maxv + dv, dv)
            rlvl[0] = -999.0
            rlvl[size(rlvl) - 1] = 999.

        show_map = 0
        if ('show_map' in keywords):
            show_map = keywords['show_map']

        stitle = ""
        if ('name' in self.attr):
            add_str = self.attr['name']
        else:
            add_str = self.category
        add_str = add_str.strip()
        stitle = stitle + ' ' + add_str

        if ('traunit' in self.attr):
            add_str = self.attr['traunit']
        else:

            add_str = self.unit
        add_str = self.unit
        add_str = add_str.strip()
        cbar_vert = 1
        if ('cbar_vert' in keywords):
            cbar_vert = keyowrds['cbar_vert']

        if (cbar_vert == 1):
            orientation = 'vertical'
        else:
            orientation = 'horizontal'

        stitle = stitle + ' (' + add_str + ')'
        if ('tau0' in self.attr):
            tau0 = self.attr['tau0']
            tau0 = 3600. * tau0  # convert to seconds
            utc = tm.tai85_to_utc(tau0)
            stitle = stitle + ' ' + utc

        if (show_map == 1):
            level = 0
            if ('level' in keywords):
                level = keywords['level']

            show_lat = 1
            if ('show_lat' in keywords):
                show_lat = keywords['show_lat']

            show_lon = 1
            if ('show_lon' in keywords):
                show_lon = keywords['show_lon']

            vals = self.data[:, :, level]
            vals = array(vals)
            vals = factor * vals

            map_proj = 'cyl'
            if ('map_proj' in keywords):
                map_proj = keywords['map_proj']
            lat_0 = 0.0
            if ('lat_0' in keywords):
                lat_0 = keywords['lat_0']

            minlat = -90.0
            if ('minlat' in keywords):
                minlat = keywords['minlat']

            maxlat = 90.0

            if ('maxlat' in keywords):
                maxlat = keywords['maxlat']

            if (self.grid.centre180 == 0):
                minlon = 0.0
                maxlon = 360.0
                lon_0 = 180
            else:
                minlon = -180
                maxlon = 180.0
                lon_0 = 0

            if ('minlon' in keywords):
                minlon = keywords['minlon']

            if ('maxlon' in keywords):
                maxlon = keywords['maxlon']

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

            if (rlon[-1] < rlon[0] + 360.0):
                rlon = resize(rlon, self.grid.ix + 1)
                rlon[-1] = rlon[0] + 360.0
                vals = squeeze(vals)
                vals = resize(vals, [self.grid.ix + 1, self.grid.jx])

            x, y = m(*meshgrid(rlon, rlat))

            if (maxv > minv):
                cs0 = m.contourf(x, y, transpose(vals), rlvl)
            else:
                cs0 = m.contourf(x, y, transpose(vals))

            m.drawcoastlines()
            m.drawcountries()
            m.drawmapboundary()
            if (show_lat == 1):
                print(minlat, maxlat)
                print(minlon, maxlon)
                m.drawparallels(
                    arange(
                        minlat,
                        maxlat + 30.0,
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
            title(stitle)
            colorbar(orientation=orientation)
        elif (dirc == 0):
            level = 0
            if ('level' in keywords):
                lvl = keywords['level']

            vals = self.data[:, :, level]
            vals = factor * vals
            print(shape(vals), shape(rlon), shape(rlat))

            if (maxv > minv):
                contourf(rlon, rlat, transpose(vals), rlvl)
            else:
                contourf(rlon, rlat, transpose(vals))
            # restart_new
            title(stitle)
            xlabel('Lon')
            ylabel('Lat')

            colorbar()
        elif (dirc == 1):
            vals = average(self.data, axis=0)
            use_log = 0
            if ('use_log' in keywords):
                use_log = keywords['use_log']
            ax = subplot(2, 1, 1)
            if (use_log == 1):
                ax.set_yscale('log')
            levels = self.get_level()
            if ('level' in keywords):
                levels = keywords['levels']

            if (maxv > minv):
                contourf(rlon, rlat, transpose(vals), rlvl)
            else:
                contourf(rlat, levels, vals)

            xlabel('Lat')
            ylabel('Level')

            title(stitle)
            colorbar()

        show()
        return 0


class tracer_info:
    def __init__(self, flnm):
        self.name = list()
        self.fullname = list()
        self.tracer = list()
        self.scale = list()
        self.unit = list()
        self.c = list()
        self.molew = list()
        file_found = False
        if (flnm.strip() != ""):
            try:
                fl = open(flnm.strip(), "r")
                file_found = True
            except IOError:
                print('Tracer info file ' + flnm + ' not found')
                file_found = False

        if (file_found):
            lines = fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0] != '#'):

                    sname, sfullname, smw, sc, stra, sscal, sunit = \
                        line[0:8], line[9:39], line[39:49], line[49:53], line[53:62], line[62:72], line[73:]
                    # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.name.append(sname.strip())
                    self.fullname.append(sfullname.strip())
                    self.molew.append(float(smw))
                    self.c.append(float(sc))
                    self.scale.append(float(sc))
                    self.tracer.append(int(stra))
                    self.unit.append(sunit.strip())

            self.tracer = array(self.tracer)

    def get_tracer_info(self, tracer_in):
        if (len(self.tracer) == 0):
            return bpch2_fill_str_val, bpch2_fill_str_val, bpch2_fill_flt_val, bpch2_fill_flt_val, bpch2_fill_str_val

        idx = where(tracer_in == self.tracer)
        # idx=compress(idx)
        if (size(idx) >= 1):
            name = self.name[idx[0]]
            fullname = self.fullname[idx[0]]
            scale = self.scale[idx[0]]
            unit = self.unit[idx[0]]
            c = self.c[idx[0]]
            molew = self.molew[idx[0]]
            return name, fullname, molew, scale, unit
        else:
            return bpch2_fill_str_val, bpch2_fill_str_val, bpch2_fill_flt_val, bpch2_fill_flt_val, bpch2_fill_str_val

    def load_tracer_info(self, flnm):
        file_found = False
        if (flnm.strip() != ""):
            try:
                fl = open(flnm.strip(), "r")
                file_found = True
            except IOError:
                print('No tracer info file: ', flnm.strip())
                file_found = False

        if (file_found):
            lines = fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0] != '#'):

                    sname, sfullname, smw, sc, stra, sscal, sunit = \
                        line[0:8], line[9:39], line[39:49], line[49:53], line[53:62], line[62:72], line[73:]
                    # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.name.append(sname.strip())
                    self.fullname.append(sfullname.strip())
                    self.molew.append(float(smw))
                    self.c.append(float(sc))
                    self.scale.append(float(sc))
                    self.tracer.append(int(stra))
                    self.unit.append(sunit.strip())

            self.tracer = array(self.tracer)


class diag_info:
    def __init__(self, flnm):
        self.category = list()
        self.comment = list()
        self.offset = list()
        file_found = False
        if (flnm.strip() != ""):
            try:
                fl = open(flnm.strip(), "r")
                file_found = True
            except IOError:
                print('Diag info file ' + flnm + ' not found')
                file_found = False
        if (file_found):
            lines = fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0] != '#'):

                    soffset, scategory, scomment = line[0:8], line[9:49], line[49:]
                    # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.offset.append(int(soffset))
                    self.category.append(scategory.strip())
                    self.comment.append(scomment.strip())

    def get_offset(self, category):
        nid = len(self.category)
        found = 0
        scate = category.strip()
        for id in range(nid):
            if (scate in self.category[id]):
                found = 1
                offset = self.offset[id]

        if (found):
            return offset
        else:
            return bpch2_fill_flt_val

    def load_diag_info(self, flnm):
        file_found = False
        if (flnm.strip() != ""):
            try:
                fl = open(flnm.strip(), "r")
                file_found = True
            except IOError:
                print('no diag info file: ', flnm.strip())
                file_found = False
        if (file_found):
            lines = fl.readlines()
            fl.close()
            for line in lines:
                # line=line.strip()
                if (line[0] != '#'):

                    soffset, scategory, scomment = line[0:8], line[9:49], line[49:]
                    # print  sname, sfullname, smw, sc, stra, sscal, sunit
                    self.offset.append(int(soffset))
                    self.category.append(scategory.strip())
                    self.comment.append(scomment.strip())


class bpch2_file_rw:
    """ the class for read and write data  """

    def __init__(self, flnm, mode, maxtracer=300,
                 tmpflnm=None,
                 do_read=0,
                 ftracerinfo=gcdf.def_tracerinfo_file,
                 fdiaginfo=gcdf.def_diaginfo_file,
                 sel_categorys=None, sel_tracers=None, sel_taus=None):
        """	initialize
        augments:
        flnm, string, the bpch file name
        mode, string,  w/r
        maxtracer, integer,  the max number of tracers
        tmpflnm,   string,   file name for  temp information
        """

        self.flnm = flnm
        self.cur_tracer = -1
        self.maxtracer = maxtracer
        self.title = None
        self.modelname = None

        self.data = list()
        self.stat = -1
        self.ntracers = 0
        self.read_mode = -1

        # grid

        self.tracerinfo = tracer_info(ftracerinfo)
        self.diaginfo = diag_info(fdiaginfo)

        if (tmpflnm is None):
            flnm_tmp = flnm.strip() + '_info_py'
        else:
            flnm_tmp = tmpflnm.strip()

        if (mode == 'r'):
            if (do_read == 0):  # only the the heads will be first
                ios, title, tracer_id, lonres, latres,\
                    ix, jx, lx,\
                    ifirst, jfirst, lfirst,\
                    halfpolar, centre180, \
                    tau0, tau1, ntracers = \
                    bpch2_rw_py.read_bpch2_head(flnm, flnm_tmp, maxtracer)

                # print bpch2_rw_py.__doc_

                if (ios > 0):
                    print('error in read data', ios)
                    self.read_mode = -1
                    self.stat = ios
                else:
                    self.title = title.strip()

                    tf = open(flnm_tmp, 'r')
                    lines = tf.readlines()
                    tf.close()

                    ic = 0

                    category = list()
                    modelname = list()
                    unit = list()
                    reserved = list()
                    if (ntracers == 1):

                        line = lines[0]
                        print(type(line))
                        print('line', line)
                        terms = line.split(',')
                        category.append(terms[0])
                        modelname.append(terms[1])
                        unit.append(terms[2])
                        reserved.append(terms[3].strip())
                    else:
                        for line in lines:
                            line = line.replace('\n', '')
                            terms = line.split(',')
                            category.append(terms[0])
                            modelname.append(terms[1])
                            unit.append(terms[2])
                            reserved.append(terms[3])

                    self.title = title
                    self.stat = 0
                    self.data = list()
                    self.ntracers = ntracers

                    """ store the data into the bpch2_data class """
                    for isp in range(ntracers):
                        data_grid = bpgrid(
                            ix[isp],
                            jx[isp],
                            lx[isp],
                            ifirst[isp],
                            jfirst[isp],
                            lfirst[isp],
                            halfpolar[isp],
                            centre180[isp],
                            lonres[isp],
                            latres[isp])

                        offset = self.diaginfo.get_offset(category[isp])
                        if (offset < 0):
                            offset = 0

                        real_id = tracer_id[isp] + offset

                        traname, trafullname, tramolew, trascale, traunit = self.tracerinfo.get_tracer_info(
                            real_id)
                        pbdata = bpch2_data(
                            tracer_id[isp],
                            category[isp],
                            unit[isp],
                            data_grid,
                            None,
                            tau0=tau0[isp],
                            tau1=tau1[isp],
                            modelname=modelname[isp],
                            reserved=reserved[isp],
                            offset=offset,
                            name=traname,
                            fullname=trafullname,
                            molew=tramolew,
                            scale=trascale,
                            traunit=traunit,
                            id=real_id)

                        self.data.append(pbdata)
                    self.read_mode = do_read
            if (do_read == 1):
                print(self.flnm)

                funit = 199
                fti, title, stat = bpch2_rw_py.open_bpch2_for_read(
                    funit, self.flnm)
                self.title = title.strip()

                if (stat != 0):
                    print('error in read in bpch2_rw_py:', stat)
                    self.stat = self.stat
                else:
                    vtracer_id, vhalfpolar, vcentre180, vni, vnj, vnl, vifirst, vjfirst, vlfirst, vlonres, vlatres,\
                        vtau0, vtau1, vmodelname, vcategory, vunit, vreserved, vdata_array, stat = bpch2_rw_py.read_bpch2_record(funit)
                    while (stat == 0):

                        data_grid = bpgrid(vni, vnj, vnl,
                                           vifirst, vjfirst, vlfirst,
                                           vhalfpolar, vcentre180,
                                           vlonres, vlatres)
                        offset = self.diaginfo.get_offset(vcategory)
                        if (offset < 0):
                            offset = 0

                        real_id = vtracer_id + offset
                        traname, trafullname, tramolew, trascale, traunit = self.tracerinfo.get_tracer_info(
                            real_id)
                        vdata = zeros([vni, vnj, vnl], float)
                        vdata[0:vni,
                              0:vnj,
                              0:vnl] = vdata_array[(vifirst - 1):(vifirst + vni - 1),
                                                   (vjfirst - 1):(vjfirst + vnj - 1),
                                                   (vlfirst - 1):(vlfirst + vnl - 1)]
                        pbdata = bpch2_data(
                            vtracer_id,
                            vcategory.strip(),
                            vunit.strip(),
                            data_grid,
                            vdata,
                            tau0=vtau0,
                            tau1=vtau1,
                            modelname=vmodelname.strip(),
                            reserved=vreserved,
                            offset=offset,
                            name=traname.strip(),
                            fullname=trafullname.strip(),
                            molew=tramolew,
                            scale=trascale,
                            traunit=traunit.strip(),
                            id=real_id)

                        self.data.append(pbdata)
                        self.ntracers = self.ntracers + 1
                        # print traname, trafullname, self.ntracers

                        vtracer_id, vhalfpolar, vcentre180, vni, vnj, vnl, vifirst, vjfirst, vlfirst, vlonres, vlatres,\
                            vtau0, vtau1, vmodelname, vcategory, vunit,\
                            vreserved, vdata_array, stat = bpch2_rw_py.read_bpch2_record(funit)
                        # print 'vtracer_id', vtracer_id, traname, trafullname,
                        # tramolew, trascale, traunit

                    if (stat == -1 or stat == 29):
                        self.stat = 0
                        self.read_mode = do_read
                    else:
                        self.stat = stat
                        self.read_mode = -1

                    stat = bpch2_rw_py.close_bpch2_file(funit)
            if (do_read == 2):
                funit = 199
                print('mode read 2', self.flnm)

                fti, title, stat = bpch2_rw_py.open_bpch2_for_read(
                    funit, self.flnm)
                self.title = title.strip()

                if (stat != 0):
                    print('error in read bpch2:', stat, self.flnm)
                    self.stat = self.stat
                else:
                    stat = bpch2_rw_py.close_bpch2_file(funit)

                    for iname in range(len(sel_categorys)):
                        cname = sel_categorys[iname]
                        tid = sel_tracers[iname]
                        tau0 = sel_taus[iname]
                        fti, title, stat = bpch2_rw_py.open_bpch2_for_read(
                            funit, self.flnm)
                        vtracer_id, vhalfpolar, vcentre180, vni, vnj, vnl, vifirst, vjfirst, vlfirst, vlonres, vlatres,\
                            vtau0, vtau1, vmodelname, vcategory, vunit, vreserved, vdata_array, stat = \
                            bpch2_rw_py.sel_bpch2_record(funit, cname, tid, tau0)

                        # print 'stat in read mode 2', stat, self.flnm, cname
                        stat = bpch2_rw_py.close_bpch2_file(funit)

                        if (stat == 0):
                            data_grid = bpgrid(vni, vnj, vnl,
                                               vifirst, vjfirst, vlfirst,
                                               vhalfpolar, vcentre180,
                                               vlonres, vlatres)
                            offset = self.diaginfo.get_offset(vcategory)
                            if (offset < 0):
                                offset = 0

                            real_id = vtracer_id + offset
                            traname, trafullname, tramolew, trascale, traunit = self.tracerinfo.get_tracer_info(
                                real_id)
                            vdata = zeros([vni, vnj, vnl], float)
                            vdata[0:vni,
                                  0:vnj,
                                  0:vnl] = vdata_array[(vifirst - 1):(vifirst + vni - 1),
                                                       (vjfirst - 1):(vjfirst + vnj - 1),
                                                       (vlfirst - 1):(vlfirst + vnl - 1)]
                            pbdata = bpch2_data(
                                vtracer_id,
                                vcategory.strip(),
                                vunit.strip(),
                                data_grid,
                                vdata,
                                tau0=vtau0,
                                tau1=vtau1,
                                modelname=vmodelname.strip(),
                                reserved=vreserved,
                                offset=offset,
                                name=traname.strip(),
                                fullname=trafullname.strip(),
                                molew=tramolew,
                                scale=trascale,
                                traunit=traunit.strip(),
                                id=real_id)

                            self.data.append(pbdata)
                            self.ntracers = self.ntracers + 1

                        else:

                            break

                        # print 'tracers ', self.ntracers, trafullname,
                        # vcategory.strip(), vtracer_id

                        stat = bpch2_rw_py.close_bpch2_file(funit)

    def open_bpch2_w(self, funit, title=""):
        self.funit = funit
        print('funit', funit)
        stat = bpch2_rw_py.open_bpch2_for_write(funit, self.flnm, title)
        self.stat = stat
        return stat

    def close_bpch2_file(self):
        if (self.funit is not None):
            stat = bpch2_rw_py.close_bpch2_file(self.funit)
            self.stat = stat
            self.funit = None

    def write_bpch2_data(self, data_id):
        stat = -1
        if (data_id < self.ntracers and self.funit is not None):
            bpdata = self.data[data_id]
            stat = bpdata.write_to_bpch2_file(self.funit)
        return stat

    def add_bpch2_data(self, bpdata, data_id=None):
        if (data_id is None):
            self.data.append(bpdata)
            self.ntracers = self.ntracers + 1
        else:
            self.data[data_id] = bpdata

    def get_data(self, categorys=None, tracers=None, taus=None, tranames=None):
        """ search the records for all the data

        """
        if (self.ntracers <= 0):
            print('tracer not found')
            return None
        sel_data = list()
        nval = 0
        array_cat = None
        array_tracer = None
        array_tau = None
        array_traname = None

        if (categorys is not None):
            array_cat = array(categorys)
            nval = size(array_cat)

        if (tracers is not None):
            array_tracer = array(tracers)
            nval = size(array_tracer)

        if (taus is not None):
            array_tau = array(taus)
            nval = size(array_tau)

        if (tranames is not None):
            array_traname = array(tranames)
            nval = size(array_traname)

        if (nval == 0):
            nval = len(self.data)
            founded = ones(nval)
            return self.data, founded

        found_data = list()
        founded = zeros(nval)

        for bpdata in self.data:
            score_cat = ones(nval)
            score_tracer = ones(nval)
            score_tau = ones(nval)
            score_traname = ones(nval)

            if (array_cat is not None):
                score_cat = where(array_cat == bpdata.category.strip(), 1, 0)
                #	print '1', score_cat
            if (array_tracer is not None):
                score_tracer = where(array_tracer == int(bpdata.ntracer), 1, 0)
                # print '2', score_tracer
            if (array_tau is not None):
                tau0, tau1 = bpdata.get_attr(['tau0', 'tau1'])
                score_tau = where(
                    logical_and(
                        array_tau >= tau0,
                        array_tau < tau1),
                    1,
                    0)

                # print '3', score_tau
                # print array_tau, tau0

            if (array_traname is not None):
                traname, tau1 = bpdata.get_attr(['name', 'tau1'])
                # print 'traname', traname
                score_traname = where(array_traname == traname.strip(), 1, 0)
                # print '4', score_traname

            final_score = score_cat * score_tracer * score_tau * score_traname
            founded = founded + final_score

            if (sum(final_score) > 0):
                if (bpdata.data is None):
                    # read in the required data
                    ix = bpdata.grid.ix
                    jy = bpdata.grid.jx
                    lz = bpdata.grid.lx

                    tracer = bpdata.ntracer
                    xtau0, xtau1 = bpdata.get_attr(['tau0', 'tau1'])
                    category = bpdata.category
                    data, dunit, ios = bpch2_rw_py.read_bpch2(
                        self.flnm, category, tracer, xtau0, ix, jy, lz)
                    print('IOS', ios)

                    if (ios > 0):
                        print('error in read data', ios)
                        return None
                    elif (ios == -1):
                        print('no data found')
                        return None
                    print('after read', size(data))
                    bpdata.update_data(data)
                found_data.append(bpdata)
        return found_data, founded

    def read_bpch2_data(self, data_id, save_data=False):
        """ read in bpch2 data record for a given tracer nu print 'ix', bpch2.ix[:]mber

        """
        if (data_id > self.ntracers):
            print('tracer not found')
            return None
        if (self.read_mode == 1):
            bpdata = self.data[data_id]
            return bpdata

        bpdata = self.data[data_id]

        ix = bpdata.grid.ix
        jy = bpdata.grid.jx
        lz = bpdata.grid.lx

        tracer = bpdata.ntracer
        tau0, tau1 = bpdata.get_attr(['tau0', 'tau1'])
        category = bpdata.category
        data, dunit, ios = bpch2_rw_py.read_bpch2(self.flnm, category, tracer,
                                                  tau0, ix, jy, lz)

        if (ios > 0):
            print('error in read data', ios)
            return None
        elif (ios == -1):
            print('no data found')
            return None

        if (save_data):
            bpdata.update_data(data)
            return bpdata
        else:
            return data

    def print_datainfo(self, data_id=None):
        if (data_id is not None and data_id < self.ntracers):
            bpdata = self.data[data_id]
            ix = bpdata.grid.ix
            jy = bpdata.grid.jx
            lz = bpdata.grid.lx
            tracer = bpdata.ntracer
            unit = bpdata.unit
            category = bpdata.category
            tau0, tau1, modelname, traname, unit, traid = bpdata.get_attr(
                ['tau0', 'tau1', 'modelname', 'name', 'traunit', 'id'])
            unit = bpdata.unit
            print(data_id, tracer, traid, traname.strip(), category.strip(), ix, jy, lz, unit.strip())

        else:
            for data_id in range(self.ntracers):

                bpdata = self.data[data_id]
                ix = bpdata.grid.ix
                jy = bpdata.grid.jx
                lz = bpdata.grid.lx
                tracer = bpdata.ntracer
                category = bpdata.category

                tau0, tau1, modelname, traname, unit, traid = bpdata.get_attr(
                    ['tau0', 'tau1', 'modelname', 'name', 'traunit', 'id'])
                unit = bpdata.unit
                print(data_id + 1, tracer, traid, traname.strip(), category.strip(), ix, jy, lz, unit.strip())

    def open_bpch2_for_write(self, title, iunit):
        state = bpch2_rw_py.open_bpch2_for_write(iunit, self.flnm, title)
        if (state == 0):
            self.funit = iunit

    def write_bpch2_head(self, funit, title):
        stat = bpch2_wr_py.write_bpch2_hdr(iunit, title)
        return stat


if (__name__ == '__main__'):

    bpch2 = bpch2_file_rw("x2ctm.bpch", 'r', do_read=1)
    print('ntrac', bpch2.ntracers)

    bpch2.print_datainfo()
    data = bpch2.read_bpch2_data(30, True)
    print(shape(data.data))
    print(data.grid.centre180)

    # data.display(0,show_map=1, map_proj='spstere', show_lat=0, boundinglat=45)
    # data.display(0,show_map=1, map_proj='ortho', show_lat=1, show_lon=1, lat_0=-90.0,  maxlat=-30)

    data.display(0, show_map=1)
    data_list, founded = bpch2.get_data(None, None, 160056.0, 'PSURF')
    print('founded', founded, len(data_list))
    for bpdata in data_list:
        traname, tau0, tau1 = bpdata.get_attr(['name', 'tau0', 'tau1'])
        print(bpdata.category, tau0, tau1, traname.strip())
        print(bpdata.data[3, 3, 0])
    test_write = True
    if (test_write):
        bpch2_w = bpch2_file_rw("restart_new", 'w')
        ifunit = 95
        bpch2_w.open_bpch2_w(ifunit, "test")
        bpch2_w.add_bpch2_data(data)
        print('I am here')
        bpch2_w.add_bpch2_data(data)
        data = bpch2.read_bpch2_data(20, True)
        print('the second read')
        bpch2_w.add_bpch2_data(data)
        print('try to write')
        bpch2_w.write_bpch2_data(0)
        bpch2_w.write_bpch2_data(1)
        bpch2_w.close_bpch2_file()
