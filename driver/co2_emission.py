from pylab import *
from numpy import *
import  geo_constant as gc
import geos_chem_def as gcfs
import bpch2_rw_v2 as brw
import bpch2_rw_py as bp
import grid_py as gm
import geos_read_input as gri
import time_module as tm
import gen_plots as gpl
import gc_dist as gcd
import numpy.linalg as nlg



LFOSSCO2     = True
LBIOCO2      = True
LOCCO2       = True
LBIOBRNCO2   = True
LBIOFUELCO2  = True
LBIONETCO2   = False 
LUSECASANEP  = False

def GET_GRID(return_edge=False):
    #  print 'get_grid'
    gm.grid_mod.compute_grid()
    # print 'try setup grid'
    IIPAR=gm.grid_mod.get_x_size()
    JJPAR=gm.grid_mod.get_y_size()

    JY=arange(1, JJPAR+1)
    YY = gm.grid_mod.get_ymids( JY )
    IX=arange(1, IIPAR+1)
    XX = gm.grid_mod.get_xmids( IX )
    sf_area=1.0E4*gm.grid_mod.get_areas(IX, JY)
    
    IX=arange(1, IIPAR+2)
    Iy=arange(1, JJPAR+2)
    xedge=gm.grid_mod.get_xedges(IX)
    yedge=gm.grid_mod.get_yedges(JY)
    if (return_edge):
        return XX, YY, xedge, yedge, sf_area
    else:
        return XX, YY



def DEFINE_FF_CO2_REGIONS(X,Y):
    REGION=zeros(shape(X))
    REGION=REGION+18

    #  Region 2 -- Western Canada      
    vx=logical_and( X >= -141.0, X < -110.0)
    vy=logical_and( Y >=   49.0,  Y <  70.0)
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=8
        
    # Region #3 -- Central Canada                 
    
    vx=logical_and(  X >= -110.0, X <  -88.0 )
    vy=logical_and( Y >=   49.0, Y <   70.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=9
    
    
    
    #4 -- Eastern Canada (1st sub-box) 
    vx=logical_and( X >= -88.0, X < -60.0 )
    vy=logical_and( Y >=  47.0, Y <  60.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=10

    #4 -- Eastern Canada (1st sub-box) 
    vx=logical_and(X >=  -60.0, X < -50.0)
    vy=logical_and(Y >=   43.0, Y <  56.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=10
        
    
    #5 -- USA
    vx=logical_and(X >= -165.0, X < -60.0)
    vy=logical_and(Y >=  26.0, Y <  49.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=11
    
    #5 -- Asia
    vx=logical_and(X >=   95.0, X < 150.0)
    vy=logical_and(Y >=   10.0, Y <  50.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=12
   

    return  REGION

def DEFINE_BB_CO2_REGIONS(X,Y):
    REGION=zeros(shape(X))
    REGION=REGION+18

    #  Region 2 -- Western Canada      
    vx=logical_and( X >= -141.0, X < -110.0)
    vy=logical_and( Y >=   49.0,  Y <  70.0)
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=13

    # Region #3 -- Central Canada                 
    
    vx=logical_and(  X >= -110.0, X <  -88.0 )
    vy=logical_and( Y >=   49.0, Y <   70.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=14
    

    
    #4 -- Eastern Canada (1st sub-box) 
    vx=logical_and( X >= -88.0, X < -60.0 )
    vy=logical_and( Y >=  47.0, Y <  60.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=15

    #4 -- Eastern Canada (2nd sub-box) 
    vx=logical_and(X >=  -60.0, X < -50.0)
    vy=logical_and(Y >=   43.0, Y <  56.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=15
        
    
    #5 -- USA
    vx=logical_and(X >= -165.0, X < -60.0)
    vy=logical_and(Y >=  26.0, Y <  49.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=11
    
    #5 -- Asia
    vx=logical_and(X >=   95.0, X < 150.0)
    vy=logical_and(Y >=   10.0, Y <  50.0 )
    idx=where(logical_and(vx, vy))
    if (size(idx)>0):
        REGION[idx]=17
    
        
    return  REGION

def define_new_regions(lon, lat, old_map, lon_idx, lat_idx, reg_vals):
    """ defined new regions following the given values"""
    
    nz=max(old_map.flat)
    lon_pos=searchsorted(lon_idx, lon)
    lat_pos=searchsorted(lat_idx, lat)
    lon_pos=where(lon_pos==size(lon_idx), 0, lon_pos)
    lat_pos=where(lat_pos==size(lat_idx), 0, lat_pos)
    # all_var=dot(lon_pos[:, newaxis], lat_pos)
    # all_var=where(all_var>0, all_var+nz, all_var)
    new_map=new_reg_val[lon_pos,:]
    new_map=new_map[:, lat_pos]
    new_map=where(new_map>0, new_map+nz, new_map)
    
    return new_map


def define_surface_co2_region_hr(x, y, xedge, yedge, \
                                 sel_id, nlon_reg, nlat_reg, \
                                 do_debug=False):

    """ the state variables define following transcom definiation but increase some new regions
        see http://www.purdue.edu/transcom/images/smoothmap2_final_2.jpg
    
    """
    
    import read_map
    
    
    lon, lat, type_map=read_map.read_map(360, 180)
    land_reg=['ice',\
              'North American Boreal',\
              'North American Temperate',\
              'South American Tropical',\
              'South American Temperate',\
              'North Africa',\
              'South Africa',\
              'Eurasia Boreal', \
              'Eurasia Temperate',\
              'Tropical Asia', \
              'Australia',  \
              'Europe']

    ocean_reg=['North Pacific Temperate', \
               'West Pacific Tropics', \
               'East Pacific Tropics', \
               'South Pacific Temperate',\
               'Northern Ocean',
               'Northen Atlantic Temperate',\
               'Atlantic Tropics', \
               'South Atlantic Temperate',\
               'South Ocean', \
               'Indian Tropical', \
               'South Indian Temperate']
    
    type_map=type_map.astype(int)
    nz=max(type_map.flat)
    type_map_sav=array(type_map)
    JJPAR=gm.grid_mod.get_y_size()
    
    JY=arange(1, JJPAR+1)
    YY = gm.grid_mod.get_ymids( JY )

    idx=where(type_map==sel_id)
    lon_reg_id=idx[0]
    lat_reg_id=idx[1]

    rlon_reg=lon[lon_reg_id]
    rlat_reg=lat[lat_reg_id]
    min_lon_reg, max_lon_reg=min(rlon_reg), max(rlon_reg)
    min_lat_reg, max_lat_reg=min(rlat_reg), max(rlat_reg)
    
    id_lon=where(logical_and(lon>=min_lon_reg, lon<=max_lon_reg))
    id_lon=squeeze(id_lon)
    id_lat=where(logical_and(lat>=min_lat_reg, lat<=max_lat_reg))
    id_lat=squeeze(id_lat)
    
    sel_id_list=list()
    areas=list()
    new_map=zeros(shape(type_map))
    for ilat in id_lat:
        sel_idx=where(lat_reg_id==ilat)
        sel_idx=squeeze(sel_idx)
        
        if (size(sel_idx)>0):
            ix, iy=lon_reg_id[sel_idx], lat_reg_id[sel_idx]
            area_y=abs(cos(pi*lat[ilat]/180.0))*size(sel_idx)
            areas.append(area_y)
        else:
            areas.append(0)
        sel_id_list.append(array(sel_idx))

    areas=array(areas)
    area_div=sum(areas)/nlat_reg
    area_sum=0.0
    nlat_div=size(id_lat)
    area_lon=zeros(size(id_lon), float)
    # print shape(area_lon)
    # print id_lon
    # print lon[id_lon]
    # print id_lat
    # print lat[id_lat]
    
    sub_sel_idx=list()
    region_count=0
    print 'nlat_div', nlat_div, area_div
    new_reg=list()
    for ilat in range(nlat_div):
        area_sum=area_sum+areas[ilat]
        # print ilat, area_sum
        
        # spread area into lon
        sel_idx=sel_id_list[ilat]
        if (size(sel_idx)>0):
            tmp_lon_pos=lon_reg_id[sel_idx]
            # print id_lon[40:50]
            # print tmp_lon_pos
            lon_pos=searchsorted(id_lon, tmp_lon_pos)
            lon_pos=squeeze(lon_pos)
            # print lon_pos
            # print id_lon[lon_pos]
            # print shape(area_lon)
            # print max(lon_pos)
            
            
            cell_area=abs(cos(pi*lat[lat_reg_id[sel_idx[0]]]/180.0))
            for ii in range(size(sel_idx)):
                area_lon[lon_pos[ii]]=area_lon[lon_pos[ii]]+cell_area
                sub_sel_idx.append(sel_idx[ii])
            
        if (area_sum>=area_div or ilat==nlat_div-1):
            sub_sel_idx=array(sub_sel_idx)
            
            # divide the along the longitude
            area_lon_div=sum(area_lon)/nlon_reg
            area_lon=add.accumulate(area_lon)
            lon_div_ar=arange(0, nlon_reg,1.0)
            lon_div_ar=area_lon_div*lon_div_ar
            lon_pos_2=searchsorted(lon_div_ar, area_lon)
            lon_pos_2=lon_pos_2-1
            
            
            ix, iy=lon_reg_id[sub_sel_idx], lat_reg_id[sub_sel_idx]
            print shape(id_lon)
            # print shape(ix), shape(iy)
            
            tmp_ix_pos=searchsorted(id_lon, ix)
            # print 'max lon_pos_2', max(lon_pos_2)
            # print lon_div_ar
            # print area_lon[tmp_ix_pos]
            # print lon_pos_2[tmp_ix_pos]
            # print lon_pos_2
            
            print 'nz', nz
            new_reg_id=lon_pos_2[tmp_ix_pos]+(region_count)*nlon_reg
            print min(new_reg_id), max(new_reg_id), region_count
            
            # new_reg_id=where(new_reg_id==0, 0, new_reg_id+nz)
            new_map[ix, iy]=new_reg_id  # -type_map[ix,iy]
            region_count=region_count+1
            area_sum=0.0
            
            area_lon=zeros(size(id_lon), float)
            sub_sel_idx=list()
    
    type_map=where(type_map>sel_id, type_map+nlon_reg*nlat_reg-1, type_map)
    type_map=type_map+new_map
                       
                              
                 
                 
                            
       
    # now we add some new regions into it
    
    new_map=read_map.regrid_map(lon, lat, xedge, yedge,type_map)
    nx=size(xedge)-1
    ny=size(yedge)-1
    new_map=new_map[0:nx, 0:ny]
    print max(new_map.flat)

    #  figure(2)
    # gpl.plot_map(type_map,lon, lat,  use_pcolor=1, cmap=cm.Paired) #Set1)
    # savefig('test_div.png')
    # savefig('test_div.ps')
    # show()
    
    if (do_debug):
        subplot(2,1,1)
        nz=max(type_map.flat)
        type_map=1.0*type_map
        
        gpl.plot_map(type_map,lon, lat,  use_pcolor=1)
        subplot(2,1,2)
        gpl.plot_map(new_map, x, y,  use_pcolor=1)
        show()

        
    # the emission error are taken transcom prior flux file 3 for basis function error while the ice one is assigned to a small value (percentage) #  
    def_err_land=[0.15, \
             0.726200, 1.500000, 1.412200,  \
             1.227300, 1.331700, 1.411400, \
             1.514300, 1.728300, 0.865600, 0.593600, 1.420900]
    #   def_err_ocean=[0.820000, \
    #            0.500000, 0.560000, 1.220000, \
    #            0.260000, 0.400000, 0.400000, \
    #            0.480000, 1.500000, 0.740000, \
    #            0.540000] # in Gt C /year
    def_err_ocean=[0.270000, \
                   0.390000, 0.370000, 0.630000, \
                   0.350000, 0.270000, 0.410000, \
                   0.550000, 0.720000, 0.480000, \
                   0.410000] # in Gt C /year
    def_err_land[sel_id]=def_err_land[sel_id]/(sqrt(nlon_reg)*sqrt(nlat_reg))
    all_def_err_land=def_err_land[0:sel_id+1]+[def_err_land[sel_id]]*(nlon_reg*nlat_reg-1)+def_err_land[sel_id+1:]
    for ireg in range(1, nlon_reg*nlat_reg):
        regname='SUB_R%2.2dS%2.2d' % (sel_id, ireg)
        new_reg.append(regname)
    print new_reg
    all_land_reg=land_reg[0:sel_id+1]+new_reg+land_reg[sel_id+1:]
    def_err=all_def_err_land+def_err_ocean
    return new_map, all_land_reg, ocean_reg, def_err


def define_surface_co2_region_cell(x, y, xedge, yedge, \
                                   sel_id, nlon_reg, nlat_reg, \
                                 do_debug=False):

    """ the state variables define following transcom definiation but increase some new regions
        see http://www.purdue.edu/transcom/images/smoothmap2_final_2.jpg
    
    """
    
    import read_map
    
    
    lon, lat, type_map=read_map.read_map(360, 180)
    land_reg=['ice',\
              'North American Boreal',\
              'North American Temperate',\
              'South American Tropical',\
              'South American Temperate',\
              'North Africa',\
              'South Africa',\
              'Eurasia Boreal', \
              'Eurasia Temperate',\
              'Tropical Asia', \
              'Australia',  \
              'Europe']

    ocean_reg=['North Pacific Temperate', \
               'West Pacific Tropics', \
               'East Pacific Tropics', \
               'South Pacific Temperate',\
               'Northern Ocean',
               'Northen Atlantic Temperate',\
               'Atlantic Tropics', \
               'South Atlantic Temperate',\
               'South Ocean', \
               'Indian Tropical', \
               'South Indian Temperate']
    
    type_map=type_map.astype(int)
    nz=max(type_map.flat)
    type_map_sav=array(type_map)
    JJPAR=gm.grid_mod.get_y_size()

    
    JY=arange(1, JJPAR+1)
    YY = gm.grid_mod.get_ymids( JY )

    new_map=read_map.regrid_map(lon, lat, xedge, yedge,type_map)
    nx=size(xedge)-1
    ny=size(yedge)-1
    new_map=new_map[0:nx, 0:ny]
    print max(new_map.flat)
    type_map=array(new_map)
    
    lon=array(xedge[0:nx])
    lat=array(xedge[0:ny])
    new_map=zeros(shape(type_map))
    idx=where(type_map==sel_id)
    idx=squeeze(idx)
    cnt_cell=0
    idx_lon=idx[0]
    idx_lat=idx[1]
    npt=size(idx_lon)
    for pt in range(npt):
        new_map[idx_lon[pt], idx_lat[pt]]=cnt_cell
        cnt_cell=cnt_cell+1
    
    print 'cnt_cell', cnt_cell
    
    type_map=where(type_map>sel_id, type_map+cnt_cell-1, type_map)
    type_map=type_map+new_map
    
                       
    if (do_debug):
        subplot(2,1,1)
        nz=max(type_map.flat)
        type_map=1.0*type_map
        
        gpl.plot_map(type_map,lon, lat,  use_pcolor=1)
        subplot(2,1,2)
        gpl.plot_map(new_map, x, y,  use_pcolor=1)
        show()

        
    # the emission error are taken transcom prior flux file 3 for basis function error while the ice one is assigned to a small value (percentage) #  
    def_err_land=[0.15, \
             0.726200, 1.500000, 1.412200,  \
             1.227300, 1.331700, 1.411400, \
             1.514300, 1.728300, 0.865600, 0.593600, 1.420900]
    #   def_err_ocean=[0.820000, \
    #            0.500000, 0.560000, 1.220000, \
    #            0.260000, 0.400000, 0.400000, \
    #            0.480000, 1.500000, 0.740000, \
    #            0.540000] # in Gt C /year
    def_err_ocean=[0.270000, \
                   0.390000, 0.370000, 0.630000, \
                   0.350000, 0.270000, 0.410000, \
                   0.550000, 0.720000, 0.480000, \
                   0.410000] # in Gt C /year
    
    def_err_land[sel_id]=def_err_land[sel_id]/(sqrt(cnt_cell))
    all_def_err_land=def_err_land[0:sel_id+1]+[def_err_land[sel_id]]*(cnt_cell-1)+def_err_land[sel_id+1:]

    new_reg=list()
    for ireg in range(1, cnt_cell):
        regname='SUB_R%2.2dS%2.2d' % (sel_id, ireg)
        new_reg.append(regname)
    print new_reg
    
    all_land_reg=land_reg[0:sel_id+1]+new_reg+land_reg[sel_id+1:]
    print len(all_land_reg)
    def_err=all_def_err_land+def_err_ocean
    # subplot(2,1,1)
    # gpl.plot_map(type_map,x, y,  use_pcolor=1,cmap=cm.jet)
    # subplot(2,1,2)
    # gpl.plot_map(new_map, x, y,  use_pcolor=1,cmap=cm.jet)
    # show()
    new_map=array(type_map)
    
    return new_map, all_land_reg, ocean_reg, def_err
    






def define_surface_co2_region(x, y, xedge, yedge, do_debug=False):

    """ the state variables define following transcom definiation
    see http://www.purdue.edu/transcom/images/smoothmap2_final_2.jpg
    
    """
    
    import read_map
    
    
    lon, lat, type_map=read_map.read_map(360, 180)
    land_reg=['ice',\
              'North American Boreal',\
              'North American Temperate',\
              'South American Tropical',\
              'South American Temperate',\
              'North Africa',\
              'South Africa',\
              'Eurasia Boreal', \
              'Eurasia Temperate',\
              'Tropical Asia', \
              'Australia',  \
              'Europe']

    ocean_reg=['North Pacific Temperate', \
               'West Pacific Tropics', \
               'East Pacific Tropics', \
               'South Pacific Temperate',\
               'Northern Ocean',
               'Northen Atlantic Temperate',\
               'Atlantic Tropics', \
               'South Atlantic Temperate',\
               'South Ocean', \
               'Indian Tropical', \
               'South Indian Temperate']
    
    type_map=type_map.astype(int)
    nz=max(type_map.flat)
    type_map_sav=array(type_map)
    JJPAR=gm.grid_mod.get_y_size()

    JY=arange(1, JJPAR+1)
    YY = gm.grid_mod.get_ymids( JY )
        
    new_map=read_map.regrid_map(lon, lat, xedge, yedge,type_map_sav)
    nx=size(xedge)-1
    ny=size(yedge)-1
    new_map=new_map[0:nx, 0:ny]

    if (do_debug):
        subplot(2,1,1)
        nz=max(type_map.flat)
        type_map=1.0*type_map
        
        gpl.plot_map(type_map,lon, lat,  use_pcolor=1)
        subplot(2,1,2)
        gpl.plot_map(new_map, x, y,  use_pcolor=1)
        show()

        
    # the emission error are taken transcom prior flux file 3 for basis function error while the ice one is assigned to a small value (percentage) #  
    def_err=[0.15, \
             0.726200, 1.500000, 1.412200,  \
             1.227300, 1.331700, 1.411400, \
             1.514300, 1.728300, 0.865600, \
             0.593600, 1.420900, 0.820000, \
             0.500000, 0.560000, 1.220000, \
             0.260000, 0.400000, 0.400000, \
             0.480000, 1.500000, 0.740000, \
             0.540000] # in Gt C /year
        
    return new_map, land_reg, ocean_reg, def_err




            
def  READ_ANNUAL_FOSSILCO2(DATA_DIR, nx, ny, CO2_DIR='CO2_200508', flnmpref='fossil03_CO2', tracer=1,category='CO2-SRCE'):
    #  3/4 only 2x25 for now...
    ext2d=bp.get_name_ext_2d()
    extres=bp.get_res_ext()
    
    FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+flnmpref.strip()+'.'+ext2d.strip()+'.'+extres.strip()
    # print 'read data from ', FILENAME
    
    
    tau0=tm.get_tau(1985,1,1)
    data, dunit, ios=bp.read_bpch2(FILENAME, category, tracer,tau0, \
                           nx,  ny, 1)
    
    if (ios==0):
        data=squeeze(data)
        return data,dunit
    else:
        print 'error in read data ', ios
        print FILENAME
        return None, ""
    
def  READ_ANNUAL_OCEANCO2(DATA_DIR, nx, ny, CO2_DIR='CO2_200508', flnmpref='ocean_CO2',  tracer=2,category='CO2-SRCE'):
    #  3/4 only 2x25 for now...
    ext2d=bp.get_name_ext_2d()
    extres=bp.get_res_ext()
    
    FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+flnmpref.strip()+'.'+ext2d.strip()+'.'+extres.strip()
    #  print 'read data from ', FILENAME
    
    tau0=tm.get_tau(1985,1,1)
    data, dunit, ios=bp.read_bpch2(FILENAME, category, tracer,tau0, \
                           nx,  ny, 1)
    
    if (ios==0):
        data=squeeze(data)
        return data,dunit
    else:
        print 'error in read data ', ios
        print FILENAME
        return None, ""

def  READ_ANNUAL_BIOFUELCO2(DATA_DIR, nx, ny, CO2_DIR='CO2_200508', flnmpref='/biofuel_CO2',  tracer=5,category='CO2-SRCE'):
    #  3/4 only 2x25 for now...
    ext2d=bp.get_name_ext_2d()
    extres=bp.get_res_ext()
    
    FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+flnmpref.strip()+'.'+ext2d.strip()+'.'+extres.strip()
    # print 'read data from ', FILENAME
    
    tau0=tm.get_tau(1985,1,1)
    data, dunit, ios=bp.read_bpch2(FILENAME, category, tracer,tau0, \
                           nx,  ny, 1)
    # print shape(data)
    # print ios
    
    if (ios==0):
        data=squeeze(data)
        return data,dunit
    else:
        print 'error in read data ', ios
        print FILENAME
        return None, ""

def  READ_ANNUAL_BIONET_CO2(DATA_DIR, nx, ny, CO2_DIR='CO2_200508', flnmpref='//net_terr_exch_CO2',  tracer=5,category='CO2-SRCE'):
    #  3/4 only 2x25 for now...
    ext2d=bp.get_name_ext_2d()
    extres=bp.get_res_ext()
    FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+flnmpref.strip()+'.'+ext2d.strip()+'.'+extres.strip()+".txt"
    # print 'read data from ', FILENAME
    
    tau0=tm.get_tau(1985,1,1)
    data, dunit, ios=bp.read_bpch2(FILENAME, category, tracer,tau0, \
                           nx,  ny, 1)
    
    if (ios==0):
        data=squeeze(data)
        return data,dunit
    
    else:
        print 'error in read data ', ios
        print FILENAME
        return None, ""


def  READ_BBIO_DAILYAVERAGE(DATA_DIR, nx, ny, \
                            MONTH, DAY, YYYY, \
                            CO2_DIR='CO2_200508/BBIO_DAILYAVG', \
                            flnmpref='/CO2.daily',  tracer=3,category='CO2-SRCE'):
    #  3/4 only 2x25 for now...
    ext2d=bp.get_name_ext_2d()
    extres=bp.get_res_ext()
    DOY=tm.day_of_year(YYYY, MONTH,DAY)
    if (DOY>59 and mod(YYYY,4)>0): 
            DOY=DOY+1

    SDOY=r'%3.3d' % DOY
    
    
    FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+flnmpref.strip()+'.'+ext2d.strip()+'.'+extres.strip()+"."+SDOY
    # print 'read data from ', FILENAME
    
    tau0=tm.get_tau(2000,MONTH,DAY)
    tau0=tau0/(3600.0)
    data, dunit, ios=bp.read_bpch2(FILENAME, category, tracer,tau0, \
                           nx,  ny, 1)
    
    if (ios==0):
        data=squeeze(data)
        return data,dunit
    else:
        print 'error in read data ', ios
        print FILENAME
        return None, ""
def   READ_BBIO_DIURNALCYCLE(DATA_DIR, nx, ny, \
                            MONTH, DAY, YYYY, HOUR,\
                            CO2_DIR='CO2_200508//BBIO_DIURNALCYCLE', \
                            flnmpref='nep',  tracer=2,category='CO2-SRCE'):
    #  3/4 only 2x25 for now...
    ext2d=bp.get_name_ext_2d()
    extres=bp.get_res_ext()
    DOY=tm.day_of_year(YYYY, MONTH,DAY)
    
    # if (DOY>59 and mod(YYYY,4)==0): 
    #         DOY=DOY+1for ireg in  self.regs:
                    
            
    SDOY=r'%3.3d' % DOY
    
    FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+flnmpref.strip()+'.'+ext2d.strip()+'.'+extres.strip()+"."+SDOY
    # print 'read data from ', FILENAME
    
    tau0=tm.get_tau(1985,MONTH,DAY, HOUR)
    tau0=tau0/(3600.0)
    data, dunit, ios=bp.read_bpch2(FILENAME, category, tracer,tau0, \
                           nx,  ny, 1)
    
    if (ios==0):
        data=squeeze(data)
        return data,dunit
    else:
        print 'error in read data ', ios
        print FILENAME
        return None, ""

def   READ_MONTH_BIOBRN_CO2(DATA_DIR, nx, ny, \
                             MONTH, YYYY,\
                             LBBSEA=True, \
                             CO2_DIR='biomass_200110', \
                             flnmpref='bioburn',  tracer=4,\
                             category='BIOBSRCE'):
    EMFACTCO2CO = 12.068
    MONTHDATES = [31, 28, 31, 30,31, 30, 31, 31,30, 31, 30, 31]

    #  3/4 only 2x25 for now...
    ext2d=bp.get_name_ext_2d()
    extres=bp.get_res_ext()
    SYYYY=r'%4.4d' % YYYY
    
    TIME=MONTHDATES[MONTH-1] * 86400.0

    if (LBBSEA):
        FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+flnmpref.strip()+'.seasonal.'+ext2d.strip()+'.'+extres.strip()
        # print 'read data from ', FILENAME
        tau0=tm.get_tau(1985,MONTH,1)
        tau0=tau0/(3600.0)
    else:
        FILENAME = DATA_DIR.strip()+'/'+CO2_DIR+'/'+'bioburn'+'.interannual.'+ext2d.strip()+'.'+extres.strip()
        #        print 'read data from ', FILENAME
        tau0=tm.get_tau(YYYY,MONTH,1)
        tau0=tau0/(3600.0)

    data, dunit, ios=bp.read_bpch2(FILENAME, category, tracer,tau0, \
                                   nx,  ny, 1)
    if (ios==0):
        data=array(data)
        # Convert from [molec CO/cm2/month] to [molec CO2/cm2/month
        factor=EMFACTCO2CO/TIME
        data=factor*data
        data=squeeze(data)
        dunit='molec/cm2/s'
        return data,dunit
    else:
        print 'error in read data ', ios
        print FILENAME
        return None, ""
    
def EMISSCO2(YYYY, MONTH, DAY, HOUR=12, do_debug=False):
    lon, lat=GET_GRID()
    XNUMOL_CO2=gc.An/(1.0e-3*gc.mco2)
    
    x, y=meshgrid(lon, lat)
    x=transpose(x)
    y=transpose(y)
    region=DEFINE_FF_CO2_REGIONS(x,y)
    menus=gri.read_geos_input()
    em=menus['EMISSIONS MENU']
    sm=menus['SIMULATION MENU']
    em_details=gri.emission_menu_parse(em)
    sm_details=gri.simulation_menu_parse(sm)
    # print sm_details
    data_dir=sm_details['data_dir']
    nx=size(lon)
    ny=size(lat)
    idx=arange(1,nx+1)
    idy=arange(1,ny+1)
    sf_area=gm.grid_mod.get_areas(idx,idy)
    sf_area=1.0e4*sf_area # convert to cm2
    LBBSEA=em_details['do_biomass_em_sea']
    TS_CHEM=em_details['timestep']
    FF_REGION=DEFINE_FF_CO2_REGIONS(x,y)
    BB_REGION=DEFINE_FF_CO2_REGIONS(x,y)
    STT=zeros([nx, ny], float)
    
    # DTSRCE = 60.0 * TS_CHEM # convert to second
    DTSRCE=1.0
    
    factor=DTSRCE/XNUMOL_CO2
    iplot=1
    # bio mass burning 
    if (LBIOBRNCO2):
        
        data, dunit=READ_MONTH_BIOBRN_CO2(data_dir, nx, ny, \
                                                 MONTH, YYYY, LBBSEA)
        STT=STT+data
        #  Convert from [molec/cm2/s] to [kg]
        if (do_debug):
            subplot(3,2,iplot)
            gpl.plot_map(factor*sf_area*data, lon, lat, pcolor=1, title='Biomass Burning')
            iplot=iplot+1
    # biosphere 
    if ( LBIOCO2 ):
        if ( LUSECASANEP ):
            data, dunit=READ_BBIO_DIURNALCYCLE(datadir, nx, ny, \
                                               MONTH, DAY, YYYY, HOUR )
        else:
            data, dunit=READ_BBIO_DAILYAVERAGE(data_dir, nx, ny,\
                                               MONTH, DAY, YYYY)
        STT=STT+data
        if (do_debug):
            subplot(3,2,iplot)
            gpl.plot_map(factor*sf_area*data, lon, lat, pcolor=1, title='Biosphere')
            iplot=iplot+1

        
        
    #  Fossil fuel emissions
    if (LFOSSCO2):
        data, dunit=READ_ANNUAL_FOSSILCO2(data_dir, nx,ny)
         #  Convert from [molec/cm2/s] to [kg]
        STT=STT+data
        

        if (do_debug):
            subplot(3,2,iplot)
            gpl.plot_map(factor*sf_area*data, lon, lat, pcolor=1, title='Fossil')
            iplot=iplot+1

    # Oceanic exchange
    
    if ( LOCCO2):
        data,dunit=READ_ANNUAL_OCEANCO2(data_dir, nx,ny)
        STT=STT+data
        if (do_debug):
            subplot(3,2,iplot)
            gpl.plot_map(factor*sf_area*data, lon, lat, pcolor=1, title='Ocean')
            iplot=iplot+1

    #  Biofuel emissions
    if (LBIOFUELCO2):
        data, dunit=READ_ANNUAL_BIOFUELCO2(data_dir, nx,ny)
        STT=STT+data
        
        if (do_debug):
            subplot(3,2,iplot)
            gpl.plot_map(factor*sf_area*data, lon, lat, pcolor=1, title='Biofuel')
            iplot=iplot+1

    #   Net terrestrial exchange
    if (LBIONETCO2):
        data, dunit=READ_ANNUAL_BIONET_CO2(DATA_DIR, nx, ny)
        STT=STT+data
        
        if (do_debug):
            subplot(3,2,iplot)
            gpl.plot_map(data, lon, lat, pcolor=1, title='Net Terrestrial')
            iplot=iplot+1
        
        
    STT=factor*sf_area*STT
    if (do_debug):
        subplot(3,2,iplot)
        gpl.plot_map(STT, lon, lat, pcolor=1, title='Total (kg/hour)')
        iplot=iplot+1
        show()

    # for test
    # print 'stt', STT[1:3]
    # print shape(STT)
    
    # STT[1:3]=1.30*STT[1:3]
    # print 'stt-enlarged', STT[1:3]
    return STT
def STT_open_bpch2_write(funit,FLNM_PREF, title='CO2_FLUX'):
    ext1=bp.get_name_ext_2d()
    ext2=bp.get_res_ext()
    full_flnm=FLNM_PREF+"."+ext1+"."+ext2
    stat=bp.open_bpch2_for_write(funit,full_flnm, title)
    return stat
def STT_write_bpch2(funit,  ntracer, lonres, latres, stt, \
                    tau0, tau1,\
                    category="CO2_FLUX", dunit='kg/hour'):

    modelname=bp.get_modelname()
    centre180=1
    halfpolar=bp.get_halfpolar()
    ifirst, jfirst, lfirst=1, 1, 1
    
    dims=shape(stt)
    data=stt
    if (size(dims)==2):
        data=reshape(data, [dims[0], dims[1], 1])

    reserved="9999999"
    
    stat = bp.write_bpch2_data(funit,modelname,category,reserved, \
                               lonres,latres,halfpolar,centre180,\
                               ntracer,dunit,tau0,tau1,\
                               ifirst,jfirst,lfirst,data)
    return stat

def STT_close_bpch2(funit):
    stat=bp.close_bpch2_file(funit)
    return stat
def STT_read_bpch2(flnm, tracer, tau0, nx, ny, nz=1,  category='CO2_FLUX'):
    data, dunit, ios=bp.read_bpch2(flnm, category, tracer,tau0, \
                           nx,  ny, nz)
    if (ios==0):
        data=squeeze(data)
        return data, dunit
    else:
        return None, ""


class transcom_co2_st:
    def __init__(self, do_debug=False, sel_id=None, \
                 nlon_reg=4, nlat_reg=4):
        """ initial the class for co2 emission """
        # model grid  
        lon, lat, xedge, yedge, areas=GET_GRID(True)
        
        
        self.lon=array(lon)
        self.nx=size(lon)

        self.lat=array(lat)
        self.ny=size(lat)

        self.lonres=lon[2]-lon[1]
        self.latres=lat[2]-lat[1]
        self.scale_to_GtC=gc.mc*3600.0*365*24.0/(1.0e12*gc.mco2) # convert to GtC/cm2/y
        self.areas=areas # areas in cm2
        # print shape(areas)
        # print self.nx, self.ny
        if (do_debug):
            figure(2)
            subplot(2,1,1)
            gpl.plot_map(areas, lon, lat, use_pcolor=1)
        self.xedge=xedge
        self.yedge=yedge

        # emission grid (variables) 
        if (sel_id==None):
            reg_map, land_reg, ocean_reg, \
                     def_err=define_surface_co2_region(lon, lat, \
                                                   xedge, yedge, do_debug=False)
            
        else:
            reg_map, land_reg, ocean_reg,\
                     def_err=define_surface_co2_region_cell(lon, lat, \
                                                   xedge, yedge,sel_id, nlon_reg, nlat_reg, \
                                                          do_debug=False)
                
        self.reg_map=reg_map # needed to convert state vector to emission map 
        self.land_reg=land_reg # names of the region  0-11 land type 12-22 ocean type 
        self.ocean_reg=ocean_reg
        self.reg_name=list()
        
        for land in land_reg:
            self.reg_name.append(land)
        for ocean in ocean_reg:
            self.reg_name.append(ocean)

        self.nreg=len(self.land_reg)+len(self.ocean_reg)
        self.def_err=array(def_err) #
        self.def_err[1:]=(1.0e12*gc.mco2/(86400*365.0*gc.mc))*self.def_err[1:] # convert to kgCO2/second
        #
        
        # the time grid which contains how many 'time period' 
        self.ntime=0
        self.timestep=0
        self.tau0=list() # in seconds since 1985.01.01.00
        self.tau1=list()
        self.yyyy=list()
        self.doy=list()
        # the state vector and errors

        self.prior=None
        self.err=None
        self.st=None
        self.err_cov=None
        self.idx=list()
        self.reg_area=list()
        
        smap=zeros(shape(self.reg_map),float)
        
        for ireg in range(self.nreg):
            idx=where(self.reg_map==ireg)
            self.idx.append(idx)
            self.reg_area.append(sum(self.areas[idx]))
            smap[idx]=ireg
        
        if (do_debug):
            figure(1)
            subplot(2,1,1)
            gpl.plot_map(smap, self.lon, self.lat,use_pcolor=1)
            subplot(2,1,2)
            gpl.plot_map(self.reg_map, self.lon, self.lat,use_pcolor=1)
            show()

        
        self.reg_area=array(self.reg_area)
        print shape(self.reg_area), shape(self.def_err)
        print 'reg_area', self.reg_area[-4:-1]
        
        # self.def_err[1:]=self.def_err[1:]/self.reg_area[1:] # to kgCO2 /cm2/s
    
        
    def gen_prior(self, timestep, ntime, yyyy0, mm0, dd0, hh0=0, mi0=0, sec0=0, \
                 update_state=True, errset=None, do_debug=False):
        """

        construction a prior surface flux used for a time period starting with yyyy0, mm0, dd0

        """
        
        
        tau0=tm.get_tau(yyyy0, mm0, dd0, hh0, mi0, sec0)
        
        self.ntime=ntime
        self.timestep=timestep
        
        self.prior=zeros([self.nreg, self.ntime], float)
        self.err=zeros([self.nreg, self.ntime], float)
        
        
        for itime in range(self.ntime):
            tau1=tau0+timestep
            self.tau0.append(tau0)
            self.tau1.append(tau1)
            tau=0.5*(tau0+tau1)
            date_new=tm.tai85_to_utc(tau0)
            yyyy, mm, dd, hh, mi, sec=tm.utc_to_time_array(date_new)
            # print yyyy, mm, dd, hh
            doy=tm.day_of_year(yyyy, mm, dd)
            self.yyyy.append(yyyy)
            self.doy.append(doy)
                

            stt=EMISSCO2(yyyy,mm, dd, hh, do_debug=False)
            stv=self.map_to_st(stt)
            # for test only
            # if (itime in [0, 2, 4,6,8]):
            #    stv[3]=1.80*stv[3]
            #    stv[5]=1.80*stv[5]
            
                
            
            # stv=1.8*stv
            self.prior[:,itime]=stv
            
            
            if (itime==1 and do_debug):
                subplot(3,1,1)
                gpl.plot_map(stt, self.lon, self.lat, use_pcolor=1)
                subplot(3,1,2)
                stt2=self.st_to_map(stv)
                gpl.plot_map(stt2, self.lon, self.lat, use_pcolor=1)
                subplot(3,1,3)
                gpl.plot_map(self.reg_map, self.lon, self.lat, use_pcolor=1)
                show()
            
            if (errset==None):
                # err_fact=sqrt(365.0*24*3600.0/self.timestep)
                err_factor=5.0
                self.err[1:,itime]=err_factor*self.def_err[1:]
                self.err[0, itime]=self.def_err[0]*stv[0]
                # if (itime==0):
                #     print self.err[:,0]
                #    print self.prior[:,0]
                
            elif (size(errset)==1):
                self.err[:,itime]=errset*stv[:]
            else:
                self.err[:,itime]=errset[:]
            
            tau0=tau1
        self.construct_err_cov()
        
        if (update_state):
            self.st=array(self.prior)
        
            
         
            
    def map_to_st(self, stt):
        stv=zeros(self.nreg, float)
        for ireg in range(self.nreg):
            idx=self.idx[ireg]
            total_co2=stt[idx]
            total_co2=sum(total_co2)
            stv[ireg]=total_co2      #/self.reg_area[ireg] # to kgCO2/cm2/s

        return stv
    
    def st_to_map(self, stv, return_grid=False, rescale=True):
        lon=self.lon
        lat=self.lat
        stt=zeros([self.nx, self.ny], float)
        for ireg in range(self.nreg):
            idx=self.idx[ireg]
	    if (rescale):
            	stt[idx]=stv[ireg]*self.areas[idx]/self.reg_area[ireg]      
	    else:
		stt[idx]=stv[ireg]
            # print ireg, stv[ireg]
        
        if (return_grid):
            return stt, lon, lat
        else:
            return stt
    
            
    def save_ensemble_bpch2(self, itime, devs,\
                            ids=1, \
                            flnm='CO2_EMISSION_EN', \
                            category='CO2_FLUX',\
                            scale_n=1.0,\
                            do_debug=False):
        
        """

        save the ensembles at one time step to one file
        
        """

        nmem=len(devs)
        memb=zeros([self.nx, self.ny, nmem], float)
        stv=self.prior[:,itime]
        nvar=size(stv)
        stvm=zeros([nvar, nmem], float)
        for im in range(nmem):
            dev=devs[im]
            dev=dev/scale_n
            new_stv=stv+dev
            stt=self.st_to_map(new_stv)
            memb[:,:,im]=stt[:,:]
            stvm[:,im]=new_stv[:]
            # if (im==0 or im==3 or im==4):
            #     print 'stv--', im, self.scale_to_GtC*new_stv[0:2], self.scale_to_GtC*dev[0:2]
        if (do_debug and itime==0):
            subplot(2,1,1)
            vals=memb[:,:,0]
            vals=self.scale_to_GtC*vals
            gpl.plot_map(vals, self.lon, self.lat, use_pcolor=1, title='Emission (priori) (Gt C/cell/y)')
            subplot(2,1,2)
            vals=self.st_to_map(self.err[:,itime])
            # vals=memb[:,:,2*10]
            vals=self.scale_to_GtC*vals
            gpl.plot_map(vals, self.lon, self.lat, use_pcolor=1, title='Emission error (Gt C/cell/y)')
            # savefig('emission_err.pdf')
            show()
        if (flnm==None):
            return memb
        else:
            ext1=bp.get_name_ext_2d()
            ext2=bp.get_res_ext()
            doy=self.doy[itime]
            yyyy=self.yyyy[itime]
            sdate=r'%4.4dD%3.3d' %(yyyy, doy)
            full_flnm=flnm+sdate+"."+ext1+"."+ext2
            funit=38
            title='ENKF_CO2'
            stat=bp.open_bpch2_for_write(funit,full_flnm, title)
            if (stat<>0):
                return stat
            
            lonres=self.lonres
            latres=self.latres
            tau0=self.tau0[itime]/3600.0  # to hours
            tau1=self.tau1[itime]/3600.0
            dunit='kg/s',
            modelname=bp.get_modelname()
            centre180=1
            halfpolar=bp.get_halfpolar()
            ifirst, jfirst, lfirst=1, 1, 1
            dims=shape(stt)
            reserved="9999999"
            ntracer=ids
            
            stat = bp.write_bpch2_data(funit,modelname,category,reserved, \
                                       lonres,latres,halfpolar,centre180,\
                                       ntracer,dunit,tau0,tau1,\
                                       ifirst,jfirst,lfirst,memb)

            category=category+'_reg'
            stat = bp.write_bpch2_data(funit,modelname,category,reserved, \
                                       lonres,latres,halfpolar,centre180,\
                                       ntracer,dunit,tau0,tau1,\
                                       ifirst,jfirst,lfirst,stvm)
            
            bp.close_bpch2_file(funit)
            if (ids==1):
                return (ids+nmem)
            else:
                return (ids+nmem-1)
            
         
    
    def save_ensemble_bin(self, itime, devs, ids=1, flnm='CO2_EMISSION_EN', category='CO2_FLUX', scale_n=sqrt(2.0)):

        """

        save the ensembles at one time step to one file
        
        """
        nmem=len(devs)
        memb=zeros([self.nx, self.ny, nmem], float)
        stv=self.prior[:,itime]
        
        for im in range(nmem):
            new_stv=stv+devs[im]
            
            
            stt=self.st_to_map(new_stv)
            memb[:,:,im]=stt[:,:]
        
        if (flnm==None):
            return memb
        else:
            ext1=bp.get_name_ext_2d()
            ext2=bp.get_res_ext()
            doy=self.doy[itime]
            yyyy=self.yyyy[itime]
            sdate=r'%4.4dD%3.3d' %(yyyy, doy)
            full_flnm=flnm+sdate+"."+ext1+"."+ext2
            funit=38
            title='ENKF_CO2'
            stat=bp.open_bpch2_for_write(funit,full_flnm, title)
            if (stat<>0):
                return stat
            
            ntracer=ids
            lonres=self.lonres
            latres=self.latres
            tau0=self.tau0[itime]/3600.0  # to hours
            tau1=self.tau1[itime]/3600.0
            dunit='kg/s',
            modelname=bp.get_modelname()
            centre180=1
            halfpolar=bp.get_halfpolar()
            ifirst, jfirst, lfirst=1, 1, 1
            dims=shape(stt)
            reserved="9999999"
            stat = bp.write_bpch2_data(funit,modelname,category,reserved, \
                                       lonres,latres,halfpolar,centre180,\
                                       ntracer,dunit,tau0,tau1,\
                                       ifirst,jfirst,lfirst,memb)
            bp.close_bpch2_file(funit)
            
        return stat
    
    def gen_def_ensemble(self, itime, en_st=0, en_end=None, do_sel=0):
        """
        generate a simple ensemble
        """
        if (en_end==None):
            en_end=self.nreg
        ensemble=list()
        dev=zeros(self.nreg,float)
        ensemble.append(dev)
        if (do_sel==0): # simply choose the every one  
            for ien in range(en_st, en_end):
                dev=zeros(self.nreg,float)
                dev[ien]=self.err[ien, itime]
                ensemble.append(array(dev))  # + err
                dev[ien]=-self.err[ien, itime]
                ensemble.append(array(dev)) # -err
            return ensemble
        elif (do_sel==1): # do SVD
            err_cov=self.err_cov[:,:,itime]
            u,w, v=nlg.svd(err_cov)
            u=u*sqrt(w)
            
            # we will use u
            for ien in range(en_st, en_end):
                dev=u[:,ien]
                ensemble.append(dev)
                dev=-u[:, ien]
                ensemble.append(dev)
            return ensemble
        
        elif (do_sel==2): # LU decomposition
            u=nlg.cholesky(err_cov)
            # we will use u
            for ien in range(en_st, en_end):
                dev=u[:,ien]
                ensemble.append(dev)
                dev=-u[:, ien]
                ensemble.append(dev)
            return ensemble
        
    def construct_err_cov(self, cor_land=900.,  cor_ocean=2000, do_debug=False):
        """ generat the error correlation between regions """
        #  we get the distance
        lon_m, lat_m=meshgrid(self.lon, self.lat)
        lon_m=transpose(lon_m)
        lat_m=transpose(lat_m)

        
        self.err_cov=zeros([self.nreg, self.nreg, self.ntime], float)
        ist=0
        print 'lang_reg', size(self.land_reg), size(self.ocean_reg)
        for i in range(len(self.land_reg)):
            px=self.idx[i]
            lonx=lon_m[px]
            latx=lat_m[px]
            wx=self.areas[px]
            wx=wx/self.reg_area[i]
            self.err_cov[i,i,:]=self.err[i,:]*self.err[i,:]
            
            
            for j in range(i+1, len(self.land_reg)):
                py=self.idx[j]
                lony=lon_m[py]
                laty=lat_m[py]
                
                dist=gcd.get_circle_distance_xy(lonx, latx, lony, laty)
                
                
                wy=self.areas[py]
                wy=wy/self.reg_area[j]
                
                schur_factor=exp(-dist/cor_land)
                # print shape(schur_factor)
                schur_factor=dot(schur_factor, wy)
                schur_factor=dot(wx, schur_factor)
                
                for itime in range(self.ntime):
                    errx=self.err[i, itime]
                    erry=self.err[j, itime]
                    self.err_cov[i, j, itime]=errx*erry*schur_factor
                    self.err_cov[j, i, itime]=self.err_cov[i,j, itime]

            ist=ist+1

        for i in range(ist, ist+len(self.ocean_reg)):
            px=self.idx[i]
            lonx=lon_m[px]
            latx=lat_m[px]
            wx=self.areas[px]
            wx=wx/self.reg_area[i]
            self.err_cov[i,i,:]=self.err[i,:]*self.err[i,:]
            
            
            for j in range(i+1, ist+len(self.ocean_reg)):
                py=self.idx[j]
                lony=lon_m[py]
                laty=lat_m[py]
                dist=gcd.get_circle_distance_xy(lonx, latx, lony, laty)
                
                wy=self.areas[py]
                wy=wy/self.reg_area[j]

                schur_factor=exp(-dist/cor_ocean)
                schur_factor=dot(schur_factor, wy)
                schur_factor=dot(wx, schur_factor)
            
                
                for itime in range(self.ntime):
                    
                    errx=self.err[i, itime]
                    erry=self.err[j, itime]
                    self.err_cov[i, j, itime]=errx*erry*schur_factor
                    self.err_cov[j, i, itime]=self.err_cov[i,j, itime]
    

        idx=[23,24,25]
        print 'self_err', diag(self.err_cov[20:, 20:, 0])
        
                        
        
    def print_reg_annual_emission(self, itime=0):

        stv=self.prior[:,  itime]
        err=self.err[:,itime]
        head_width=30
        print '='*80
        sdate=r'%4.4dD%3.3d' % (self.yyyy[itime],self.doy[itime])
        print '----------CO2 EMISSION ON '+sdate+'-------------'
        print '      region name:   emission (Gt C/y)    err (Gt C/y)'
        print '='*80
        temi=0.0
        terr=0.0
        for ireg in range(self.nreg):
            reg_name=self.reg_name[ireg]
            reg_name=reg_name.strip()
            lrg=len(reg_name)
            line_head=reg_name+' '*(head_width-lrg)
            emi=self.scale_to_GtC*stv[ireg] # self.reg_area[ireg]
            emi_err=self.scale_to_GtC*err[ireg] # *self.reg_area[ireg]
            temi=temi+emi
            terr=terr+emi_err**2
            sline=r'%13.5f, %13.5f' % (emi, emi_err)
            line=line_head+': '+sline
            print line
        print '-'*80
        reg_name='total'
        lrg=len(reg_name)
        line_head=reg_name+' '*(head_width-lrg)
        sline=r'%13.5f, %13.5f' % (temi, sqrt(terr))
        line=line_head+': '+sline
        print line
        print '-'*80
        
        print '-'*30+'Error covariance'+'-'*30
        err_mat=self.scale_to_GtC*self.scale_to_GtC*self.err_cov[:,:,itime]
        print 'I am here', diag(err_mat[20:, 20:])
        print  err_mat[10:, 23]
        
        subplot(2,1,1)
        pcolor(err_mat[3:,3:], vmin=0.0, vmax=6.5)
        # xlim([0,22])
        # ylim([0,22])
        title('Error covariance (GtC/y)^2')
        colorbar()
        subplot(2,1,2)
        u,w,v=nlg.svd(err_mat)
        plot(w)
        xlim([0,22])
        title('Singular values')
        # savefig('err_cov.ps')
        show()
            
        
        # print array2string(err_mat, precision=3)
        
        print '='*80
        
                          
            
        


if (__name__=="__main__"):
    co2=transcom_co2_st(sel_id=3)
    # starting time 
    yyyy=2003
    mm=3
    dd=18
    timestep=8.0*24.0*3600.0 #  8-days
    ntime=12
    pos=list()
    ipos=0
    ftt=open("ens_pos.dat", "w")
    
    co2.gen_prior(timestep, ntime, yyyy,mm, dd)
    co2.print_reg_annual_emission()
    
    ids=1
    for i in range(ntime):
        
        devs=co2.gen_def_ensemble(i)
        nmem=len(devs)
        co2.save_ensemble_bpch2(i, devs)
        ipos=ipos+nmem
        pos.append(ipos)
        line=r'%4.4d %4.4d %4.4d %3.3d %12.5f %12.5f' % (ipos, nmem, co2.yyyy[i], co2.doy[i], co2.tau0[i]/3600.0, co2.tau1[i]/3600)
        ftt.write(line+'\n')
        ids=ids+nmem
    
    ftt.close()

    
        
    






























    
              
    
        
