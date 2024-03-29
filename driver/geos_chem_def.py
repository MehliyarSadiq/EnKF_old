# GEOS-Chem definitions
import time_module as tm

# GEOS-Chem related directories
run_path='/scratch/local/lfeng/oco2_project/enkf_oco2_inv_792/' # 
data_path=run_path+'/enkf_output/'     # model output
rerun_datapath=run_path+'/enkf_rerun/' # restart directory
hm_data_path=run_path+'/surface_hm/'   # model output
sat_hm_data_path=run_path+'/b340_hm/'  # model output

# observation related
# obs_path=run_path+'/oco_obs/'            
inv_path=run_path+'/oco_inv/'           # inversion directory
oco_orbit_path=run_path+'/gosat_orbit/' # OCO orbit 
oco_ak_path=run_path+'/gosat_ak/'       # OCO averaging kernel

# for multiple observation
view_type='new_iss' # ???
view_mode= "glint"  # 'glint', 'nadir' or 'n16g16'
obs_path=run_path+'/'+view_type+'_obs/' # observation directory

def_input_file=run_path+'/input.geos'   # GEOS-Chem model input file

#  starting time for geos chem simulation
st_yyyy, st_mm, st_dd=2009,1,1
st_doy=tm.day_of_year(st_yyyy, st_mm, st_dd) # day of year

# time resultion in day
use_fixed_mod_err=True  # ???
temporal_resolution=32  # temporal resolution for inversions in days

# ??? the ending time for geos chem simulation
inv_lag_window=4
# inversion lag window = temporal_resolution*n_inv_window

ntime_geos_chem=inv_lag_window 
# 2*inv_lag_window if used the moving window # the last day of geos-chem simulation is given by ntime * temporal_resolution
inv_ntime=inv_lag_window

# spatial resolution for land and ocean
n_land_lat_div, n_land_lon_div=3,3
n_ocean_lat_div,n_ocean_lon_div=2,2 # or 3,3

# number of regions???
region_num=11*n_land_lat_div*n_land_lon_div+11*n_ocean_lat_div*n_ocean_lon_div+1
num_used_en=max(340, region_num) # number of ensembles

output_daily_obs=True
new_restart=True

restart_file=data_path+'restart.jan2003.kalman.borealasia'  # first restart file/initial condition

do_short= False # if choose do short will start a run with only 5 tagged tracers


######  the following section  is  for inversion options #####
do_update=True    # force the system to read in the data again

xy_obs_file=data_path+'/'+view_mode+'_.2003D001_N02.nc' # if use old data provide the file name

select_obs=False # use all observation 
add_rnd_err= False # add random error
sate_err_scale=1.0
obs_err_scale=1.0    

# model errors
model_land_err=1.0 # 1.5 # ppm
model_ocean_err=1.0  # 1.5 # ppm
model_err=1.0  #.5 # ppm
cor_len=200.0
tcor_len=1.0
tcor_factor=1.0


# tracer info file needs to initialize
def_tracerinfo_file=run_path+"tracerinfo.dat" # ???
def_diaginfo_file=run_path+"diaginfo.dat"     # ???

# a-priori error
# rescale the error covariance.
# if necessary, you may use it to
# convert the default digonal error covariance to a more complicated one.
xnorm=0.3


# run control ???
do_debug=False    
do_retry=False
do_init_run=True
maxne=inv_lag_window*region_num
start_step=0


# tell which instrument will be ???????

# view_type_list=['new_aqua', 'new_landsat']
# view_mode_list=["glint", "glint"]
# view_station_list=['satellite', 'satellite']

view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=list()
for view_type in view_type_list:
    opath='./'+view_type+'_obs/'
    obs_path_list.append(opath)

# obs_path_list=['./new_trmm_obs/']
hm_update_list=[True]
mean_update_list=[True]
obs_use_simulation=False
obs_simulation_path='./obs_aqua_simulation/'
add_mean_flux_pb=0.0 #1.6


# view_type_list=['new_gosat','new_gosat']
# view_mode_list=["nadir", "glint" ]
# view_station_list=['satellite', 'satellite']
# obs_path_list=list()
# for view_type in view_type_list:
#     opath='./'+view_type+'_obs/'
#     obs_path_list.append(opath)

# obs_path_list=['./gosat_obs/']
# hm_update_list=[True, True]
# mean_update_list=[True, True]
# obs_use_simulation=False
# obs_simulation_path='./obs_aqua_simulation/'

# view_type_list=['esrl_flask',  'esrl_in_situ', 'esrl_tower']
# view_mode_list=["", "", ""]
# view_station_list=['surface', 'surface', 'surface']
# obs_path_list=['./esrl_obs/flask/','./esrl_obs/in-situ/', './esrl_obs/tower/']
# hm_update_list=[False, False, False]
# mean_update_list=[True, True, True] # ,True,True,True]

view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./obs_2014/event_nc/']
hm_update_list=[False]
mean_update_list=[True]


view_type_list=['TCC']
view_mode_list=[""]
view_station_list=['tccon']
obs_path_list=['./tccon/']
hm_update_list=[False]
mean_update_list=[True]


view_type_list=['esrl_flask',  'TCC']
view_mode_list=["", ""]
view_station_list=['surface', 'tccon']
obs_path_list=['./obs_2014/event_nc/', './tccon/']
hm_update_list=[False, False]
mean_update_list=[True, True]

view_type_list=['esrl_flask', 'TCC']
view_mode_list=["", ""]
view_station_list=['surface', 'tccon']
obs_path_list=['./obs_2014/event_nc/','./tccon_2014/'] # different
hm_update_list=[False, False]
mean_update_list=[True, True]

view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./obs_2014/event_nc/']
hm_update_list=[True] # different
mean_update_list=[True]

view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=[ './new_b340_r3']
hm_path_list=['./b340_hm/']
hm_update_list=[False]
mean_update_list=[True]


view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=[ './gosat_obs_v5/'] # different
hm_path_list=['./uol5_hm/']        # different
hm_update_list=[False]
mean_update_list=[True]


view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=[ './gosat_uol_mod2/'] # different
hm_path_list=['./uol4_hm/']
hm_update_list=[False]
mean_update_list=[True]


view_type_list=['gosat_v29', 'esrl_flask']
view_mode_list=["", ""]
view_station_list=['satellite', 'surface']
obs_path_list=[ './oco2_10s/', './new_insitu_obs/']
hm_path_list=['./oco2_hm/', './surface_hm/']
hm_update_list=[False, False]
mean_update_list=[True, True]

view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./new_insitu_obs/']
hm_update_list=[False]
mean_update_list=[True]
hm_path_list=['./surface_hm/']

view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=[ './oco2_10s/']
hm_path_list=['./oco2_hm/']
hm_update_list=[False]
mean_update_list=[True]

view_type_list=['esrl_flask',  'esrl_tower', 'gosat_v29']
view_mode_list=["", "", ""]
view_station_list=['surface', 'surface', 'satellite']
obs_path_list=['./obs_2016/flask/','./obs_2016/tower/', './gosat_obs_v6/']
hm_update_list=[False, False, False]
mean_update_list=[True, True, True]
hm_path_list=['./surface_hm/', './surface_hm/', './uol_hm_v6/']

view_type_list=['esrl_flask', 'esrl_tower', 'gosat_v29']
view_mode_list=["", "", ""]
view_station_list=['surface', 'surface','satellite']
obs_path_list=['./obs_2016/flask/', './obs_2016/tower/', './oco2_10s_v10/']
hm_path_list=['./surface_hm/', './surface_hm/', './oco2_hm/']
hm_update_list=[True, True, True]
mean_update_list=[True, True, True]

view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=[ './oco2_10s_v10/']
hm_path_list=['./oco2_hm/']
hm_update_list=[True]
mean_update_list=[True]

view_type_list=['esrl_flask',  'gosat_v29']
view_mode_list=["", "", ""]
view_station_list=['surface', 'satellite']
obs_path_list=['./obs_2016/flask/','./oco2_10s_v10/']
hm_update_list=[False, False]
mean_update_list=[True, True]
hm_path_list=['./surface_hm/', './oco2_hm/']

view_type_list=['esrl_flask',  'flask_prep']
view_mode_list=["", ""]
view_station_list=['surface', 'flask_prep']
obs_path_list=['./obs_2016/flask/','./surface_prep_obs/']
hm_update_list=[True, True]
mean_update_list=[True, True]
hm_path_list=['./surface_hm/', './surface_prep_hm/']

view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./obs_2016/flask/']
hm_update_list=[False]
mean_update_list=[True]
hm_path_list=['./surface_hm/']

view_type_list=['esrl_flask', 'esrl_tower']
view_mode_list=["", ""]
view_station_list=['surface', 'surface']
obs_path_list=['./obs_2016/flask/', './obs_2016/tower/']
hm_update_list=[False, False]
mean_update_list=[True, True]
hm_path_list=['./surface_hm/', './surface_hm/']

view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=[ './gosat_obs_v7.1/']
hm_path_list=['./uol_hm_v7.1/']
hm_update_list=[False]
mean_update_list=[True]

view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./obs_2016/flask/']
hm_update_list=[True]
mean_update_list=[True]
hm_path_list=['./surface_hm/']

instrument_id={'gv':9001,\
               'ersl_flask':1001, \
               'ersl_in_situ':1002, \
               'ersl_tower':1003, \
               'ce_stn':1004, \
               'new_gosat_nadir':2001, 'new_gosat_glint':2002, \
               'new_aqua_nadir':3001,  'new_aqua_glint':3002, \
               'new_trmm_nadir':4001, 'new_trmm_glint':4002,\
               'new_landsat_nadir':5001, 'new_landsat_glint':5002}

add_ersl=False

nbias=0        # ??????
mod_bias=[0.0]
bias_err=[0.0]

