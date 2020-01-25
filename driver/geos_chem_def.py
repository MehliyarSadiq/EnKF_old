# default directories
import time_module as tm
run_path='/scratch/local/lfeng/oco2_project/enkf_oco2_inv_792/'
data_path=run_path+'/enkf_output/' # the model output
rerun_datapath=run_path+'./enkf_rerun/'
hm_data_path=run_path+'/surface_hm/' # the model output
sat_hm_data_path=run_path+'/b340_hm/' # the model output


# rerun_datapath=run_path+'./enkf_rerun_pb/'

# obs_path=run_path+'/oco_obs/'            # observations
inv_path=run_path+'/oco_inv/'
def_input_file=run_path+'/input.geos' # GEOS-CHEM model input file
oco_orbit_path=run_path+'/gosat_orbit/' # OCO orbit 
oco_ak_path=run_path+'/gosat_ak/'       # OCO averaging kernel
#  starting time for geos chem simulation
view_type='new_iss'
view_mode= "glint"
# for multipe observation



obs_path=run_path+'/'+view_type+'_obs/'

st_yyyy, st_mm, st_dd=2009,1,1
st_doy=tm.day_of_year(st_yyyy, st_mm, st_dd)

# time resultion in day
use_fixed_mod_err=True
temporal_resolution=32  #  the temporal resolution for inversions in days

# the ending time for geos chem simulation
# inversion lag window =temporal_resolution*n_inv_window
inv_lag_window=4 #

ntime_geos_chem=inv_lag_window # 2 *inv_lag_window if used the moving window # t#the last day of geos-chem simulation is given by ntime * temporal_resolution
inv_ntime=inv_lag_window

# spatial resolution

n_land_lat_div, n_land_lon_div=3,3
# n_ocean_lat_div,n_ocean_lon_div=3,3
n_ocean_lat_div,n_ocean_lon_div=2,2

region_num=11*n_land_lat_div*n_land_lon_div+11*n_ocean_lat_div*n_ocean_lon_div+1
num_used_en=max(340, region_num) # the restart file
output_daily_obs=True

new_restart=True
restart_file=data_path+'restart.jan2003.kalman.borealasia'  # 

#
# if choose do short will start a run with only 5 tagged tracers

do_short= False
# view modes 
# view_mode= 'glint'  # 'glint' # 'nadir', 'n16g16'

######  the following section  is  for inversion options ####
# observation options
do_update=True    # force the system to read in the data again
# if use old data provide the file name
xy_obs_file=data_path+'/'+view_mode+'_.2003D001_N02.nc'

select_obs=False # True # False  # use all observation 
add_rnd_err= False # True # True # False  # True # True # True
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
def_tracerinfo_file=run_path+"tracerinfo.dat"
def_diaginfo_file=run_path+"diaginfo.dat"

#  a-priori error
#  rescale the error covariance.
# if necessary, you may use it to
# convert the default digonal error covariance to a more complicated one.
xnorm=0.3





# run control 
do_debug=False    
do_retry=False
do_init_run=True
maxne=inv_lag_window*region_num
start_step=0


# tell which instrument will be 

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



view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./obs_2014/event_nc/']
hm_update_list=[False]
mean_update_list=[True]

view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./obs_2014/event_nc/']
hm_update_list=[False]
mean_update_list=[True]



view_type_list=['esrl_flask', 'TCC']
view_mode_list=["", ""]
view_station_list=['surface', 'tccon']
obs_path_list=['./obs_2014/event_nc/','./tccon_2014/']
hm_update_list=[False, False]
mean_update_list=[True, True]

view_type_list=['esrl_flask']
view_mode_list=[""]
view_station_list=['surface']
obs_path_list=['./obs_2014/event_nc/']
hm_update_list=[True]
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
obs_path_list=[ './gosat_obs_v5/']
hm_path_list=['./uol5_hm/']
hm_update_list=[False]
mean_update_list=[True]


view_type_list=['gosat_v29']
view_mode_list=[""]
view_station_list=['satellite']
obs_path_list=[ './gosat_uol_mod2/']
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

nbias=0
mod_bias=[0.0]
bias_err=[0.0]









