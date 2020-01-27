###! /home/lfeng/tmp_cp/ecmwf/python2.5/bin/python2.5

"""
        Originally written by Dr. Liang Feng (University of Edinburgh), circa 2009
        Orginal environment: python 2.5
        Rewritten and commented by Mehliyar Sadiq, 2020
        First working version: 2020-xxxx
"""


# default modules in Python:
import sys # constants, functions and methods of the Python interpreter
import os # operating system dependent functionality
import time as systime # functions for working with time

# modules written by Liang:
# definitions of paths, time, resolution, inversion options, instruments, etc.
import geos_chem_def as gcdf

#import input_geos_gen as igg # create the new input to drive the ensemble run
#import co2_emission as co2em # ???
#import restart_gen as rg # ???
import time_module as tm # time conversion


### START ###
print('*'*70)
print('*'*20+'CO2 ENSEMBLE RUN DRIVER'+'*'*30)
print('*'*70)

print('======>Step 1: Generate co2 emission data<======')

#co2=co2em.transcom_co2_st()

# cound the following definitions be moved to geos_chem_def.py??? 
# starting time
yyyy, mm, dd=2003,1,1

temp_res=8                     # 8 days 
timestep=temp_res*24.0*3600.0  # 8-days in seconds
ntime=12                       # 12 cycles? 12 ensembles?
pos=list()
ipos=0                         # ???
gmt=systime.gmtime()           # GMT time. UTC time
a_mst=[1, 186, 370]            # ??? model start?
a_mend=[185, 369, 553]         # ??? model end?
# a_mst=[1]                    # testing?
# a_mend=[185]
new_restart=True
nrun=len(a_mst)                # number of runs?

ftt=open(gcdf.data_path+"ens_pos.dat", "w") # is ens_pos_dat a temp file to record information?
line='geos_chem run at %4.4d%2.2d%2.2d, %2.2d:%2.2d:%2.2d' % (gmt[0], gmt[1], gmt[2], gmt[3], gmt[4], gmt[5]) # echos real time, as starting point
print(line)
ftt.write(line+'\n')
line=r'temp_res: %4.4d  nstep: %4.4d' % (temp_res, ntime)
ftt.write(line+'\n')
line='mem_st mem_end year_st  year_end day_st day_end flnm'
ftt.write(line+'\n')

for irun in range(nrun):
    mst=a_mst[irun]
    mend=a_mend[irun]
    print('-'*10+'starting year, month, day:',  yyyy, mm, dd)
    print('-'*10+ 'timestep (day), ntime',  timestep/(24.0*3600.0),  ntime)
    if (irun==0):
        print('-'*10+'generat  prior' +'-'*10)
        co2.gen_prior(timestep, ntime, yyyy,mm, dd)
        ids=1 # the first one is the mean state
        # the prefix  for co2 emission 
        co2flnm=gcdf.data_path+'/'+'CO2_EMISSION_EN' 
        print('-'*10+'generate ensemble member' +'-'*10)
        for i in range(ntime):
            devs=co2.gen_def_ensemble(i)
            nmem=len(devs)
            co2.save_ensemble_bpch2(i, devs, ids, flnm=co2flnm)
            ipos=ipos+nmem
            pos.append(ipos)
            if (ids==1):
                ids=ids+nmem
            else:
                ids=ids+nmem-1


    line=r'%4.4d %4.4d %4.4d %3.3d %3.3d %3.3d' % (mst, mend, co2.yyyy[0], co2.yyyy[-1], co2.doy[0], co2.doy[-1])
    line=line+ ' '+co2flnm
    print(line)
    ftt.write(line+'\n')
    
    print('======>Step 2: Generate input file<======')

    # igg.create_new_input_file(co2.yyyy, co2.doy, member_start=mst, \
    #                      member_end=mend, co2flnm=co2flnm, time_start='20030328')
    igg.create_new_input_file(co2.yyyy, co2.doy, member_start=mst, \
                              member_end=mend, co2flnm=co2flnm)
    
    os.remove('input.geos')
    os.rename('input.geos.new', 'input.geos')
    
    
    print('======>Step 3: Generate restart file<=====')
    
    ntracers=mend-mst+1
    rsf=rg.geos_chem_restart_file('restart.jan2003.kalman.borealasia')
    yst=co2.yyyy[0]
    dst=co2.doy[0]
    tau0=co2.tau0[0]/3600.0
    tst=tm.doy_to_utc(dst, sec=0, yyyy=yst)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '')
    enaf=r'EN%4.4d-EN%4.4d' % (mst, mend)
    full_restart_name='restart.'+enaf+'.'+tst[0:8]+'00'
    real_ntracers=ntracers
    
    if (new_restart):
        new_lon, new_lat=co2em.GET_GRID()
        rsf.copy_restart_file(1, full_restart_name,real_ntracers, tau0, do_regrid=True, new_lon=new_lon, new_lat=new_lat)
    
    print('======>Step 4: Launch geos-chem<======')
    os.system('sh ./rungeos.sh')
ftt.close()
        


