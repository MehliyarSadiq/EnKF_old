import  bpch2_rw_v2 as brw
import bpch2_rw_py as bp
from pylab import *
import gen_plots as gpl
from numpy import *
import time_module as tm
class geos_chem_restart_file:
    def __init__(self, org_flnm=None):
        self.ntracers=0
        if (org_flnm!=None):
             bpch2=brw.bpch2_file_rw(org_flnm, 'r', do_read=1)
             print(bpch2.stat)
             self.org_flnm=org_flnm
             self.bpch2=bpch2
        else:
             self.org_flnm=None
             self.bpch2=None
        
    def copy_restart_file(self,sel_idx, newfile, ntracers=1,newtau=None, do_regrid=False, new_lon=None, new_lat=None):
        """ copy to new_restart file and shift the time to a new time"""

        data_list, found=self.bpch2.get_data(None, sel_idx, None, None)
        if (found[0]==0):
            print('No record found')
            return None
        
        bpdata=data_list[0]
        bpch2_w=brw.bpch2_file_rw(newfile, 'w')
        ifunit=95
        bpch2_w.open_bpch2_w(ifunit,'test')
        
        if (newtau!=None):
            tau1, tau0=bpdata.get_attr(['tau0', 'tau1'])
            tau1=newtau+tau1-tau0
            tau0=newtau
            bpdata.set_attr('tau0', tau0)
            bpdata.set_attr('tau', tau1)
        sl=r'write restart data for %d tracers into file ' % ntracers
        print('='*80)
        print(sl+newfile.strip())
        if (do_regrid):
            print('I am do regrid', size(new_lon), size(new_lat))
            
            bpdata.regrid_data(new_lon, new_lat)
        
        for i in range(ntracers):
            bpdata.ntracer=i+1
            traname=bpdata.get_attr(['name'])
            print(i, traname[0].strip(), bpdata.category.strip(), bpdata.ntracer, bpdata.unit.strip(), shape(bpdata.data), max(bpdata.data.flat), min(bpdata.data.flat))
            
            
            bpdata.write_to_bpch2_file(ifunit)
            
        bpch2_w.close_bpch2_file()
        print('-'*30+'Done'+'-'*30)
        print('='*80)

    def mod_restart_file(self,sel_idx, newfile, ntracers=1, \
                         newtau=None, fixed_value=None, do_regrid=False,  \
                         new_lon=None, new_lat=None):
        """ copy to new_restart file and shift the time to a new time"""

        data_list, found=self.bpch2.get_data(None, sel_idx, None, None)
        if (found[0]==0):
            print('No record found')
            return None
        
        bpdata=data_list[0]
        bpch2_w=brw.bpch2_file_rw(newfile, 'w')
        ifunit=95
        bpch2_w.open_bpch2_w(ifunit,'test')
        
        if (newtau!=None):
            tau1, tau0=bpdata.get_attr(['tau0', 'tau1'])
            tau1=newtau+tau1-tau0
            tau0=newtau
            bpdata.set_attr('tau0', tau0)
            bpdata.set_attr('tau', tau1)
        sl=r'write restart data for %d tracers into file ' % ntracers
        print('='*80)
        print(sl+newfile.strip()
        if (do_regrid):
            print('I am doing regrid', size(new_lon), size(new_lat))
            
            bpdata.regrid_data(new_lon, new_lat)

        
        bpdata.ntracer=1
        bpdata.write_to_bpch2_file(ifunit)
        traname=bpdata.get_attr(['name'])
        print(0, traname[0].strip(), bpdata.category.strip(), bpdata.ntracer, bpdata.unit.strip(), shape(bpdata.data), max(bpdata.data.flat), min(bpdata.data.flat))
        if (fixed_value!=None):
            bpdata.data=bpdata.data-bpdata.data+fixed_value
        
        for i in range(1, ntracers):
            bpdata.ntracer=i+1
            traname=bpdata.get_attr(['name'])
            print(i+1, traname[0].strip(), bpdata.category.strip(), bpdata.ntracer, bpdata.unit.strip(), shape(bpdata.data), max(bpdata.data.flat), min(bpdata.data.flat))
            
            
            bpdata.write_to_bpch2_file(ifunit)
            
        bpch2_w.close_bpch2_file()
        print('-'*30+'Done'+'-'*30)
        print('='*80)
        
    
    def adjust_restart_file(self, restart_flnms, enst, enend, nt_time, tma):
        """ adjust restart file after read in """
        data_list=list()
        id_list=list()
        ifile=0
        file_id=list()
        for flnm in restart_flnms:
            bpch2=brw.bpch2_file_rw(flnm, 'r', do_read=1)
            en_start=enst[ifile]
            en_end=enend[ifile]
            en_id_time=en_start
            id=en_id_time
            file_ids=list()
            data_array=list()
            for idata in bpch2.data:
                data_array.append(idata.data)
                data_list.append(idata)
                id_list.append(id)
                id=id+1
                if (id-en_id_time==en_end):
                    en_id_time=en_id+nt_time
                    id=en_start+nt_time

        sort_id=argsort(id_list)
        data_list=array(data_array)
        nx, ny, nen=shape(data_array)
        data_list=reshape(data_array, [1, nx*ny, nen])
        data_list=squeeze(data_array)
        sorted_data_list[:,] =data_list[:, sort_id]
        sorted_data_list=dot(sorted_data_list, tma)
        data_list[:,sort_id]=data_list[:,:]
        
        for ifile in range(len(restart_files)):
            idx=file_id[ifile]
            bpch2_w=brw.bpch2_file_rw(restart_files[ifile], 'w')
            ifunit=95
            bpch2_w.open_bpch2_w(ifunit,'test')
           
            for id in idx:
                bpdata=data_list[id]
                new_data=data_list[:,id]
                new_data=reshape(new_data, [nx, ny])
                bpdata.update_data(new_data)
                bpdata.write_to_bpch2_file(ifunit)
            bpch2_w.close_bpch2_w(ifunit)

        # write data back 
        
     
if (__name__=="__main__"):
    import co2_emission_std as co2em
    import time_module as tm
    new_lon, new_lat=co2em.GET_GRID()
    rsf=geos_chem_restart_file('restart.jan2003.kalman.borealasia')
    real_ntracers=28
    yyyy=2006
    mm=1
    dd=1
    doy=tm.day_of_year(yyyy=yyyy, mm=mm, dd=dd)
    tau0=tm.doy_to_tai85(doy, sec=0, yyyy=yyyy)
    tau0=tau0/3600.0
    full_restart_name='restart.20060101'
    rsf.mod_restart_file(1, full_restart_name,real_ntracers, tau0, \
                         fixed_value=1.0e-8, do_regrid=True, new_lon=new_lon, new_lat=new_lat)
