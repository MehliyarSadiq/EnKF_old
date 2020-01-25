import geos_chem_def as dfs
class menu_item():
    pass

class input_menu:
    def __init__(self):
        self.dict={}
    def getvals(self, names, short_names=None):
        vals={}
        for it in range(len(names)):
            valname=names[it]
            if (short_names==None):
                short_name=valname.strip()
            else:
                short_valname=short_names[it]

            try:
                val=self.dict[valname.strip()]
            except KeyError:
                val=""
                
            if (val=='T'):
                val=True
            elif (val=='F'):
                val=False
            else:
                if (val.isdigit()):
                    val=int(val)
                else:
                    try:
                        tmp_val=float(val)
                        val=tmp_val
                    except ValueError:
                        val=val.strip()
                
            vals.update({short_valname:val})
        return vals

        
    
def read_geos_input(flnm=dfs.def_input_file):
    FIRSTCOL = 26
    menus_names=['SIMULATION MENU', 'EMISSIONS MENU']
    menu_list={}
    
    fin=open(flnm.strip(), 'r')
    head_name=""
    lines=fin.readlines()
    fin.close()
    read_val=False
    for line in lines[1:]:
        if (line[0]=='#' or line[0]=='-'):
            pass
        elif (line[0:3]=='%%%'):
            # read in section name
            if (read_val):
                menu_list.update({head_name:read_menu})
            line=line.strip()
            
            terms=line.split('%%%')
            head_name=terms[1].strip()
            
            if (head_name in menus_names):
                read_val=True
                read_menu=input_menu()
            else:
                read_val=False
        elif (line.strip=='END OF FILE'):
            break
        else:
            if (read_val):
                
                title,val=line.split(':')
                title=title.strip()
                val=val.strip()
                read_menu.dict.update({title:val})
            else:
                pass
    return menu_list
def emission_menu_parse(em):
    names=['Turn on emissions?', 'Emiss timestep (min)', 'Include anthro emiss?', '=> Scale 1985 to year',\
           '=> Use EMEP emissions?', 'Include biofuel emiss?', 'Include biogenic emiss?', '=> Use MEGAN inventory?',\
           'Include biomass emiss?', '=> Seasonal biomass?', '=> Scaled to TOMSAI?', '=> Use aircraft NOx?', \
           '=> Use lightning NOx', '=> Use soil NOx', 'Use ship SO2 emissions?', 'Use EPA/NEI99 emissions?:', 'Use AVHRR-derived LAI?']

    short_names=['do_em', 'timestep', 'do_anthro_em', '=> anthro_em_scale_year',\
           'do_anthor_em_emep', 'do_biofuel_em?', 'do_biogenic_em', 'biogenic_em_megan',\
           'do_biomass_em', 'do_biomass_em_sea', 'do_biomass_em_scale', 'do_nox_aircraft', 'do_nox_lightning', 'do_nox_soil', 'do_ship_so2_em', 'do_epa_nei99_em', 'do_avhrr_lai']
    
    em_details=em.getvals(names, short_names)
    return em_details
    
   
def simulation_menu_parse(sm):
    names=['Start YYYYMMDD, HHMMSS', 'End   YYYYMMDD, HHMMSS', \
           'Run directory', 'Input restart file',\
           'Make new restart file?', 'Output restart file(s)',  \
           'Root data directory', '=> GCAP       subdir',       \
           '=> GEOS-1     subdir', '=> GEOS-STRAT subdir',     \
           '=> GEOS-3     subdir', ' => GEOS-4     subdir',    \
           '=> GEOS-5     subdir', 'Dir w/ 1x1 emissions etc',  \
           'Temporary directory', 'Unzip met fields?',\
           'Wait for met fields?', 'Global offsets I0, J0']

    short_names=['start', 'end', \
                 'run_dir', 'input_restart', \
                 'do_new_restart',    \
                 'out_restart_files', \
                 'data_dir', 'gcap_dir', \
                 'geos1_dir', 'geos_strat',\
                 'geos3_dir','geos4_dir', \
                 'geos5_dir', 'all_1_1_dir',\
                 'tmp_dir', 'do_unzip_met', \
                 'do_wait_met', 'glb_offset']
    
                 
    sm_details=sm.getvals(names, short_names)
    return sm_details
    
    
if (__name__=='__main__'):
    menus=read_geos_input()
    print(menus.keys())
    em=menus['EMISSIONS MENU']
    sm=menus['SIMULATION MENU']
    em_details=emission_menu_parse(em)
    print(em_details)
    sm_details=simulation_menu_parse(sm)
    print(sm_details)
    
    


