from numpy import *
import time_module as tm
import geos_chem_def as gcdf
def create_new_input_file(run_step,\
                          YYYY,\
                          DOY,\
                          pbuse,\
                          pb_start=1,\
                          member_start=1, \
                          member_end=12, \
                          co2flnm='CO2_EMISSION_EN',\
                          tmpfile="input.geos.temp", \
                          newfile="input.geos.new" , \
                          em_step=None, \
                          em_doy=None,\
                          em_yyyy=None,\
                          em_pbuse=None,\
                          do_bk_run=4,\
                          data_path=gcdf.data_path,\
                          run_path=gcdf.run_path,\
                          **keywords \
                          ):
    
    """ create the new input to drive the ensemble run
    YYYY: years of the time serier
    DOY:  doys of the time serier
    em_doy: the date need emission
    """
    
    ntime=size(YYYY)
    if ('time_start' in keywords): 
        ts=keywords['time_start']
        syyyy=ts[0:4]
        yst=int(syyyy)
        dst=ts[4:7]
        dst=int(dst)
        
        # mmst=int(smm)
        # sdd=ts[6:8]
        # ddst=int(sdd)
        # dst=tm.day_of_year(yst, mmst, ddst)
    else:
        yst=YYYY[0]
        dst=DOY[0]

    tst=tm.doy_to_utc(dst, sec=0, yyyy=yst)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '')
    print tst
    if (em_step==None):
        iday0=DOY[0]
        for iday in DOY[1:]:
            if (iday>iday0):
                em_step=iday-iday0
                break
            else:
                iday0=iday

    if ('time_end' in keywords): 
        te=keywords['time_end']
        syyyy=te[0:4]
        yed=int(syyyy)
        smm=te[4:6]
        mmst=int(smm)
        sdd=te[6:8]
        ddst=int(sdd)
        ded=tm.day_of_year(yst, mmst, ddst)
    else:
        yed=YYYY[1] #-1
        ded=DOY[1]  #-1
    
    # the endtime should be + timestep
    # yed, ded=tm.next_doy(yed, ded, em_step)
    
    tend=tm.doy_to_utc(ded, sec=0, yyyy=yed)
    tend=tend.replace('-', '')
    tend=tend.replace(':', '')
    print tend

    fin=open(tmpfile, "r")
    lines=fin.readlines()
    print  len(lines)
    fin.close()
    fout=open(newfile, "w")
    section_start=0
    colwidth=25
    line_count=0
    ntracer=member_end-member_start+1
    # enafix=r'EN%4.4d-EN%4.4d' % (member_start, member_end)
    enafix=r'ST%3.3d.EN%4.4d-EN%4.4d' % (run_step, member_start, member_end)
    
    for line in lines:
        if ("SIMULATION MENU" in line):
            section_start=1
            line_count=0
            fout.write(line)
        elif ("TRACER MENU" in line):
            section_start=2
            line_count=0
            fout.write(line)
        elif ("ND51 MENU" in line):
            section_start=3
            line_count=0
            fout.write(line)
        elif("OUTPUT MENU" in line):
            section_start=10
            line_count=0
            fout.write(line)
        elif ("ENSEMBLE MENU" in line):
            section_start=4
            line_count=0
            fout.write(line)
        elif ("DIAGNOSTIC MENU" in line):
             section_start=5
             line_count=0
             fout.write(line)
        elif("GAMAP MENU" in line):
            section_start=6
            line_count=0
            fout.write(line)
        elif("PROD & LOSS MENU" in line):
            section_start=7
            line_count=0
            fout.write(line)
        elif ("ND50 MENU" in line):
            section_start=8
            line_count=0
            fout.write(line)
            
        elif (section_start>0):
            line_count=line_count+1
            line_head=line[0:colwidth]
            line_left=line[colwidth:]
            if (section_start==1):
                if (line_count==1):
                    new_line=line_head+" "+tst+"   "
                    fout.write(new_line+"\n")
                elif (line_count==2):
                    line_head=line[0:colwidth]
                    new_line=line_head+" "+tend+"   "
                    fout.write(new_line+"\n")
                elif (line_count==3):
                    print line_left
                    line_left=line_left.replace('$RUNPATH', run_path)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==4):
                    print line_left
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==6):  
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==13):
                    line_left=line_left.replace('$RUNPATH', run_path)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line[0:2]=='--'):
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            elif (section_start==2): # tracer menu 
                if (line_count==1):
                    fout.write(line)
                elif(line_count==2):
                    line_head=line[0:colwidth]
                    snum=r'%4d' % (ntracer)
                    new_line=line_head+" "+snum+"   "+"\n"
                    fout.write(new_line)
                elif (line_count==3):
                    fout.write(line)
                elif (line_count==4):
                    line_head=line[0:colwidth]
                    tracer_no=line[colwidth:colwidth+5]
                    line_left=line[colwidth+5:]
                    tracer_id=1
                    tracers=list()
                    
                    for itracer in range(member_end-member_start+1):
                        st1=r'Tracer #%3.3d ENSEM' % (itracer+member_start)
                        lst1=len(st1)
                        new_head=st1[:]+' '*(colwidth-1-lst1)+":"
                        st1=r'%5d' % (tracer_id)
                        lst1=len(st1)
                        new_tracer_no=(5-lst1)*' '+st1
                        new_line=new_head+new_tracer_no+line_left
                        fout.write(new_line)
                        tracers.append(tracer_id)
                        tracer_id=tracer_id+1
                        
                elif (line[0:2]=='--'):
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    pass
            elif (section_start==3): # ND51 menu 
                if (line_count==1):
                    fout.write(line)
                elif (line_count==2):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==4):
                    line_left=""

                    for itracer in [1, member_end-member_start+1]:
                        stracer=r' %d' % (itracer)
                        line_left=line_left+stracer
                                            
                    line_left=line_left+' '+'199' # 196 198 199 200 201 output other information
                    new_line=line_head+line_left+' \n'
                    print new_line
                    # tx=raw_input()
                    
                    fout.write(new_line)
                    
                    
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            elif (section_start==4): #ensemble  menu
                if (line_count==1):
                    st1=r' %d  ' % member_start
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                elif (line_count==2):
                    st1=r' %d  ' % member_end
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                elif (line_count==3):
                    st1=r' %d  ' % do_bk_run
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)

                elif (line_count==4):
                    st1=r' %d ' % ntime
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                elif (line_count==5):
                    st1=""
                    if (em_doy==None):
                        em_doy=DOY
                    for doy in em_doy:
                        sdoy=r'%d' % doy
                        st1=st1+' '+sdoy
                    st1=st1+' '
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                elif (line_count==6):
                    st1=""
                    if (em_yyyy==None):

                        for yyyy in YYYY:
                            syyyy=r'%d' % yyyy
                            st1=st1+' '+syyyy
                        
                    else:

                        for yyyy in em_yyyy:
                            syyyy=r'%d' % yyyy
                            st1=st1+' '+syyyy
                    
                    st1=st1+' '
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                
                elif (line_count==7):
                    st1=r' %d  ' % pb_start
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                
                elif (line_count==8):
                    st1=""
                    if (em_pbuse==None):
                        for pbuse in pbuse:
                            spbuse=r'%d' % pbuse
                            st1=st1+' '+spbuse
                        
                    else:
                        
                        for pbuse in em_pbuse:
                            spbuse=r'%d' % pbuse
                            st1=st1+' '+spbuse
                    
                    st1=st1+' '
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                
                
                elif (line_count==9): # 'flnm'
                    new_line=line_head+' '+co2flnm.strip()+'\n'
                    fout.write(new_line)
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
            elif (section_start==47): #P/L  menu
                line_head=line[0:colwidth]
                tracer_no=line[colwidth:colwidth+5]
                line_left=line[colwidth+5:]
                tracer_id=1
                tracers=list()
                
                if (line_count<=3):
                    #  print line_count, line
                    fout.write(line)
                elif (line_count==4):
                    st1=r'%d' % (member_end-member_start+1+5)
                    new_line=line_head+' '+st1+' \n'
                    fout.write(new_line)
                elif (line_count==5):
                    fout.write(line)
                    for itracer in range(1, member_end-member_start+1+5):
                        st1=r'Prod/Loss Family #%d' % (itracer+member_start)
                        lst1=len(st1)
                        new_head=st1[:]+' '*(colwidth-1-lst1)+":"
                        st1=r'%5d' % (tracer_id)
                        lst1=len(st1)
                        new_tracer_no=(5-lst1)*' '+st1
                        line_left='PCO: CO%3.3d' % (itracer)
                        new_line=new_head+' '+line_left+'\n'
                        fout.write(new_line)
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    pass
                    
            elif(section_start==5): #diagnostic  menu
                if (line_count==1):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            elif(section_start==6): #gamap  menu
                if (line_count==1):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==2):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            
            elif (section_start==8): # ND50 menu 
                if (line_count==1):
                    fout.write(line)
                elif (line_count==2):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==4):
                    line_left=""
                    for itracer in [1, member_end-member_start+1]:
                        stracer=r' %d' % (itracer)
                        line_left=line_left+stracer
                                            
                    line_left=line_left+' '+'399'  # output other information
                    new_line=line_head+line_left+' \n'
                    fout.write(new_line)
                    
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            
            elif (section_start==10): #output menu
                
                if (line_count==1):
                    id_cnt=0
                    yyyy_list=list()
                    mm_list=list()
                    dd_list=list()
                    
                    for idoy in DOY:
                        cyyyy, cmm, cdd=tm.doy_to_time_array(idoy, YYYY[id_cnt])
                        yyyy_list.append(cyyyy)
                        mm_list.append(cmm)
                        dd_list.append(cdd)
                        id_cnt=id_cnt+1
                    
                    yyyy_list=array(yyyy_list)
                    mm_list=array(mm_list)
                    dd_list=array(dd_list)
                    
                       
                if ((line_count>=1) and (line_count<13)):
                    line_left=line_left.replace('3', '0')
                    sel_mm=where(mm_list==line_count)
                    # print size(sel_mm), sel_mm
                    if (size(sel_mm)>0):
                        if (size(sel_mm)==1):
                            sel_mm=[squeeze(sel_mm)]
                        else:
                            sel_mm=squeeze(sel_mm)
                            
                        for isel in sel_mm:
                                            
                            sel_dd=dd_list[isel]
                            ispace=0
                            # skip the line
                            
                            for ichar in line_left:
                                if (ichar<>' '):
                                    break
                                else:
                                    ispace=ispace+1
                            
                            sel_dd=sel_dd+ispace
                            
                            print line_left
                            
                            if (sel_dd>=len(line_left)):
                                sel_dd=len(line_left)-1
                        
                            if (sel_dd==1):
                                tline_1='3'
                            
                            else:
                                tline_1=line_left[0:sel_dd-1]
                                tline_1=tline_1+'3'
                            
                                tline_1=tline_1+line_left[sel_dd:]
                            line_left=tline_1
                            
                    new_line=line_head+line_left
                    fout.write(new_line)        
                elif (line[0:2]=='--'):  # the section end line          
                    line_count=0
                    section_start=0
                    fout.write(line)

            else:
                fout.write(line)
                    
            
            

        else:
            fout.write(line)
    
    fout.close()
    
    print 'reach the end'
    #  tx2=raw_input()
    
if (__name__=="__main__"):
    create_new_input_file([2003, 2003, 2003], [1, 9, 17])


                  
                
