import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from scipy import signal
import os
from scipy.stats import linregress
mpl.use('tkagg')


##########################################-- GET ALL FILES EXCEPT FILES' VS30 IS "NONE" --##########################################
####################################################################################################################################

def getFiles(data_folder): 

    os.chdir(data_folder)
    folder_list=[]
    file_list=[]
    for i in os.listdir():
        if i.startswith("._")==False:
            folder_list.append(i)

    for j in folder_list:
        os.chdir(data_folder + "/" + j)
        for k in os.listdir():
            if k.startswith("._")==False:
                file_list.append(j+ "/" +k)
    
    for l in range(0, len(file_list)):
        with open(data_folder + "/" + file_list[l], "r", encoding="utf-8") as file:
            test_value="VS30_M/S: "

            file.seek(file.read().find(test_value)+len(test_value))
            Vs30_value=file.readline()
            Vs30_value=Vs30_value[:-1]
            Vs30_value=[''.join(Vs30_value)]
            Vs30_value=Vs30_value[0]
            if Vs30_value=="None":
                os.remove(data_folder + "/" + file_list[l])
                file_list.remove(data_folder + "/" + file_list[l])

    return file_list


##################################################-- FIND TIME INTERVAL OF DATA --##################################################
####################################################################################################################################

def findInterval(file_use_quotationmark_with_folder):
    interval=[]
    with open(file_use_quotationmark_with_folder, "r", encoding="utf-8") as file:

        test_value="INTERVAL_S: "

        file.seek(file.read().find(test_value)+len(test_value))
        interval=file.readline()
        interval=interval[:-1]
        interval=[''.join(interval)]
        dt=interval[0]
        dt=float(dt)

    return dt

    
############################################-- RAW ACCELERATION RECORDS OF EARTHQUAKE --############################################
#################################################################################################################################### 
  
def rawAccRec(file_use_quotationmark_with_folder):
    ag_txt=np.loadtxt(file_use_quotationmark_with_folder, skiprows=64)
    groundacc=ag_txt[:]
    ags=groundacc.flatten("C")
    
    return ags


#########################################################-- TIME SERIES --##########################################################
####################################################################################################################################

def timeSer(file_use_quotationmark_with_folder):
    ags=rawAccRec(file_use_quotationmark_with_folder)
    t_amount=len(ags)
    dt=findInterval(file_use_quotationmark_with_folder)
    initial_time=dt
    end_time=initial_time+dt*t_amount
    t=np.arange(initial_time, end_time, dt)

    return t


#############################-- RAW VELOCITY RECORDS TRANSFORMED WITH CONSTANT ACCELERATION METHOD --###############################
#################################################################################################################################### 

def rawVelRec(file_use_quotationmark_with_folder):
    vgs=list()
    ag_txt=np.loadtxt(file_use_quotationmark_with_folder, skiprows=64)
    groundacc=ag_txt[:]
    ags=groundacc.flatten("C")
    t_amount=len(ags)
    dt=findInterval(file_use_quotationmark_with_folder)
    vgs.insert(0,0)
    
    for i in range(1,t_amount):
        vgs.append(vgs[i-1]+(((ags[i-1]+ags[i])/2)*dt))
    
    return vgs


#############################-- RAW DISPLACEMENT RECORDS TRANSFORMED WITH CONSTANT VELOCITY METHOD --###############################
####################################################################################################################################

def rawDispRec(file_use_quotationmark_with_folder):
    vgs=rawVelRec(file_use_quotationmark_with_folder)
    ags=rawAccRec(file_use_quotationmark_with_folder)
    t_amount=len(ags)
    dt=findInterval(file_use_quotationmark_with_folder)
    dgs=list()
    dgs.insert(0,0)

    for i in range(1,t_amount):
        dgs.append(dgs[i-1]+(((vgs[i-1]+vgs[i])/2)*dt))
    
    return dgs


#################################-- BASELINE CORRECTION FOR RAW ACCELERATION DATA OF EARTHQUAKE --##################################
####################################################################################################################################

def BaselineCorrection(file_use_quotationmark_with_folder):
    ags=rawAccRec(file_use_quotationmark_with_folder)
    newAgs=signal.detrend(ags)
    t_amount=len(newAgs)
    dt=findInterval(file_use_quotationmark_with_folder)
    newVgs=list()
    newVgs.insert(0,0)
    for i in range(1,t_amount):
        newVgs.append(newVgs[i-1]+(((newAgs[i-1]+newAgs[i])/2)*dt))

    newDgs=list()
    newDgs.insert(0,0)
    for j in range(1,t_amount):
        newDgs.append(newDgs[j-1]+(((newVgs[j-1]+newVgs[j])/2)*dt))
    
    return newAgs, newVgs, newDgs


########################################################-- ARIAS INTENSTY --########################################################
####################################################################################################################################

def AriasIntensity(file_use_quotationmark_with_folder):
    
    g=981
    newAgs, newVgs, newDgs=BaselineCorrection(file_use_quotationmark_with_folder)
    dt=findInterval(file_use_quotationmark_with_folder)
    cumSum=[np.cumsum(np.abs(newAgs)*np.abs(newAgs))]
    AI=[]
    for i in cumSum:
        AI.append((np.pi/(2*g))*i*dt)
    AI=AI[0]

    return AI


#################################################-- CUMULATIVE ABSOLUTE VELOCITY --#################################################
####################################################################################################################################

def CumAbsVel(file_use_quotationmark_with_folder):

    newAgs, newVgs, newDgs=BaselineCorrection(file_use_quotationmark_with_folder)
    dt=findInterval(file_use_quotationmark_with_folder)
    cumSum=[np.cumsum(np.abs(newAgs*dt))]
    CAV=[]
    for i in cumSum:
        CAV.append(i)
    CAV=CAV[0]
    
    return CAV


#####################################################-- SIGNIFICANT DURATION --#####################################################
####################################################################################################################################

def SigDur(file_use_quotationmark_with_folder):
    t=timeSer(file_use_quotationmark_with_folder)
    AI=AriasIntensity(file_use_quotationmark_with_folder)
    
    AI_05=AI[-1]*0.05
    idx_05=np.argmin(np.abs(AI-AI_05))

    AI_95=AI[-1]*0.95
    idx_95=np.argmin(np.abs(AI-AI_95))

    t_05=t[idx_05]
    t_95=t[idx_95]

    SD_05_95=t_95-t_05

    return t_05, AI_05, t_95, AI_95, SD_05_95


########################################################-- VS30 --##################################################################
####################################################################################################################################

def Vs30(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="VS30_M/S: "

        file.seek(file.read().find(test_value)+len(test_value))
        vS30=file.readline()
        vS30=vS30[:-1]
        vS30=[''.join(vS30)]
        vS30=vS30[0]
        vS30=float(vS30)

    return vS30


#####################################################-- EPICENTRAL DISTANCE --######################################################
####################################################################################################################################

def Repi(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="EPICENTRAL_DISTANCE_KM: "

        file.seek(file.read().find(test_value)+len(test_value))
        rEpi=file.readline()
        rEpi=rEpi[:-1]
        rEpi=[''.join(rEpi)]
        rEpi=rEpi[0]
        rEpi=float(rEpi)
    
    return rEpi


##########################################################-- STREAM CODE --#########################################################
####################################################################################################################################

def streamCode(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="STREAM: "

        file.seek(file.read().find(test_value)+len(test_value))
        stream_code=file.readline()
        stream_code=stream_code[:-1]
        stream_code=[''.join(stream_code)]
        stream_code=stream_code[0]
    
    return stream_code


####################################################-- STATION LATITUDE DEGREE --###################################################
####################################################################################################################################

def stationLatitude(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="STATION_LATITUDE_DEGREE: "

        file.seek(file.read().find(test_value)+len(test_value))
        s_Latitude=file.readline()
        s_Latitude=s_Latitude[:-1]
        s_Latitude=[''.join(s_Latitude)]
        s_Latitude=s_Latitude[0]
    
    return s_Latitude


####################################################-- STATION LONGITUDE DEGREE --##################################################
####################################################################################################################################

def stationLongitude(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="STATION_LONGITUDE_DEGREE: "

        file.seek(file.read().find(test_value)+len(test_value))
        s_Longitude=file.readline()
        s_Longitude=s_Longitude[:-1]
        s_Longitude=[''.join(s_Longitude)]
        s_Longitude=s_Longitude[0]
    
    return s_Longitude

##########################################################-- EVENT DEPTH --#########################################################
####################################################################################################################################

def eventDepth(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="EVENT_DEPTH_KM: "

        file.seek(file.read().find(test_value)+len(test_value))
        e_depth=file.readline()
        e_depth=e_depth[:-1]
        e_depth=[''.join(e_depth)]
        e_depth=e_depth[0]
    
    return e_depth

######################################################-- SITE CLASSIFICATION --#####################################################
####################################################################################################################################

def siteClass(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="SITE_CLASSIFICATION_EC8: "

        file.seek(file.read().find(test_value)+len(test_value))
        s_class=file.readline()
        s_class=s_class[:-1]
        s_class=[''.join(s_class)]
        s_class=s_class[0]
    
    return s_class

########################################################-- EVENT DATE TIME --#######################################################
####################################################################################################################################

def dateTime(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS: "

        file.seek(file.read().find(test_value)+len(test_value))
        d_time=file.readline()
        d_time=d_time[:-1]
        d_time=[''.join(d_time)]
        d_time=d_time[0]
    
    return d_time

###################################################-- PEAK GROUND ACCELERATION --###################################################
####################################################################################################################################

def PGA(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="PGA_CM/S^2: "

        file.seek(file.read().find(test_value)+len(test_value))
        pga=file.readline()
        pga=pga[:-1]
        pga=[''.join(pga)]
        pga=pga[0]
    
    return pga


#######################################################-- MOMENT MAGNITUDE --#######################################################
####################################################################################################################################

def momentMagnitude(path):
    with open(path, "r", encoding="utf-8") as file:
        test_value="MAGNITUDE_W: "

        file.seek(file.read().find(test_value)+len(test_value))
        mW=file.readline()
        mW=mW[:-1]
        mW=[''.join(mW)]
        mW=mW[0]
        mW=float(mW)

    return mW


#####################################################-- PLOT EARTHQUAKE DATA --#####################################################
####################################################################################################################################

def plotData(path):
    files=getFiles(path)
    for i in range(0, len(files),3):
        pgaE=PGA(path+ "/" + files[i])
        pgaN=PGA(path+ "/" + files[i+1])
        pgaZ=PGA(path+ "/" + files[i+2])
        s_latitude=stationLatitude(path+ "/" + files[i])
        s_longitude=stationLongitude(path+ "/" + files[i])
        e_depth=eventDepth(path+ "/" + files[i])
        s_class=siteClass(path+ "/" + files[i])
        d_time=dateTime(path+ "/" + files[i])
        mw=momentMagnitude(path+ "/" + files[i])
        rEpi=Repi(path+ "/" + files[i])
        vS30=Vs30(path+ "/" + files[i])
        ags=rawAccRec(path+ "/" + files[i])
        ags1=rawAccRec(path+ "/" + files[i+1])
        ags2=rawAccRec(path+ "/" + files[i+2])
        vgs=rawVelRec(path+ "/" + files[i])
        vgs1=rawVelRec(path+ "/" + files[i+1])
        vgs2=rawVelRec(path+ "/" + files[i+2])
        dgs=rawDispRec(path+ "/" + files[i])
        dgs1=rawDispRec(path+ "/" + files[i+1])
        dgs2=rawDispRec(path+ "/" + files[i+2])
        t=timeSer(path+ "/" + files[i])
        AI=AriasIntensity(path+ "/" + files[i])
        AI1=AriasIntensity(path+ "/" + files[i+1])
        AI2=AriasIntensity(path+ "/" + files[i+2])
        t_05, AI_05, t_95, AI_95, SD_05_95=SigDur(path+ "/" + files[i])
        t1_05, AI1_05, t1_95, AI1_95, SD1_05_95=SigDur(path+ "/" + files[i+1])
        t2_05, AI2_05, t2_95, AI2_95, SD2_05_95=SigDur(path+ "/" + files[i+2])
        newAgs, newVgs, newDgs=BaselineCorrection(path+ "/" + files[i])
        newAgs1, newVgs1, newDgs1=BaselineCorrection(path+ "/" + files[i+1])
        newAgs2, newVgs2, newDgs2=BaselineCorrection(path+ "/" + files[i+2])
        CAV=CumAbsVel(path+ "/" + files[i])
        CAV1=CumAbsVel(path+ "/" + files[i+1])
        CAV2=CumAbsVel(path+ "/" + files[i+2])
        fig, axs = plt.subplots(4,3,figsize=(20,10))
        fig.subplots_adjust(hspace=0.35, wspace=0.23, left=0.08, bottom=0.07, right=0.96, top=0.93)

        with open(path+ "/" + files[i], "r", encoding="utf-8") as file:
            test_value="STATION_CODE: "
            

            file.seek(file.read().find(test_value)+len(test_value))
            station_code=file.readline()
            station_code=station_code[:-1]
            station_code=[''.join(station_code)]
            station_code=station_code[0]
        
        fig.suptitle("Aegean Earthquake" + "  -" + "Station Code " + station_code + "-" , fontsize=20, color="red")

        axs[0,0].plot(t,ags, "b", linewidth=2, label="Raw Data")
        axs[0,0].plot(t, newAgs, "r", linewidth=0.5, label="Baseline Correction") 
        axs[0,0].set_xlabel("Time (sec)")
        axs[0,0].set_ylabel("Acceleration (cm/sec^2)")
        stream_code=streamCode(path+ "/" + files[i])
        axs[0,0].legend(title=stream_code, fontsize="x-small")
        axs[0,0].grid(True)

        axs[1,0].plot(t, ags1, 'b', linewidth=2, label="Raw Data")
        axs[1,0].plot(t, newAgs1, 'r', linewidth=0.5, label="Baseline Correction")
        axs[1,0].set_xlabel("Time (sec)")
        axs[1,0].set_ylabel("Acceleration (cm/sec^2)")
        stream_code=streamCode(path+ "/" + files[i+1])
        axs[1,0].legend(title=stream_code, fontsize="x-small")
        axs[1,0].grid(True)
        
        axs[2,0].plot(t, ags2, 'b', linewidth=2, label="Raw Data")
        axs[2,0].plot(t, newAgs2, 'r', linewidth=0.5, label="Baseline Correction")
        axs[2,0].set_xlabel("Time (sec)")
        axs[2,0].set_ylabel("Acceleration (cm/sec^2)")
        stream_code=streamCode(path+ "/" + files[i+2])
        axs[2,0].legend(title=stream_code, fontsize="x-small")
        axs[2,0].grid(True)
        
        axs[0,1].plot(t, vgs, 'b', linewidth=2, label="Raw Data")
        axs[0,1].plot(t, newVgs, 'r', linewidth=1, label="Baseline Correction")
        axs[0,1].set_xlabel("Time (sec)")
        axs[0,1].set_ylabel("Velocity (cm/sec)")
        stream_code=streamCode(path+ "/" + files[i])
        axs[0,1].legend(title=stream_code, fontsize="x-small")
        axs[0,1].grid(True)
        
        axs[1,1].plot(t, vgs1, 'b', linewidth=2, label="Raw Data")
        axs[1,1].plot(t, newVgs1, 'r', linewidth=1, label="Baseline Correction")
        axs[1,1].set_xlabel("Time (sec)")
        axs[1,1].set_ylabel("Velocity (cm/sec)")
        stream_code=streamCode(path+ "/" + files[i+1])
        axs[1,1].legend(title=stream_code, fontsize="x-small")
        axs[1,1].grid(True)
        
        axs[2,1].plot(t, vgs2, 'b', linewidth=2, label="Raw Data")
        axs[2,1].plot(t, newVgs2, 'r', linewidth=1, label="Baseline Correction")
        axs[2,1].set_xlabel("Time (sec)")
        axs[2,1].set_ylabel("Velocity (cm/sec)")
        stream_code=streamCode(path+ "/" + files[i+2])
        axs[2,1].legend(title=stream_code, fontsize="x-small")
        axs[2,1].grid(True)
        
        axs[0,2].plot(t, dgs, 'b', linewidth=2, label="Raw Data")
        axs[0,2].plot(t, newDgs, 'r', linewidth=1, label="Baseline Correction")
        axs[0,2].set_xlabel("Time (sec)")
        axs[0,2].set_ylabel("Displacement (cm)")
        stream_code=streamCode(path+ "/" + files[i])
        axs[0,2].legend(title=stream_code, fontsize="x-small")
        axs[0,2].grid(True)
        
        axs[1,2].plot(t, dgs1, 'b', linewidth=2, label="Raw Data")
        axs[1,2].plot(t, newDgs1, 'r', linewidth=1, label="Baseline Correction")
        axs[1,2].set_xlabel("Time (sec)")
        axs[1,2].set_ylabel("Displacement (cm)")
        stream_code=streamCode(path+ "/" + files[i+1])
        axs[1,2].legend(title=stream_code, fontsize="x-small")
        axs[1,2].grid(True)
        
        axs[2,2].plot(t, dgs2, 'b', linewidth=2, label="Raw Data")
        axs[2,2].plot(t, newDgs2, 'r', linewidth=1, label="Baseline Correction")
        axs[2,2].set_xlabel("Time (sec)")
        axs[2,2].set_ylabel("Displacement (cm)")
        stream_code=streamCode(path+ "/" + files[i+2])
        axs[2,2].legend(title=stream_code, fontsize="x-small")
        axs[2,2].grid(True)

        axs[3,0].plot(t,AI,'b', label=streamCode(path+ "/" + files[i]))
        axs[3,0].plot( t_05, AI_05, 'b', marker='o',label="AI05 and AI95")
        axs[3,0].plot(t_95, AI_95, 'b', marker='o')
        axs[3,0].plot(t,AI1,'r', label=streamCode(path+ "/" + files[i+1]))
        axs[3,0].plot( t1_05, AI1_05, 'r', marker='o', label="AI05 and AI95")
        axs[3,0].plot(t1_95, AI1_95, 'r', marker='o')
        axs[3,0].plot(t,AI2,'m', label=streamCode(path+ "/" + files[i+2]))
        axs[3,0].plot( t2_05, AI2_05, 'm', marker='o', label="AI05 and AI95")
        axs[3,0].plot(t2_95, AI2_95, 'm', marker='o')
        axs[3,0].set_xlabel("Time (sec)")
        axs[3,0].set_ylabel("Arias Intensity (cm/sec)")
        axs[3,0].legend(fontsize="small")
        axs[3,0].grid(True)

        axs[3,1].plot(t,CAV,'b', label=streamCode(path+ "/" + files[i]))
        axs[3,1].plot(t,CAV1,'r', label=streamCode(path+ "/" + files[i+1]))
        axs[3,1].plot(t,CAV2,'m', label=streamCode(path+ "/" + files[i+2]))
        axs[3,1].set_xlabel("Time (sec)")
        axs[3,1].set_ylabel("Cumulative Absolute Velocity (cm/sec)")
        axs[3,1].legend(fontsize="small")
        axs[3,1].grid(True)

        axs[3,2].text(1, 11, 'Station Latitude: ' + str(s_latitude))
        axs[3,2].text(1, 10, 'Station Longitude: ' + str(s_longitude))
        axs[3,2].text(1, 9, 'Event Date Time: ' + str(d_time))
        axs[3,2].text(1, 8, 'Moment Magnitude (Mw): ' + str(mw))
        axs[3,2].text(1, 7, 'Event Depth: ' + str(e_depth) + ' km')
        axs[3,2].text(1, 6, 'Vs30: ' + str(vS30) + ' m/s')
        axs[3,2].text(1, 5, 'Site Classification: ' + str(s_class))
        axs[3,2].text(1, 4, 'PGA_EW: ' + str(pgaE) + r' $cm/s^2$')
        axs[3,2].text(1, 3, 'PGA_NS: ' + str(pgaN) + r' $cm/s^2$')
        axs[3,2].text(1, 2, 'PGA_U: ' + str(pgaZ) + r' $cm/s^2$')
        axs[3,2].text(1, 1, 'Epicentral Distance: ' + str(rEpi) + ' km')
        axs[3,2].plot()
        axs[3,2].axis([0, 12, 0, 12])

        # plt.show()
        
        # plt.savefig("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/Plots" + "/" + "Aegean Earthquake" + "  -" + "Station Code " + station_code + "-"  + ".png")


#########################################-- PLOT SIGNIFICANT DURATION - DISTANCE AND VS30 --########################################
####################################################################################################################################

def plotSigDur(path):
    rEpiList=list()
    vS30List=list()
    SD_05_95List=list()
    SD_05_95List1=list()
    SD_05_95List2=list()
    SD_05_95tList=list()

    files=getFiles(path)
    for i in range(0, len(files),3):
        rEpi=Repi(path+ "/" + files[i])
        rEpiList.append(rEpi)
        vS30=Vs30(path+ "/" + files[i])
        vS30List.append(vS30)
        t_05, AI_05, t_95, AI_95, SD_05_95=SigDur(path+ "/" + files[i])
        SD_05_95List.append(SD_05_95)
        t1_05, AI1_05, t1_95, AI1_95, SD1_05_95=SigDur(path+ "/" + files[i+1])
        SD_05_95List1.append(SD1_05_95)
        t2_05, AI2_05, t2_95, AI2_95, SD2_05_95=SigDur(path+ "/" + files[i+2])
        SD_05_95List2.append(SD2_05_95)
        SD_05_95tList.append((SD_05_95+SD1_05_95)/2)

    fig, axs = plt.subplots(3,3,figsize=(20,10))
    fig.subplots_adjust(hspace=0.23, wspace=0.19, left=0.04, right=0.98, bottom=0.06, top=0.95)

    from scipy.stats import linregress
    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList,SD_05_95List)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    axs[0,0].scatter(rEpiList, SD_05_95List, label="SD5-95")
    axs[0,0].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    axs[0,0].set_ylabel("Significant Duration 5-95 (sec)")
    axs[0,0].set_xlabel("Epicenteral Distance (km)")
    stream_code=streamCode(path+ "/" + files[0])
    axs[0,0].legend(title=stream_code, fontsize="x-small")
    axs[0,0].grid(True)
    axs[0,0].set_title("R - Ds5-95", fontsize=20, color="red")

    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList,SD_05_95List1)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    axs[1,0].scatter(rEpiList, SD_05_95List1, label="SD5-95")
    axs[1,0].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    axs[1,0].set_ylabel("Significant Duration 5-95 (sec)")
    axs[1,0].set_xlabel("Epicenteral Distance (km)")
    stream_code=streamCode(path+ "/" + files[1])
    axs[1,0].legend(title=stream_code, fontsize="x-small")
    axs[1,0].grid(True)

    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList,SD_05_95List2)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    axs[2,0].scatter(rEpiList, SD_05_95List2, label="SD5-95")
    axs[2,0].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    axs[2,0].set_ylabel("Significant Duration 5-95 (sec)")
    axs[2,0].set_xlabel("Epicenteral Distance (km)")
    stream_code=streamCode(path+ "/" + files[2])
    axs[2,0].legend(title=stream_code, fontsize="x-small")
    axs[2,0].grid(True)

    gradient, intercept, r_value, p_value, std_err=linregress(vS30List,SD_05_95List)
    x1=np.linspace(0,1200,5000)
    y1=gradient*x1+intercept
    axs[0,1].scatter(vS30List, SD_05_95List, label="SD5-95")
    axs[0,1].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Vs30" + "+" + str(round(intercept,3)))
    axs[0,1].set_ylabel("Significant Duration 5-95 (sec)")
    axs[0,1].set_xlabel("Vs30 (m/s)")
    stream_code=streamCode(path+ "/" + files[0])
    axs[0,1].legend(title=stream_code, fontsize="x-small")
    axs[0,1].grid(True)
    axs[0,1].set_title("Vs30 - Ds5-95", fontsize=20, color="red")

    gradient, intercept, r_value, p_value, std_err=linregress(vS30List, SD_05_95List1)
    x1=np.linspace(0,1200,5000)
    y1=gradient*x1+intercept
    axs[1,1].scatter(vS30List, SD_05_95List1, label="SD5-95")
    axs[1,1].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Vs30" + "+" + str(round(intercept,3)))
    axs[1,1].set_ylabel("Significant Duration 5-95 (sec)")
    axs[1,1].set_xlabel("Vs30 (m/s)")
    stream_code=streamCode(path+ "/" + files[1])
    axs[1,1].legend(title=stream_code, fontsize="x-small")
    axs[1,1].grid(True)

    gradient, intercept, r_value, p_value, std_err=linregress(vS30List,SD_05_95List2)
    x1=np.linspace(0,1200,5000)
    y1=gradient*x1+intercept
    axs[2,1].scatter(vS30List, SD_05_95List2, label="SD5-95")
    axs[2,1].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Vs30" + "+" + str(round(intercept,3)))
    axs[2,1].set_ylabel("Significant Duration 5-95 (sec)")
    axs[2,1].set_xlabel("Vs30 (m/s)")
    stream_code=streamCode(path+ "/" + files[2])
    axs[2,1].legend(title=stream_code, fontsize="x-small")
    axs[2,1].grid(True)

    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList,SD_05_95tList)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    axs[0,2].scatter(rEpiList, SD_05_95tList, label="SD5-95")
    axs[0,2].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    axs[0,2].set_ylabel("Significant Duration 5-95 (sec)")
    axs[0,2].set_xlabel("Epicenteral Distance (km)")
    axs[0,2].legend(fontsize="x-small")
    axs[0,2].grid(True)

    gradient, intercept, r_value, p_value, std_err=linregress(vS30List,SD_05_95tList)
    x1=np.linspace(0,1200,5000)
    y1=gradient*x1+intercept
    axs[1,2].scatter(vS30List, SD_05_95tList, label="SD5-95")
    axs[1,2].plot(x1,y1, '-r', label="SD5-95=" + str(round(gradient,3)) + "Vs30" + "+" + str(round(intercept,3)))
    axs[1,2].set_ylabel("Significant Duration 5-95 (sec)")
    axs[1,2].set_xlabel("Vs30 (m/s)")
    axs[1,2].legend(fontsize="x-small")
    axs[1,2].grid(True)

    # plt.show()

    # plt.savefig("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/Plots" + "/" + "R-SD5-95 and VS30-SD5-95" + ".png")


#################################################-- MODIFIED MERCALLI INTENSITY  --#################################################
####################################################################################################################################

def modifiedMercalliIntensity(path):
    files=getFiles(path)
    SList=[]
    a=5.919
    b=0.844
    r=-0.997
    s=-0.105
    vS30List=list()
    AIList=list()
    AI1List=list()
    AItList=list()
    AI2List=list()
    rEpiList=list()
    rEpiList1=list()
    MMI=[]
    for i in range(0, len(files),3):
        vS30=Vs30(path+ "/" + files[i])
        vS30List.append(vS30)
        AI=AriasIntensity(path+ "/" + files[i])
        AI=AI[-1]
        AIList.append(AI)
        AI1=AriasIntensity(path+ "/" + files[i+1])
        AI1=AI1[-1]
        AI1List.append(AI1)
        AI2=AriasIntensity(path+ "/" + files[i+2])
        AI2=AI2[-1]
        AI2=np.log10(AI2)
        AI2List.append(AI2)
        rEpi=Repi(path+ "/" + files[i])
        rEpiList1.append(rEpi)
        rEpi=np.log10(rEpi)
        rEpiList.append(rEpi)
    for j in range(0, len(vS30List)):
        if vS30List[j]<=360:
            SList.append(1)
        elif vS30List[j]>360:
            SList.append(0)

        AItList.append(np.log10((AIList[j]*AI1List[j])**0.5))
        MMI.append(a+b*AItList[j]+r*rEpiList[j]+s*SList[j])
    
    return MMI


##############################################-- MODIFIED MERCALLI INTENSITY GRAPH --###############################################
####################################################################################################################################

def MMIGraph(path):
    files=getFiles(path)
    rEpiList1=list()
    MMI=modifiedMercalliIntensity(path)
    for i in range(0, len(files),3):
        rEpi=Repi(path+ "/" + files[i])
        rEpiList1.append(rEpi)

    from scipy.stats import linregress
    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList1,MMI)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    fig, ax =plt.subplots()
    fig.suptitle("MMI - R GRAPH ", fontsize=14, color="red")
    ax.plot(rEpiList1,MMI, "yo", label="MMI")
    ax.plot(x1,y1, '-r', label="MMI=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    ax.set_xlabel("Epicenteral Distance (km)")
    ax.set_ylabel("MMI(AI)")
    ax.set_ylim(3,7)
    plt.legend()
    # plt.show()
    # plt.savefig("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/Plots" + "/" + "MMI" + ".png")
    


#################################################-- PREDICTION OF ARIAS INTENSITY --################################################
####################################################################################################################################

def predictionAI(path):
    files=getFiles(path)
    Fn=1
    Fr=0
    a1=5.48044
    a2=0.39645
    a3=-0.10684
    a4=-2.15658
    a5=0.39202
    a6=9
    a7=-0.01572
    a8=-0.09282
    a9=-1.02908
    vS30List=list()
    rEpiList=list()
    predictedAI=list()
    MwList=list()
    for i in range(0, len(files),3):
        vS30=Vs30(path+ "/" + files[i])
        if vS30>1000:
            vS30List.append(100)
        else:
            vS30List.append(vS30)

        Mw=momentMagnitude(path+ "/" + files[i])
        MwList.append(Mw)
        rEpi=Repi(path+ "/" + files[i])
        rEpiList.append(rEpi)
    
    for j in range(0, len(vS30List)):
        preAI=a1+a2*(MwList[j]-6.5)+a3*(8.5-MwList[j])**2+(a4+a5*(MwList[j]-6.5))*np.log((((rEpiList[j])**2)+(a6**2))**0.5)+a7*Fn+a8*Fr+a9*np.log(vS30List[j]/750)
        predictedAI.append(np.exp(preAI)*100)
    
    # print(predictedAI)
    return predictedAI 


###########################################-- PREDICTION OF CUMULATIVE ABSOLUTE VELOCITY --#########################################
####################################################################################################################################

def predictionCAV(path):
    files=getFiles(path)
    Fn=1
    Fr=0
    a1=9.01362
    a2=0.7516
    a3=-0.03713
    a4=-0.93573
    a5=0.10417
    a6=9
    a7=0.02738
    a8=-0.16629
    a9=-0.644
    vS30List=list()
    rEpiList=list()
    predictedCAV=list()
    MwList=list()
    for i in range(0, len(files),3):
        vS30=Vs30(path+ "/" + files[i])
        if vS30>1000:
            vS30List.append(100)
        else:
            vS30List.append(vS30)

        Mw=momentMagnitude(path+ "/" + files[i])
        MwList.append(Mw)
        rEpi=Repi(path+ "/" + files[i])
        rEpiList.append(rEpi)
    
    for j in range(0, len(vS30List)):
        preCAV=a1+a2*(MwList[j]-6.5)+a3*(8.5-MwList[j])**2+(a4+a5*(MwList[j]-6.5))*np.log((((rEpiList[j])**2)+(a6**2))**0.5)+a7*Fn+a8*Fr+a9*np.log(vS30List[j]/750)
        predictedCAV.append(np.exp(preCAV))
    
    # print(predictedCAV)
    return predictedCAV 


############################-- FIND MODIFIED MERCALLI INTENSITY BY USING PREDICTED ARIAS INTENSITY --###############################
####################################################################################################################################

def predictedModifiedMercalliIntensity(path):
    files=getFiles(path)
    SList=[]
    a=5.919
    b=0.844
    r=-0.997
    s=-0.105
    vS30List=list()
    AIpreList=list()
    rEpiList=list()
    rEpiList1=list()
    predictedMMI=[]
    for i in range(0, len(files),3):
        rEpi=Repi(path+ "/" + files[i])
        rEpiList1.append(rEpi)
        rEpi=np.log10(rEpi)
        rEpiList.append(rEpi)
        vS30=Vs30(path+ "/" + files[i])
        vS30List.append(vS30)
    for j in range(0, len(rEpiList1)):
        AIpreList=predictionAI(path)
        if vS30List[j]<=360:
            SList.append(1)
        elif vS30List[j]>360:
            SList.append(0)

        predictedMMI.append(a+b*np.log10(AIpreList[j])+r*rEpiList[j]+s*SList[j])

    return predictedMMI
    

###########################-- MODIFIED MERCALLI INTENSITY GRAPH BY USING PREDICTED ARIAS INTENSITY --###############################
####################################################################################################################################
    
def predictedMMIGraph(path):
    files=getFiles(path)
    rEpiList1=list()
    predictedMMI=predictedModifiedMercalliIntensity(path)
    for i in range(0, len(files),3):
        rEpi=Repi(path+ "/" + files[i])
        rEpiList1.append(rEpi)
    
    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList1,predictedMMI)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    fig, ax =plt.subplots()
    fig.suptitle("Predicted MMI - R GRAPH ", fontsize=14, color="red")
    ax.plot(rEpiList1, predictedMMI, "yo", label="MMI")
    ax.plot(x1,y1, '-r', label="MMI=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    ax.set_xlabel("Epicenteral Distance (km)")
    ax.set_ylabel("Predicted MMI(AI)")
    ax.set_ylim(3,7)
    plt.legend()
    # plt.show()
    # plt.savefig("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/Plots" + "/" + "PredictedMMI" + ".png")

    


################################################-- AI AND CAV RESIDUAL VALUES --####################################################
####################################################################################################################################

def residuals(path):
    files=getFiles(path)
    predictedAI=predictionAI(path)
    predictedCAV=predictionCAV(path)
    AIList=list()
    AI1List=list()
    AItList=list()
    CAVList=[]
    CAV1List=[]
    CAVtList=[]
    residualsAIList=[]
    residualsCAVList=[]
    for i in range(0, len(files),3):
        AI=AriasIntensity(path+ "/" + files[i])
        AI=AI[-1]
        AIList.append(AI)
        AI1=AriasIntensity(path+ "/" + files[i+1])
        AI1=AI1[-1]
        AI1List.append(AI1)
        CAV=CumAbsVel(path+ "/" + files[i])
        CAV=CAV[-1]
        CAVList.append(CAV)
        CAV1=CumAbsVel(path+ "/" + files[i+1])
        CAV1=CAV1[-1]
        CAV1List.append(CAV)
    for j in range(0, len(AI1List)):
        AItList.append((AIList[j]*AI1List[j])**0.5)
        residualAI=np.log(AItList[j]/predictedAI[j])
        residualsAIList.append(residualAI)
        CAVtList.append((CAVList[j]*CAV1List[j])**0.5)
        residualCAV=np.log(CAVtList[j]/predictedCAV[j])
        residualsCAVList.append(residualCAV)
    
    return [residualsAIList, residualsCAVList]


###########################################-- R-RESIDUAL AND VS30-RESIDUAL GRAPHS --################################################
####################################################################################################################################

def residualsGraph(path):
    files=getFiles(path)
    rEpiList=list()
    vS30List=list()
    for i in range(0, len(files),3):
        rEpi=Repi(path+ "/" + files[i])
        rEpiList.append(rEpi)
        vS30=Vs30(path+ "/" + files[i])
        vS30List.append(vS30)

    [residualAI, residualCAV]=residuals(path)

    fig, axs = plt.subplots(2,2,figsize=(20,10))
    fig.subplots_adjust(hspace=0.23, wspace=0.19, left=0.04, right=0.98, bottom=0.06, top=0.95)

    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList,residualAI)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    axs[0,0].scatter(rEpiList, residualAI, label="Residual AI")
    axs[0,0].plot(x1,y1, '-r', label="Residual AI=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    axs[0,0].set_xlabel("Epicenteral Distance (km)")
    axs[0,0].set_ylabel("Residual AI")
    axs[0,0].set_ylim(-3,3)
    axs[0,0].set_title("Residual AI Graphs", fontsize=12, color="red")
    axs[0,0].legend()

    gradient, intercept, r_value, p_value, std_err=linregress(vS30List,residualAI)
    x1=np.linspace(0,1200,5000)
    y1=gradient*x1+intercept
    axs[1,0].scatter(vS30List, residualAI, label="Residual AI")
    axs[1,0].plot(x1,y1, '-r', label="Residual AI=" + str(round(gradient,3)) + "Vs30" + "+" + str(round(intercept,3)))
    axs[1,0].set_xlabel("Vs30 (cm/s)")
    axs[1,0].set_ylabel("Residual AI")
    axs[1,0].set_ylim(-3,2)
    axs[1,0].legend()

    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList,residualCAV)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    axs[0,1].scatter(rEpiList, residualCAV, label="Residual CAV")
    axs[0,1].plot(x1,y1, '-r', label="Residual CAV=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    axs[0,1].set_xlabel("Epicenteral Distance (km)")
    axs[0,1].set_ylabel("Residual CAV")
    axs[0,1].set_ylim(-2,1)
    axs[0,1].set_title("Residual CAV Graphs", fontsize=12, color="red")
    axs[0,1].legend()

    gradient, intercept, r_value, p_value, std_err=linregress(vS30List,residualCAV)
    x1=np.linspace(0,1200,5000)
    y1=gradient*x1+intercept
    axs[1,1].scatter(vS30List, residualCAV, label="Residual CAV")
    axs[1,1].plot(x1,y1, '-r', label="Residual CAV=" + str(round(gradient,3)) + "Vs30" + "+" + str(round(intercept,3)))
    axs[1,1].set_xlabel("Vs30 (cm/s)")
    axs[1,1].set_ylabel("Residual CAV")
    axs[1,1].set_ylim(-2,1)
    axs[1,1].legend()
    # plt.show()
    # plt.savefig("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/Plots" + "/" + "Residual" + ".png")


#############################################-- INSTRUMENTAL INTENSITY RESIDUAL --##################################################
####################################################################################################################################

def residualIntensity(path):
    files=getFiles(path)
    rEpiList=list()
    vS30List=list()
    recInt=modifiedMercalliIntensity(path)
    predInt=predictedModifiedMercalliIntensity(path)
    residualMMIList=[]
    for i in range(0, len(files),3):
        rEpi=Repi(path+ "/" + files[i])
        rEpiList.append(rEpi)
        vS30=Vs30(path+ "/" + files[i])
        vS30List.append(vS30)
    
    for j in range(0, len(rEpiList)):
        residualMMI=np.log(recInt[j]/predInt[j])
        residualMMIList.append(residualMMI)
        
    
    fig, (ax1, ax2)=plt.subplots(2,figsize=(10,10))
    fig.suptitle("Residual Intensity ", fontsize=14, color="red")

    gradient, intercept, r_value, p_value, std_err=linregress(rEpiList,residualMMIList)
    x1=np.linspace(0,100,500)
    y1=gradient*x1+intercept
    ax1.scatter(rEpiList, residualMMIList, label="Residual Intensity")
    ax1.plot(x1,y1, '-r', label="Residual Intensity=" + str(round(gradient,3)) + "Repi" + "+" + str(round(intercept,3)))
    ax1.set_xlabel("Epicenteral Distance (km)")
    ax1.set_ylabel("Residual Intensity")
    ax1.set_ylim(-1,1)
    ax1.legend()

    gradient, intercept, r_value, p_value, std_err=linregress(vS30List,residualMMIList)
    x1=np.linspace(0,1200,5000)
    y1=gradient*x1+intercept
    ax2.scatter(vS30List, residualMMIList, label="Residual Intensity")
    ax2.plot(x1,y1, '-r', label="Residual Intensity=" + str(round(gradient,3)) + "Vs30" + "+" + str(round(intercept,3)))
    ax2.set_xlabel("Vs30 (cm/s)")
    ax2.set_ylabel("Residual Intensity")
    ax2.set_ylim(-1,1)
    ax2.legend()

    # plt.show()
    # plt.savefig("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/Plots" + "/" + "ResidualIntensity" + ".png")


    


# residualIntensity("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data")
# residualsGraph("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data")
# plotSigDur("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data")
# modifiedMercalliIntensity("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data")
# BaselineCorrection("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data/GuestUser_RawAcc_20201030115124_0905/20201030115124_0905_mp_RawAcc_E.asc")
# plotData("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data")    
# Repi("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data/GuestUser_RawAcc_20201030115124_0905/20201030115124_0905_mp_RawAcc_N.asc")
# print(AriasIntensity("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data/GuestUser_RawAcc_20201030115124_0905/20201030115124_0905_mp_RawAcc_E.asc")[-1])






