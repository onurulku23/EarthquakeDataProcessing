from scipy import signal
import numpy as np

ag_txt=np.loadtxt("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data/GuestUser_RawAcc_20201030115124_0905/20201030115124_0905_mp_RawAcc_E.asc", skiprows=64)
groundacc=ag_txt[:]
ags=groundacc.flatten("C")
t_amount=len(ags)
dt=0.01
dt=float(dt)
initial_time=0
end_time=initial_time+dt*t_amount
t=np.arange(initial_time, end_time, dt)
print(np.mean(ags))
detrend_sdata=signal.detrend(ags)
print(np.mean(detrend_sdata))