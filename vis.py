import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns

hours_series = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\time.dat")

exp_coszen_series = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\exp_coszen.dat")
ground_T_series = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\ground_T.dat")
canopy_T_series = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\canopy_T.dat")

Isoprene_emi_ts = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\emi_iso.dat")
Apenene_emi_ts = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\emi_alp.dat")

OH_data = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\Conc_OH.dat")
H2O_data = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\Conc_H20.dat")
H2SO4_data = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\Conc_H2SO4.dat")
ELVOC_data = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\Conc_ELVOC.dat")

ApineneConc_data = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\emi_alp_Conc.dat")
IsopreneConc_data = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\emi_iso_Conc.dat")

O3_Conc = np.genfromtxt("C:\\Users\\stranger\\Documents\\Exercises\\Day4\\MeteoChem\\Output\\O3_Conc.dat")

#plotting exp_coszen, ground and canopy's T
plt.figure()
plt.subplot(3,1,1, )
plt.plot(exp_coszen_series)
plt.tick_params(labelbottom=False)
plt.grid()
plt.ylabel("exp_coszen")
plt.subplot(3,1,2)
plt.plot(canopy_T_series)
plt.tick_params(labelbottom=False)
plt.grid()
plt.ylabel("Canopy T")
plt.subplot(3,1,3)
plt.plot(ground_T_series)
plt.grid()
plt.xlabel("Day")
plt.ylabel("Ground T")
plt.suptitle("Exp_coszen, canopy's and ground's T")


#plotting emission's timeseries
plt.figure()
plt.plot(hours_series,Isoprene_emi_ts, label="isoprene")
plt.plot(hours_series,Apenene_emi_ts, label="a-pinene")
plt.xlim(xmin=3.0)
plt.title("Isoprene and a-pinene time series")
plt.xlabel("Day")
plt.ylabel("Emissions rate")
plt.grid()
plt.legend()

#plot gases
plt.figure()
plt.subplot(2,2,1)
plt.title("OH")
plt.plot(OH_data[71:,5], label="50 m")
plt.plot(OH_data[71:,1], label="10 m")
plt.grid()
plt.legend()

plt.subplot(2,2,2)
plt.title("HO2")
plt.plot(H2O_data[71:,5], label="50 m")
plt.plot(H2O_data[71:,1], label="10 m")
plt.grid()
plt.legend()

plt.subplot(2,2,3)
plt.title("H2SO4")
plt.plot(H2SO4_data[71:,5], label="50 m")
plt.plot(H2SO4_data[71:,1], label="10 m")
plt.grid()
plt.legend()

plt.subplot(2,2,4)
plt.title("ELVOC")
plt.plot(ELVOC_data[71:,5], label="50 m")
plt.plot(ELVOC_data[71:,1], label="10 m")
plt.grid()
plt.legend()
plt.suptitle("Smth here")
plt.subplots_adjust(hspace=0.4,wspace=0.3)

#Gas time-series
plt.figure()
plt.subplot(1,2,1)
plt.title("Apinene_Conc")
plt.plot(ApineneConc_data[71:,5], label="50 m")
plt.plot(ApineneConc_data[71:,1], label="10 m")
plt.grid()
plt.legend()
plt.subplot(1,2,2)
plt.title("Isoprene_Conc")
plt.plot(IsopreneConc_data[71:,5], label="50 m")
plt.plot(IsopreneConc_data[71:,1], label="10 m")
plt.grid()
plt.legend()


#Ozone plot
plt.figure()
plt.title("Ozone")
plt.plot(O3_Conc[71:, 5], label="50M")
plt.plot(O3_Conc[71:, 1], label="10M")
plt.grid()
plt.legend()


#SHOW THEM ALL
plt.show()