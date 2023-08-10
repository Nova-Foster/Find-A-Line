'''
Created by Nova Foster (They/Them) 26/07/2023
'''

import pandas as pd                                               #Use pandas to read csv file
import matplotlib.pyplot as py
import numpy as np
Boltzmann = 8.617333262e-5

'''
NIST Values reference: Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2022). NIST Atomic Spectra Database (ver. 5.10), [Online]. Available: https://physics.nist.gov/asd [2023, July 29]. National Institute of Standards and Technology, Gaithersburg, MD. DOI: https://doi.org/10.18434/T4W30F
NIST File contains:
- Element
- Spectral Number
- Observed wavelength in air (nm) from 250nm to 950nm
- Observed wv. uncertainty
- Intensity
- Probability in terms of 10^8 s^-1 (closer to 1 is more likely)
- Accuracy
- Initial energy (eV)
- Final energy (eV)
- Statistical weighting of intial state
- Statisitcal weighting of excited state
- Type (mostly NA)
- TP Reference
- Line Reference
'''
NIST_Data = pd.read_csv("Lines_eV.csv")

#List of elements to be used in compare, not done in function for efficiency
main_air_elements = {"H", "He", "Ar", "N", "O"}
all_air_elements = {"H", "He", "Ar", "N", "O","Ne", "Kr", "Xe", "I"}

def compare(Observed,margin=0.5,source="all"):                    #Compare observed wavelength (nm) to NIST values within range +-margin. Source being elements to select from (Main Air, All Air or All)
 #Change format of inputs
  observed_float = float(Observed)
  margin_float = float(margin)
  source=source.lower()

 #Select all wavelength from NIST that fit the range
  lines = NIST_Data[NIST_Data['obs_wl_air(nm)'].between(observed_float-margin_float,observed_float+margin_float)]

 #Select the soures from that range
  #Main Air: H, He, Ar, N, O
  if(source=="main air"):
    lines = lines[ lines['element'].isin(main_air_elements)]
  #All Air: H, He, Ar, N, O, Ne, Kr, Xe, I
  elif(source=="all air"):
    lines = lines = lines[ lines['element'].isin(all_air_elements)]
  #Else covers all possible sources

  return lines


def temp_using_2line(matched,intensity):
  '''
  TODO:
  - Loop so each pair of lines is used (If the elements match)
  '''

 #Select values for each of the two lines and assign them to a new structure. iloc to get them as numerical values not data frames
  line1 = NIST_Data[NIST_Data["obs_wl_air(nm)"]==matched[0]].iloc[0]
  line2 = NIST_Data[NIST_Data["obs_wl_air(nm)"]==matched[1]].iloc[0]

 #Calculate each part of the equation
  prefactor = (float(line2['Ek(eV)']) - float(line1['Ek(eV)'])) /Boltzmann
  numerator = intensity[0]*line1['obs_wl_air(nm)']*line2['Aki(10^8 s^-1)']*line2['g_k']
  denomonator = intensity[1]*line2['obs_wl_air(nm)']*line1['Aki(10^8 s^-1)']*line1['g_k']

 #Calculate temperature using the equation
  Temp = prefactor * np.log(numerator/denomonator)**(-1)

  return Temp

def line(x,m,c):   #Line function used for curve_fit in bolt_line
  return m*x+c

def boltz_line(matched,intensity):
  from scipy.optimize import curve_fit

  #Load values for matched lines into seperate array for easier handling
  lines = NIST_Data[NIST_Data["obs_wl_air(nm)"].isin(matched)]

  #Calculate x and y values for the plot
  x_data = lines['Ek(eV)']
  y_data = np.log( (intensity*lines['obs_wl_air(nm)'])/(lines['g_k']*lines['Aki(10^8 s^-1)']))

  #Format the data
  x_data, y_data = zip(*sorted(zip(x_data, y_data)))   #Combine, sort then split x and y data so they are plotted correctly. This does convert the dtype from dframe to tuple
  x_data = np.asarray(x_data,dtype=float)   #Convert the two to numpy arrays as floats
  y_data = np.asarray(y_data,dtype=float)
  py.plot(x_data,y_data,".",label="Data")

  
  #Fit a line to x and y data. Converting to numpy arrays as pandas can't be used
  #Pop[0] is gradient of line, pop[1] is y intercept
  #pcov covaraince so np.sqrt(pcov.diagonal()[:]) is uncertainty
  pop, pcov = curve_fit(line,x_data,y_data)

  #Plot the fit
  x_range = np.linspace(np.min(x_data), np.max(x_data), 1000)
  py.plot(x_range, line(x_range, pop[0], pop[1]), "r-", label="Fit")

  py.legend()
  py.show()


  # Calculate temperature from the slope
  Temp = (-1/pop[0])/Boltzmann


  return Temp


def boltz_fit(matched,intensity):
  '''
  TODO:
  - Implement this using curvefit
  - Test using Boltzmann dist. with noise and added lines
  - Look at various windows to see what works
  '''

  return Temp


def convert_imd(file_name):
  with open(file_name, "rb") as file:      #Open the file reading as non-text
    if file:
        #Read first 2 bytes, least important byte first, covert to int from hex
        VER = int.from_bytes(file.read(2), byteorder='little', signed=False) 
        W = int.from_bytes(file.read(2), byteorder='little', signed=False)
        H = int.from_bytes(file.read(2), byteorder='little', signed=False)

        #Read the rest of the data
        m_s32data = np.fromfile(file, dtype=np.int32, count=W * H)

        #Conver the data into floats, with error handling for 0
        fdata = np.where(m_s32data != 0, m_s32data / 1000.0, 0.0)


  #Create 2d array
  rows = np.split(fdata,1040,axis=0)                                         #Split the data into each 1040 pixel row
  data_2d_array = np.stack(rows,axis=0)                                      #Stack each row on top of eachother

  return data_2d_array

def spark_plot(data_2d, wavelengths=0):

  #Colour plot of whole image
  py.imshow(data_2d,interpolation="nearest",origin="lower")
  py.colorbar()
  py.show()


  #Horizontal profile
  hori = np.sum(data_2d,axis=0)
  py.plot(np.arange(0,len(hori),1),hori)
  py.show()

  #Verticle profile
  vert = np.sum(data_2d,axis=1)
  py.plot(np.arange(0,len(vert),1),vert)
  py.show()
  return


def auto_peaks(data_2d,rel_prom=0.01,rel_width=0.15):
  import scipy.signal as sg

  #Generate horizontal profile
  hori = np.sum(data_2d,axis=0)
  py.plot(np.arange(0,len(hori),1),hori)
  
  #Use scipy to find peaks based on local maxima
  peak_x_cords,_ = sg.find_peaks(hori)                                       #Generate X coords

  peak_prominence,_,_ = sg.peak_prominences(hori,peak_x_cords)               #Find how prominent each peak is based on surroundings
  prominent_x = peak_x_cords[peak_prominence>rel_prom*max(peak_prominence)]
  
 
  peak_width,_,_,_ = sg.peak_widths(hori,prominent_x)
  Selected_peaks = prominent_x[peak_width>rel_width*max(peak_width)]

  #Draw a line for each prominent peak
  for i in range(len(Selected_peaks)):
    py.plot( np.linspace(Selected_peaks[i],Selected_peaks[i],1000), np.linspace(min(hori),max(hori),1000), "r-",alpha=0.2)

  py.show()
  return

data = convert_imd( "800V - 250Hz - 500usmm - 15um vertical.imd")
#spark_plot( data )
auto_peaks( data )

'''
lambdas = [453.0785,453.9695,510.5541,515.3235,521.8202,522.0070,529.2517]
intens = np.array([800,800,1500,2000,1650,2500,1650])
intens = intens / 2500


print(intens)

print(boltz_line(lambdas,intens))


print(temp_using_2line([510.5541,515.3235],[0.55,0.9]))
'''
'''
while True:
  print("Please enter the observed wavelength(nm), margin to search (observed +-) and lines to check (Main Air, All Air or all)")
  observed = input()
  margin = input()
  source = input()

  print( compare(observed, margin, source) )
  '''