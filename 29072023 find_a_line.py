'''
Created by Nova Foster (They/Them) 26/07/2023
'''

import pandas as pd                                               #Use pandas to read csv file

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
  import numpy as np
  Boltzmann = 8.617333262e-5

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
  import matplotlib.pyplot as py
  import numpy as np
  from scipy.optimize import curve_fit

  Temp=-1
  #Load values for matched lines into seperate array for easier handling
  lines = NIST_Data[NIST_Data["obs_wl_air(nm)"].isin(matched)]

  #Calculate x and y values for the plot & plot them
  x_data = lines.iloc[:]['Ek(eV)']
  y_data = np.log( (intensity*lines.iloc[:]['obs_wl_air(nm)'])/(lines.iloc[:]['g_k']*lines.iloc[:]['Aki(10^8 s^-1)']))

  #Format the data
  x_data, y_data = zip(*sorted(zip(x_data, y_data)))   #Combine, sort then split x and y data so they are plotted correctly. This does convert the dtype from dframe to tuple
  x_data = np.asarray(x_data,dtype=float)   #Convert the two to numpy arrays as floats
  y_data = np.asarray(y_data,dtype=float)
  py.plot(x_data,y_data,".")

  
  #Fit a line to x and y data. Converting to numpy arrays as pandas can't be used
  #Pop[0] is gradient of line, pop[1] is y intercept
  #pcov covaraince so np.sqrt(pcov.diagonal()[:]) is uncertainty
  pop, pcov = curve_fit(line,x_data,y_data)
  print(pop)
  print(np.sqrt(pcov.diagonal()[:]))


  #Plot the fit
  x_min = np.min(x_data)
  x_max = np.max(x_data)

  x_range = np.linspace(float(x_min),float(x_max),1000)
  y_fit_vals = line(x_range,pop[0],pop[1])
  py.plot(x_range,y_fit_vals,"r-")
  py.show()


  # Calculate temperature from the slope
  k = 8.617333262145 * 10**-5  # Boltzmann constant in eV/K
  slope = pop[0]
  weighted_avg_energy_diff = -1 * slope  # Negative because of the formula
  Temp = -weighted_avg_energy_diff / (k * slope)

  return Temp


def boltz_fit(matched,intensity):
  '''
  TODO:
  - Implement this using curvefit
  - Test using Boltzmann dist. with noise and added lines
  - Look at various windows to see what works
  '''

  return Temp

import numpy as np
lambdas = [453.0785,453.9695,510.5541,515.3235,521.8202,522.0070,529.2517]
intens = np.array([800,800,1500,2000,1650,2500,1650])
intens = intens / 2500


print(intens)

print(boltz_line(lambdas,intens))

'''
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