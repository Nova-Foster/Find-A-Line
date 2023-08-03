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

def compare(Observed,margin=0.5,source="all"):                    #Compare observed wavelength (nm) to NIST values within range +-margin. Source being elements to select from (Main Air, All Air or All)
  source=source.lower()
#Select all wavelength from NIST that fit the range
  lines = NIST_Data[NIST_Data['obs_wl_air(nm)'].between(float(Observed)-float(margin),float(Observed)+float(margin))]

#Select the soures from that range
  #Main Air: H, He, Ar, N, O
  if(source=="main air"):
    lines = lines[ (lines['element'] =="H") | (lines['element'] =="He") | (lines['element'] =="Ar") | (lines['element'] =="N") | (lines['element'] =="O")]
  #All Air: H, He, Ar, N, O, Ne, Kr, Xe, I
  elif(source=="all air"):
    lines = lines[ (lines['element'] =="H") | (lines['element'] =="He") | (lines['element'] =="Ar") | (lines['element'] =="N") | (lines['element'] =="O")
    | (lines['element'] =="Ne") | (lines['element'] =="Kr") | (lines['element'] =="Xe")| (lines['element'] =="I")]
  #Else covers all possible sources

  return lines

def temp_using_2line(matched,intensity):
  '''
  Todo:
  - Add actual function
  - Comment
  - Check that it works using values in paper
  '''
  import numpy as np
  Boltzmann = 8.617333262e-5

#Select values for each of the two lines and assign them to a new structure
  line1 = NIST_Data[NIST_Data["obs_wl_air(nm)"]==matched[0]]
  line2 = NIST_Data[NIST_Data["obs_wl_air(nm)"]==matched[1]]

 
#Select values using iloc[0][THING] to select single values
  prefactor = (float(line2.iloc[0]['Ek(eV)']) - float(line1.iloc[0]['Ek(eV)'])) /Boltzmann
  numerator = intensity[0]*line1.iloc[0]['obs_wl_air(nm)']*line2.iloc[0]['Aki(10^8 s^-1)']*line2.iloc[0]['g_k']
  denomonator = intensity[1]*line2.iloc[0]['obs_wl_air(nm)']*line1.iloc[0]['Aki(10^8 s^-1)']*line1.iloc[0]['g_k']

  print(line2 )


  Temp = prefactor * np.log(numerator/denomonator)**(-1)

  return Temp


def boltz_fit():
  '''
  Todo:
  - Get lines selected
  - Look at paper to see relation
  - Curvefit a BZ distr. to those lines
  - Estimate temperature

  '''
  return Temp

'''
print(temp_using_2line([510.5541,515.3235],[0.55,0.9]))

'''
while True:
  print("Please enter the observed wavelength(nm), margin to search (observed +-) and lines to check (Main Air, All Air or all)")
  observed = input()
  margin = input()
  source = input()

  print( compare(observed, margin, source) )