'''
Created by Nova Foster (They/Them) 26/07/2023
'''


import numpy as np
import pandas as pd


NIST_Data = pd.read_csv("Lines.csv")

def compare(Observed,margin=0.5,source="All"):                    #Compare observed wavelength (nm) to NIST values within range +-margin. Source being elements to select from (Main Air, All Air or All)

#Select all wavelength from NIST that fit the range
  lines = NIST_Data[NIST_Data['obs_wl_air(nm)'].between(Observed-margin,Observed+margin)]

#Select the soures from that range
  #Main Air: H, He, Ar, N, O
  if(source=="Main Air"):
    lines = lines[ (lines['element'] =="H") | (lines['element'] =="He") | (lines['element'] =="Ar") | (lines['element'] =="N") | (lines['element'] =="O")]
  #All Air: H, He, Ar, N, O, Ne, Kr, Xe, I
  elif(source=="All Air"):
    lines = lines[ (lines['element'] =="H") | (lines['element'] =="He") | (lines['element'] =="Ar") | (lines['element'] =="N") | (lines['element'] =="O")
    | (lines['element'] =="Ne") | (lines['element'] =="Kr") | (lines['element'] =="Xe")| (lines['element'] =="I")]
  #Else covers all possible sources

  return lines

print(compare(712.219,5,"Main Air"))
