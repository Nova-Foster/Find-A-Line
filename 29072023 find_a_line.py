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

def compare(Observed,margin=0.5,source="main air"):                    #Compare observed wavelength (nm) to NIST values within range +-margin. Source being elements to select from (Main Air, All Air or All)
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

def convert_imd(file_name):                #Create 2d array from imd: 1392 columns of wavelength & 1040 rows for time
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

  data_2d_array = np.flip(data_2d_array,0)      #Flip data so time 0 is at index 0
  return data_2d_array

def cont_wavelengths(calibration):         #Determine the wavelength value for each pixel based on the calibration file
#Working based on calibration being 2d array: 0 = pixel, 1 = wavelength

  #Find wavelength values per pixel between calibration point 1 and 2
  known_points = np.size(calibration,1)
  wavelength_whole = np.zeros(known_points)

  current_wavelength = calibration[1][0]
  next_wavelength = calibration[1][1]
  pixel_distance = int(calibration[0][1] - calibration[0][0])

  #Find linearity between first two measurements
  temp_scale = np.linspace(current_wavelength,next_wavelength,pixel_distance)  

  #Apply to pixel before first measurement: assumes same lineraity as between point 1 and 2
  wavelength_change = temp_scale[1] - temp_scale[0]
  start_val = current_wavelength - (wavelength_change*calibration[0][0])
  start_scale = np.linspace(start_val,current_wavelength,int(calibration[0][0]))

  #Create array to be output: currently has wavelength values from pixel 0 to calibration pixel 2
  wavelength_whole = np.append(start_scale,temp_scale)

  #Loop for each calibration value
  for i in range(known_points):
    try:                                                                                   #In a try statement as it will try to access an index past the limit
      current_wavelength = calibration[1][i+1]                                             #Current wavelength is +1 as first calibration value has already been used
      next_wavelength = calibration[1][i+2]
      pixel_distance = calibration[0][i+2] - calibration[0][i+1]                           #Find distance in pixels between the two to be checked
      temp_scale = np.linspace(current_wavelength,next_wavelength,int(pixel_distance))     #Create the region for current wavelength to next

      wavelength_whole = np.append(wavelength_whole,temp_scale)                            #Append that to what has already been done
    except:
      break                                                                                #Exit the for loop, should be at the last calibrated value                                       

  #Continue the scale for last value to the final pixel
  wavelength_change = wavelength_whole[-1] - wavelength_whole[-2]                          #Amount of wavelength per pixel
  end_val = next_wavelength + (wavelength_change*(1392-calibration[0][-1]))                       #Final value of the scale
  end_scale = np.linspace(calibration[1][-1],end_val,1392-int(calibration[0][-1]))       #Create scale between final calibration and the end
  wavelength_whole = np.append(wavelength_whole,end_scale)                                 #Append to pre-existing scale
  

  py.plot(np.linspace(0,1391,1392),wavelength_whole,lw=1)
  py.xlabel("Pixel")
  py.ylabel("Wavelength(nm)")
  py.show()
  return wavelength_whole

def add_calibration(data,wavelengths=0,time=0):
  if wavelengths.any() !=0:
    values = np.vstack([wavelengths,data])
  return values

def background_subtraction(data,center=0,grating=0):
  back_ground = convert_imd("20230810 Full sweep grating 3\Sequence test 1\center 300 Dark.imd")
  subtracted = data - back_ground
  return subtracted


def spark_plot(data_2d):
  data = data_2d[1:,:]
  wavelengths = data_2d[0]
  #Colour plot of whole image
  im_show = py.figure()
  py.imshow(data,interpolation="nearest",origin="lower")
  py.xlabel("x pixel")
  py.ylabel("y pixel")
  py.colorbar()
  py.show(block=False)


  #Horizontal profile
  horizontal = py.figure()
  hori = np.sum(data,axis=0)
  py.plot(wavelengths,hori,".")
  py.show(block=False)

  #Verticle profile
  verticle = py.figure()
  vert = np.sum(data,axis=1)
  py.plot(np.arange(0,len(vert),1),vert)
  py.show(block=False)
  return


def auto_peaks(data_2d,strict=True):
  import scipy.signal as sg

  #Load prom. and width. retrctions based on strict: all done by eye so not perfect
  if strict ==False:
    rel_prom=0.075
    rel_width=0.2
  else:
    rel_prom = 0.15
    rel_width=0.25

  data = data_2d[1:,:]
  wavelengths = data_2d[0]

  #Generate horizontal profile
  hori = np.sum(data,axis=0)
  auto = py.figure()
  py.plot(wavelengths,hori)
  
  #Use scipy to find peaks based on local maxima
  peak_x_cords,_ = sg.find_peaks(hori)                                       #Generate X coords

  peak_prominence,_,_ = sg.peak_prominences(hori,peak_x_cords)               #Find how prominent each peak is based on surroundings
  prominent_x = peak_x_cords[peak_prominence>rel_prom*max(peak_prominence)]
  
 
  peak_width,_,_,_ = sg.peak_widths(hori,prominent_x)
  Selected_peaks = prominent_x[peak_width>rel_width*max(peak_width)]


  returned_val = np.zeros((len(Selected_peaks),3))
  for i in range(len(Selected_peaks)):
    current_wavelength = wavelengths[Selected_peaks[i]]

    #Format data to be returned
    returned_val[i][0] = Selected_peaks[i]
    returned_val[i][1] = current_wavelength
    returned_val[i][2] = hori[wavelengths==current_wavelength]


    #Draw a line at each wavelength
    py.plot( np.linspace(current_wavelength,current_wavelength,1000), np.linspace(min(hori),max(hori)*1.1,1000), "r-",alpha=0.2)

  py.show(block=False)

  #Format data to be returned: index: wavelength

  return returned_val

def plot_3d(intensity):

  data = intensity[1:,]
  x_coords = intensity[0]
  y_coords = np.linspace(0,5*0.718,1040)
  x_mesh, y_mesh = np.meshgrid(x_coords, y_coords)
  threed_plot = py.figure()
  ax = py.axes(projection="3d")

  ax.plot_surface(x_mesh,y_mesh,data,cmap="turbo")
  ax.set_xlabel("Wavelength (nm)")
  ax.set_ylabel("Time (ms)")
  ax.set_zlabel("Intensity")
  py.show(block=False)
  return 0

def remove_zeros(data):
  removed = np.clip(data,a_min=0,a_max=1e9)
  return removed

def integrated_line_image(data):
  import matplotlib as mpl
  import matplotlib.cm as cm
  wavelengths = data[:1,]
  intensity = data[1:,]
  hori = np.sum(intensity,axis=0)
  integrated = py.figure()

  #Convert relative data to colours
  relative_data = hori/np.amax(hori)

  for i in range(1392):
      py.plot( np.linspace(wavelengths[0][i],wavelengths[0][i],1000), np.linspace(0,1,1000),color=str(relative_data[i]))
  
  py.xlabel("Wavelength (nm)")
  py.ylabel("Relative intensity (a.u.)")
  py.show(block=False)
  return 0

test = np.array([[414,600,843,996,1010,1081,1098,1113,1119,1152,1167,1170,1215,1234,1256,1284,1293,1328],[253.652,300,365.015,404.656,407.783,427.397,431.958,435.833,437.612,446.369,450.235,452.186,462.42,469.804,473.415,479.262,480.702,491.651]])
cont_test = cont_wavelengths(test)
print(cont_test[0])
print(cont_test[-1])


data = convert_imd( "this one.IMD")
data_subd = background_subtraction(data,300,3)
values = add_calibration(data_subd,cont_test)
values_to_zeros = remove_zeros(values)

spark_plot(values_to_zeros)
plot_3d(values_to_zeros)
integrated_line_image(values_to_zeros)

temp = auto_peaks(values_to_zeros)


sha = np.shape(temp)
for i in range(sha[0]):
  print(compare(temp[i][1],1))

input()


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


'''
TODO
- BB plot
- Make boltz line work
- add axis for time into the 2d array   <- need slope values for each speed
# 3d plot across whole thing
# image for each wavelength across image with relative intensities
- stitch the two images together
- Make autopeaks record intensity and wavelength
- adjust temp calcs so it uses the recorded wavelengths and the known
- plug auto peaks into temp calcs
- plot of one wavelengths as it progresses across time, do this based on autopeaks?

- Open scaling file, Can't figure out file format   <- not really needed as long as the format is input correctly
'''

'''
b'\x02\x00X\x02\x02nmX\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xe0?\x0c\x00\x9e\x01\x9e\x01K\x03\xbe\x9f\x1a/\xdd\xb4o@K\x03\xe4\x03\n'
b'\x02\x00X\x02\x02nmX\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xe0?\x17\x00J\x00J\x00\xdd\x00\n'
b'\x02\x00X\x02\x02nmX\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xe0?\x10\x00\x1e\x00\x1e\x00i\x00\xac\x1cZd;h\x81@i\x00j\x00H\xe1z\x14\xae\x07\x82@j\x00\x91\x00J\x0c\x02+\x87\x18\x82@\x91\x001\x02T\xe3\xa5\x9b\xc4X\x82@1\x02X\x02\xd3Mb\x10X\xc4\x85@X\x02\xbb\x02\xb2\x9d\xef\xa7\xc6\x15\x86@\xbb\x02\x00\x03=\n'

Would make sense to be wavlengths as 500 center is much longer due to having more lines & 900 being more than 300 due to all values being bigger
b'\xd7\xa3p=\xd0v@\xe4\x03\xf2\x03\xd1"\xdb\xf9~Jy@\xf2\x03Y\x04J\x0c\x02+\x87|y@Y\x04\x92\x04\x17\xd9\xce\xf7S={@\x92\x04\xbf\x04\xe5\xd0"\xdb\xf9B|@\xbf\x04\xd2\x04\x1f\x85\xebQ\xb8\xe6|@\xd2\x04\xe8\x04\xbe\x9f\x1a/\xdd\\}@\xe8\x04\x04\x05q=\n'
b'\xd7\xa3p=\xd0v@\xdd\x00\xea\x00\xd1"\xdb\xf9~Jy@\xea\x003\x01J\x0c\x02+\x87|y@3\x01D\x01J\x0c\x02+\x87\xbcz@D\x01S\x01\x17\xd9\xce\xf7S\xffz@S\x01Y\x01\x17\xd9\xce\xf7S={@Y\x01z\x01o\x12\x83\xc0\xcaY{@z\x01\x8a\x01\x96C\x8bl\xe7\xe5{@\x8a\x01\x8d\x015^\xbaI\x0c$|@\x8d\x01\x95\x01\xe5\xd0"\xdb\xf9B|@\x95\x01\xbc\x01\x1f\x85\xebQ\xb8\xe6|@\xbc\x01\xce\x01d;\xdfO\x8dQ}@\xce\x01\xd7\x01\xbe\x9f\x1a/\xdd\\}@\xd7\x01\xe6\x01q=\n
b'\xd7\xa3p\xef\x86@\x00\x03#\x03j\xbct\x93\x18s\x87@#\x03&\x039\xb4\xc8v\xbe\xb3\x87@&\x03(\x03\n'

b'\xd7\xa3\x96}@\x04\x05\r\x05\xd5x\xe9&1\xf4}@\r\x050\x05\xac\x1cZd;\x0b~@0\x050\x05#\xdb\xf9~j\xba~@333333,@'
b'\xd7\xa3\x96}@\xe6\x01\x0c\x02\xac\x1cZd;\x0b~@\x0c\x02.\x02%\x06\x81\x95CE~@.\x02\xfb\x02\xd7\xa3p=\n
b'\xd7\xa3p=\xc1\x87@(\x038\x03\x0c\x02+\x87\x16\xdc\x87@8\x03I\x03\x0c\x02+\x87\x16\xe4\x87@I\x03T\x03\x12\x83\xc0\xca\xa1\x0b\x88@T\x03\xe2\x03^\xbaI\x0c\x02#\x88@\xe2\x03\xcf\x04h\x91\xed|?\\\x89@\xcf\x04\xcf\x04L7\x89A`\x82\x8c@333333,@''
'''