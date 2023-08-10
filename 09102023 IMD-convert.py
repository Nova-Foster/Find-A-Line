import numpy as np
import matplotlib.pyplot as py

with open("800V - 250Hz - 500usmm - 15um vertical.imd", "rb") as file:      #Open the file reading as non-text
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

#Plot the data
py.imshow(data_2d_array,interpolation="nearest",origin="lower")
py.colorbar()
py.show()



#Horizontal profile
hori = np.sum(data_2d_array,axis=0)
py.plot(np.arange(0,len(hori),1),hori)
py.show()