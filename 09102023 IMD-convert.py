import numpy as np
import codecs


stream_reader = codecs.open("20230810 - 500usmm - 800V - 250Hz -00000.IMD","rb")     #Create stream reader and open non text file
raw = stream_reader.readlines()                                                    #Read all of the data and assign to raw
stream_reader.close()                                                              #Close stream reader

raw = str(raw)                                                                     #Convert data into a string so it can be split
raw_split = raw.split("\\")

#Header is first 6 values
header = raw_split[:6]
# 0 = '[b"'
# 1 = '01'
# 2 = '01p'
# 3 = '05'
# 4 = '10'
# 5 = '04'

#Data is the rest
data = raw_split[6:]

#Need to concatenate into pairs: i.e. data[1] & data[0] in a single value
#That joined data 




#Start of data
# 6 = '00'
# 7 = ....