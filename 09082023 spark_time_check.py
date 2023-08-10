'''
Created by Nova Foster (They/Them) 09/08/2023
'''
import numpy as np
import scipy as sp
Timebase = 0
Frequency = 0
pixels = 1040

print("Please input the frequency (Hz): ")
Frequency = float(input())

print("Please input the timebase (ms): ")
Timebase = float(input())
spark_gap = (1/Frequency)*1000                               #Time between sparks
pixel_in_time = Timebase/pixels                          #How much time one pixel is worth (ms pixel^-1)
spark_gap_inpixels = int(np.ceil(spark_gap/pixel_in_time))    #How many pixels between sparks, rounded up

#loop for each starting pixel & see how many sparks are seen for each


'''
Need to have half a spark time, timebase
The number of pixels between sparks is the number of combinations
'''
region_end = int(np.ceil(pixels+(0.5*spark_gap_inpixels)))        #Has to be rounded otherwise there is a chance of indexing by flaot, should only be 1 more spark at worse
#Garunteed to have a pixel on pixel i (imitating when the ddg is started)
#Work out where the corresponding spark pixels are by adding and subbing spark_gap_pixels until it throws an error
#save how many spark pixels there are in that range
sparks = np.zeros([spark_gap_inpixels,region_end])
no_sparks = region_end//spark_gap_inpixels
for i in range(spark_gap_inpixels):


    count = i
    for j in range(no_sparks):
        sparks[i][count] = 1
        count = count+spark_gap_inpixels
    ''' Every spark gap in pixels, have a 1 starting from i
    use booleans with mod?, int(x==true) should convert to 1 and 0'''

    '''
    copy last 1040 values from with_sparks
    add them all up & divide by 1040, floor it so that it is always rounded down
    '''
useful_region = np.arange(int(np.ceil(0.5*spark_gap_inpixels)),region_end,1)
spark_in_region = np.take(sparks,useful_region,axis=1)
spark_sum = np.sum(spark_in_region,axis=1)

nos,counts = np.unique(spark_sum,return_counts=True)
total = dict(zip(nos,counts))


print(total)
