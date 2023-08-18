'''
Created by Nova Foster (They/Them) 16/08/2023
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

freq = 1000
time_per = (1/freq)*1000
time_base = 0.2*0.7
#time_per_pixel = time_base/1040
ddg_per = 3.000125
no_sparks = 4
fraction = 0.025
times_to_check = 40

no_sparks = int( np.ceil(time_base+(fraction*times_to_check)) ) + 1   #No. sparks that will be interacted with in total based + 1 so the plot doesn't change axis
#sparks_between = ddg_per / time_per

sparks = np.arange(0,no_sparks*time_per,time_per)



def frame(interval):
    fig.clear()
    ys = np.linspace(1,1,no_sparks)

    plt.plot(sparks,ys,".",label="Sparks")

    offset = interval * fraction
    window = np.linspace(offset,time_base+offset,1040)
    window_y = np.linspace(1,1,1040)
    plt.plot(window,window_y,"r-")
    plt.xlabel("Time(ms)")

    max_window = max(window)
    min_window = min(window)
    hit[interval] =  min_window < sparks.any() < max_window




hit = np.zeros(times_to_check)

for i in range(times_to_check):
    offset = i * fraction
    window = np.linspace(offset,time_base+offset,1040)
    max_window = max(window)
    min_window = min(window)

    for j in range(no_sparks):
        hit[i] =  (min_window < sparks[j] < max_window) + hit[i]

print(hit)


fig, ax = plt.subplots()
ani = anim.FuncAnimation(fig=fig, func=frame, frames=times_to_check, interval=300)
plt.show()