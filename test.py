from matplotlib import colors
import numpy as np
import matplotlib.pyplot as py

x_coords = wavelengths
y_coords = times[start_y:stop_y]
intensity = intensity[start_y:stop_y,:]
z_vals = np.where(intensity<5,float("Nan"),intensity)

threed_lines_plot = py.figure()
ax = py.axes(projection="3d")

plot_colors = py.get_cmap("rainbow",(stop_y-start_y)*10)


for i in range(len(y_coords)):
  current_time = np.linspace(y_coords[i],y_coords[i],1392)
  current_intensity = z_vals[i]
  current_color = plot_colors(i)
  ax.plot(x_coords,current_time,current_intensity,color=current_color)

ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Time (μs)")
ax.set_zlabel("Intensity")
py.show()