import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Define a list of colors for the custom colormap (smooth transition between colors)
list_color = [[16/256, 70/256, 128/256], [49/256, 124/256, 183/256], 
                 [109/256, 173/256, 209/256], [182/256, 215/256, 232/256],
                 [233/256, 241/256, 244/256]]

# Create the custom colormap
custom_cmap = LinearSegmentedColormap.from_list("smooth_cmap", list_color, N=100)

# Generate some data for demonstration (a 2D array)
data = np.random.rand(10, 10)  # 10x10 random data

# Create the plot using the custom colormap
plt.figure(figsize=(6, 5))
plt.imshow(data, cmap=custom_cmap, interpolation='nearest')

# Add a colorbar to show the mapping of values to colors
plt.colorbar()

# Set a title
plt.title("Custom Smooth Colormap", fontsize=16)

# Show the plot
plt.show()