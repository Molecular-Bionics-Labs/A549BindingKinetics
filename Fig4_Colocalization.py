import matplotlib.pyplot as plt
import numpy as np

# Original data
ves = np.array([0.015, 0.12, 0.20, 0.31, 0.44]) 
mic = np.array([0, 0.202, 0.284, 0.354, 0.412]) 
tube = np.array([0.0406, 0.104, 0.160, 0.228, 0.310])

# EEA1 data
ves_eea1 = np.array([0.0116, 0.0378, 0.0438, 0.0524, 0.110])
mic_eea1 = np.array([0.00101, 0.00604, 0.0101, 0.0186, 0.0368])
tube_eea1 = np.array([0, 0.0307, 0.0504, 0.0806, 0.00504, 0.147])

# Combine data for boxplot (main data and EEA1 colocalization data)
data = [ves, ves_eea1, mic, mic_eea1, tube, tube_eea1]
colors = ['black', 'gray', 'black', 'gray', 'black', 'gray']  # Alternating colors

# Create figure and axes
fig, ax = plt.subplots(figsize=(8, 6))  # Set figure size

# Loop through each dataset and plot the boxplot
for i, (dataset, color) in enumerate(zip(data, colors), start=1):
    ax.boxplot(dataset, 
               positions=[i],       # Set x-position for each box
               vert=True,           # Vertical boxplots
               patch_artist=True,   # Allows filling the boxes with color
               widths=0.5,          # **Make the boxes wider**
            #    showfliers=False,     # **Don't show outliers (ignore them)**
               
               # Box properties (color and border thickness)
               boxprops=dict(facecolor='white', color=color, linewidth=3),
               
               # Whisker line properties (extend from box)
               whiskerprops=dict(color=color, linewidth=3),
               
               # Cap properties (horizontal lines at the ends of whiskers)
               capprops=dict(color=color, linewidth=3),
               
               # Median line properties (middle line in the box)
               medianprops=dict(color=color, linewidth=3))

# Set x-tick positions and labels
ax.set_xticks(range(1, 7))
ax.set_xticklabels(['ves', 'ves_eea1', 'mic', 'mic_eea1', 'tube', 'tube_eea1'])

# Set labels and title
ax.set_ylabel("Normalized Values")  # Label for the y-axis
ax.set_title("Box Plot Ignoring Outliers")  # Chart title

# Save and show plot
plt.savefig("/Users/ttomis9651/Downloads/colocalization_without_outliers.svg", dpi=800)  # Save as high-quality SVG
plt.show()  # Display the plot
