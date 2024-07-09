import pandas as pd
import matplotlib.pyplot as plt

####### File path to the CSV data from Relative_silencing_error.py
file_path = "/path/relative_silencing_error.csv"

####### Read the data from CSV file
df = pd.read_csv(file_path)

####### Extract data from the DataFrame
x_values = df['silencing rates']
y1_values = df['relative error X']
y2_values = df['relative error Y']

####### Plot the data using a line graph
plt.plot(x_values, y1_values, label='Autonomous (x)', color='b', marker='o', linewidth=3)
plt.plot(x_values, y2_values, label='Nonautonomous (y)', color='g', marker='o', linewidth=3)

####### Customize the plot
plt.xlabel('Silencing rates', fontsize=14, fontweight='bold')  ####### Set x-axis label
plt.ylabel('Relative Errors', fontsize=14, fontweight='bold')  ####### Set y-axis label
plt.xscale('log')  ####### Set x-axis to logarithmic scale

####### Set x-axis ticks, including intermediate points
ticks = [2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4]
tick_labels = [r'$2 \times 10^{-6}$', r'$5 \times 10^{-6}$', r'$1 \times 10^{-5}$', r'$2 \times 10^{-5}$', r'$5 \times 10^{-5}$', r'$1 \times 10^{-4}$']
plt.xticks(ticks, tick_labels)

####### Remove minor ticks
plt.minorticks_off()

####### Add legend
plt.legend()
plt.legend(prop={'weight': 'bold'})

####### Rotate x-axis labels by 0 degrees (default)
plt.xticks(rotation=0)

####### Enable grid
plt.grid(True)

####### Show the plot
plt.show()
