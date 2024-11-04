import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.font_manager import FontProperties

######## Define font properties with larger, bold font for the legend
legend_font = FontProperties(weight='bold', size=13)

######## Read CSV file
### No recombination
file_path = "/paths/pop_distribution_0_recombination.csv"

### With recombination
# file_path = "/paths/pop_distribution_0.01_recombination.csv"
# file_path = "/paths/pop_distribution_0.1_recombination.csv"
# file_path = "/paths/pop_distribution_1_recombination.csv"

df = pd.read_csv(file_path)

####### Parameters
alpha = 2.0
beta1 = 1.0
beta2 = 0.5
delta1 = 0.5
delta2 = 1.5
V = 500

# Concatenate columns X_1, X_2 and Y_1, Y_2
concatenated_X = pd.concat([df['X_1'], df['X_2']])
concatenated_Y = pd.concat([df['Y_1'], df['Y_2']])

####### Scale the data by V
data_X_V = [x / V for x in concatenated_X]
data_Y_V = [x / V for x in concatenated_Y]


####### Calculate the mean of the scaled data
mean_X = np.mean(concatenated_X)
mean_Y = np.mean(concatenated_Y)

####### Analytic mean calculations
analytic_X_mean = beta1 / beta2
analytic_Y_mean = (beta1 * alpha - delta1 * delta2 - beta1 * delta1) / (beta2 * delta1)

####### Estimated covariance from the data
cov_matrix = np.cov(concatenated_X, concatenated_Y)
covariance = cov_matrix[0][1]

####### Analytic covariance calculation
analytic_covariance = (alpha * beta1 / (beta2 * delta1)) / V

####### Eq. 16 variance estimation
Eq16_es_X_variance = (mean_X * (1 + (covariance / mean_Y) + (mean_X / mean_Y))) / V**2
Eq16_es_Y_variance = (covariance * (1 + (covariance / mean_X))) / V**2

####### Eq. 16 analytic variance calculation
Eq16_an_X_variance = (beta1 * (2 * alpha * beta1 - delta1 * delta2) / (
            beta2 * (alpha * beta1 - delta1 * (beta1 + delta2)))) / V
Eq16_an_Y_variance = (alpha * beta1 * (alpha + delta1) / (beta2 * delta1 ** 2)) / V

####### Eq. 33 analytic variance calculation
Eq33_an_X_variance = (beta1 / beta2 + ((2 * delta1**2 * delta2 + beta1**2) / alpha**2 * beta2) + beta1**2 * (alpha + delta1) / (
            beta2 * (alpha * beta1 - delta1 * (beta1 + delta2)))) / V
Eq33_an_Y_variance = (alpha * beta1 * (alpha + delta1) / (beta2 * delta1 ** 2)) / V

####### Calculate mean and standard deviation for various estimations
# Data
mu_X_data, sigma_X_data = mean_X, np.std(data_X_V)
mu_Y_data, sigma_Y_data = mean_Y, np.std(data_Y_V)

# Eq. 16 estimation
mu_X_v, sigma_X_est_v = mean_X, np.sqrt(Eq16_es_X_variance)
mu_Y_v, sigma_Y_est_v = mean_Y, np.sqrt(Eq16_es_Y_variance)

# Eq. 16 analytic
mu_X_an_v, sigma_X_an_v = analytic_X_mean, np.sqrt(Eq16_an_X_variance)
mu_Y_an_v, sigma_Y_an_v = analytic_Y_mean, np.sqrt(Eq16_an_Y_variance)

# Eq. 33 analytic
mu_X_an_k, sigma_X_an_k = analytic_X_mean, np.sqrt(Eq33_an_X_variance)
mu_Y_an_k, sigma_Y_an_k = analytic_Y_mean, np.sqrt(Eq33_an_Y_variance)

####### Print variance and mean values
print("Eq. 16 X Variance:", Eq16_an_X_variance)
print("Eq. 16 Y Variance:", Eq16_an_Y_variance)
print("\nEq. 33 X Variance:", Eq33_an_X_variance)
print("Eq. 33 Y Variance:", Eq33_an_Y_variance)
print("\ndata X variance:", sigma_X_data ** 2)
print("data Y variance:", sigma_Y_data ** 2)
print("\nEq. 16 estimated X variance:", Eq16_es_X_variance)
print("Eq. 16 estimated Y variance:", Eq16_es_Y_variance)
print("\nEstimated covariance", covariance / V**2)
print("Analytic covariance:", analytic_covariance)
print("\ndata mean X:", mean_X / V)
print("data mean Y:", mean_Y / V)


####### Plot histogram and normal distribution for X
plt.figure(figsize=(15, 6))
plt.subplot(1, 2, 1)
# plt.hist(data_X, bins=20, density=True, alpha=0.6, color='b', label='Histogram')
plt.hist(data_X_V, bins=20, density=True, alpha=0.6, color='b')

x = np.linspace(min(data_X_V), max(data_X_V), 100)
plt.plot(x, norm.pdf(x, mu_X_v/V, sigma_X_est_v), color='r', linewidth=3, label='Eq. (16) estimated')
plt.plot(x, norm.pdf(x, mu_X_an_k, sigma_X_an_k), color='orange', linewidth=3, label='Eq. (33)')
plt.plot(x, norm.pdf(x, mu_X_an_v, sigma_X_an_v), color='k', linewidth=3, linestyle=':', label='Eq. (16)')
plt.xlabel('Copy number (scaled)', fontsize=20, fontweight='bold')
plt.ylabel('Density', fontsize=20, fontweight='bold')
plt.title('Autonomous TE (x)', fontsize=20, fontweight='bold')
plt.legend(prop=legend_font)  # Apply font properties


# Set font size and make tick labels bold
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')

####### Plot histogram and normal distribution for Y
plt.subplot(1, 2, 2)
# plt.hist(data_Y, bins=20, density=True, alpha=0.6, color='g', label='Histogram')
plt.hist(data_Y_V, bins=20, density=True, alpha=0.6, color='g')
y = np.linspace(min(data_Y_V), max(data_Y_V), 100)
plt.plot(y, norm.pdf(y, mu_Y_v/V, sigma_Y_est_v), color='r', linewidth=3, label='Eq. (16) estimated')
plt.plot(y, norm.pdf(y, mu_Y_an_k, sigma_Y_an_k), color='orange', linewidth=3, label='Eq. (33)')
plt.plot(y, norm.pdf(y, mu_Y_an_v, sigma_Y_an_v), color='k', linewidth=3, linestyle=':', label='Eq. (16)')
plt.xlabel('Copy number (scaled)', fontsize=20, fontweight='bold')
plt.ylabel('Density', fontsize=20, fontweight='bold')
plt.title('Nonautonomous TE (y)', fontsize=20, fontweight='bold')
plt.legend(prop=legend_font)  # Apply font properties

# Set font size and make tick labels bold
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')

####### Adjust layout and show the plot
plt.tight_layout()
plt.show()

