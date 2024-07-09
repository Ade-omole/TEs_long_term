import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

####### File path to the CSV data
file_path = "/Users/adekanmiomole/Julia_programming/gillespie_data_10000_distribution.csv"

####### Parameters
alpha = 2.0
beta1 = 1.0
beta2 = 0.5
delta1 = 0.5
delta2 = 1.5
V = 500

####### Read the data from CSV file
df = pd.read_csv(file_path)

####### Extract data from the DataFrame
data_X = df['x']
data_Y = df['y']

####### Scale the data by V
data_X_V = [x * V for x in data_X]
data_Y_V = [x * V for x in data_Y]

####### Calculate the mean of the scaled data
mean_X = np.mean(data_X_V)
mean_Y = np.mean(data_Y_V)

####### Analytic mean calculations
analytic_X_mean = beta1 / beta2
analytic_Y_mean = (beta1 * alpha - delta1 * delta2 - beta1 * delta1) / (beta2 * delta1)

####### Estimated covariance from the data
cov_matrix = np.cov(data_X_V, data_Y_V)
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
mu_X_data, sigma_X_data = mean_X, np.std(data_X)
mu_Y_data, sigma_Y_data = mean_Y, np.std(data_Y)

# Eq. 16 estimation
mu_X_v, sigma_X_est_v = mean_X, np.sqrt(Eq16_es_X_variance)
mu_Y_v, sigma_Y_est_v = mean_Y, np.sqrt(Eq16_es_Y_variance)

# Eq. 16 analytic
mu_X_an_v, sigma_X_an_v = analytic_X_mean, np.sqrt(Eq16_an_X_variance)
mu_Y_an_v, sigma_Y_an_v = analytic_Y_mean, np.sqrt(Eq16_an_Y_variance)

# Eq. 33 analytic
mu_X_an_k, sigma_X_an_k = analytic_X_mean, np.sqrt(Eq33_an_X_variance)
mu_Y_an_k, sigma_Y_an_k = analytic_Y_mean, np.sqrt(Eq33_an_Y_variance)

####### Print variance calculations
print("Eq. 16 X Variance:", Eq16_an_X_variance)
print("Eq. 16 Y Variance:", Eq16_an_Y_variance)
print("\nEq. 33 X Variance:", Eq33_an_X_variance)
print("Eq. 33 Y Variance:", Eq33_an_Y_variance)
print("\nEstimated covariance", covariance / V**2)
print("Analytic covariance:", analytic_covariance)

####### Plot histogram and normal distribution for X
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.hist(data_X, bins=20, density=True, alpha=0.6, color='b', label='Histogram')
x = np.linspace(min(data_X), max(data_X), 100)
plt.plot(x, norm.pdf(x, mu_X_v/V, sigma_X_est_v), color='r', linewidth=3, label='Eq. 16 estimated distribution')
plt.plot(x, norm.pdf(x, mu_X_an_k, sigma_X_an_k), color='orange', linewidth=3, label='Eq. 33 distribution')
plt.plot(x, norm.pdf(x, mu_X_an_v, sigma_X_an_v), color='k', linewidth=3, linestyle=':', label='Eq. 16 distribution')
plt.xlabel('Copy number (scaled)', fontsize=16, fontweight='bold')
plt.ylabel('Density', fontsize=16, fontweight='bold')
plt.title('Autonomous TE (x)', fontsize=16, fontweight='bold')
plt.legend()
plt.legend(prop={'weight': 'bold'}) # make the legend text thicker

####### Plot histogram and normal distribution for Y
plt.subplot(1, 2, 2)
plt.hist(data_Y, bins=20, density=True, alpha=0.6, color='g', label='Histogram')
y = np.linspace(min(data_Y), max(data_Y), 100)
plt.plot(y, norm.pdf(y, mu_Y_v/V, sigma_Y_est_v), color='r', linewidth=3, label='Eq. 16 estimated distributionn')
plt.plot(y, norm.pdf(y, mu_Y_an_k, sigma_Y_an_k), color='orange', linewidth=3, label='Eq. 33 distribution')
plt.plot(y, norm.pdf(y, mu_Y_an_v, sigma_Y_an_v), color='k', linewidth=3, linestyle=':', label='Eq. 16 distribution')
plt.xlabel('Copy number (scaled)', fontsize=16, fontweight='bold')
plt.ylabel('Density', fontsize=16, fontweight='bold')
plt.title('Nonautonomous TE (y)', fontsize=16, fontweight='bold')
plt.legend()
plt.legend(prop={'weight': 'bold'})

####### Adjust layout and show the plot
plt.tight_layout()
plt.show()
