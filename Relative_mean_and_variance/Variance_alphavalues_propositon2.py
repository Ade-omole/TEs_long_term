import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

####### Parameters
beta1 = 0.3
beta2 = 0.4
delta1 = 0.2
delta2 = 0.1
V = 500
alphas = [0.35, 0.40, 0.45, 0.5, 0.55, 0.6]
alphas_fine = np.arange(0.35, 0.61, 0.01)  # Finer alpha range from 0.35 to 0.6 with interval of 0.01

####### Lists to store results for original alphas
data_X_variance_list = []
data_Y_variance_list = []
estimated_X_variance_list = []
estimated_Y_variance_list = []
analytic_X_variance_list = []
analytic_Y_variance_list = []
analytic_covariance_list = []
data_covariance_list = []

####### Iterate through each alpha value in original alphas
for alpha in alphas:
    ####### Update file path
    file_path = f"/path/variance_alpha_values_{alpha}.csv"

    ####### Read data
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

    ####### Variances from the data
    cov_matrix = np.cov(data_X_V, data_Y_V)
    data_X_variance = cov_matrix[0][0] / V ** 2
    data_Y_variance = cov_matrix[1][1] / V ** 2

    ####### Eq. 16 estimated variances
    covariance = cov_matrix[0][1]
    estimated_X_variance = (mean_X * (1 + (covariance / mean_Y) + (mean_X / mean_Y))) / V ** 2
    estimated_Y_variance = (covariance * (1 + (covariance / mean_X))) / V ** 2

    ####### Analytic variances
    analytic_X_variance = (beta1 * (2 * alpha * beta1 - delta1 * delta2) / (
                beta2 * (alpha * beta1 - delta1 * (beta1 + delta2)))) / V
    analytic_Y_variance = (alpha * beta1 * (alpha + delta1) / (beta2 * delta1 ** 2)) / V

    ####### Analytic covariance
    analytic_covariance = (alpha * beta1 / (beta2 * delta1)) / V

    ####### Estimated covariance from data
    data_covariance = cov_matrix[0][1] / V ** 2

    ####### Append results to the lists
    data_X_variance_list.append(data_X_variance)
    data_Y_variance_list.append(data_Y_variance)
    estimated_X_variance_list.append(estimated_X_variance)
    estimated_Y_variance_list.append(estimated_Y_variance)
    analytic_X_variance_list.append(analytic_X_variance)
    analytic_Y_variance_list.append(analytic_Y_variance)
    analytic_covariance_list.append(analytic_covariance)
    data_covariance_list.append(data_covariance)

####### Lists to store results for finer alphas
analytic_X_variance_list_fine = []

####### Compute only analytic X variance for finer alphas
for alpha in alphas_fine:
    analytic_X_variance = (beta1 * (2 * alpha * beta1 - delta1 * delta2) / (
                beta2 * (alpha * beta1 - delta1 * (beta1 + delta2)))) / V
    analytic_X_variance_list_fine.append(analytic_X_variance)

####### Plotting
plt.figure(figsize=(12, 7))

####### Plot X variances
plt.plot(alphas_fine, analytic_X_variance_list_fine, 'b-', label='Autonomous Analytic Variance', linewidth=3)
plt.plot(alphas, estimated_X_variance_list, 'b^', markerfacecolor='none', label='Autonomous Estimated Variance', markersize=10, linewidth=2)
plt.plot(alphas, data_X_variance_list, 'bo', label='Autonomous Data Variance', markersize=10, linewidth=2)

####### Plot Y variances
plt.plot(alphas, analytic_Y_variance_list, 'g-', label='Nonautonomous Analytic Variance', linewidth=3)
plt.plot(alphas, estimated_Y_variance_list, 'g^', markerfacecolor='none', label='Nonautonomous Estimated Variance', markersize=10, linewidth=2)
plt.plot(alphas, data_Y_variance_list, 'go', label='Nonautonomous Data Variance', markersize=10, linewidth=2)

####### Draw a vertical line at alpha = 0.47
plt.axvline(x=0.46667, color='k', linestyle='--', linewidth=2)

####### Add shaded regions
plt.axvspan(0.35, 0.3926, color='gray', alpha=0.2)
plt.axvspan(0.3926, 0.46667, color='yellow', alpha=0.2)

####### Add labels and legend
plt.xlabel(r'$\alpha$ (number of complex)', fontsize=16, fontweight='bold')
plt.ylabel('Variance', fontsize=14, fontweight='bold')
plt.legend(fontsize=12, loc='upper left')
# plt.legend(prop={'weight': 'bold'})
plt.grid(True)
plt.show()
