import numpy as np
import pandas as pd
import os

####### Function to calculate average data variance and plot the bar chart
def calculate_and_plot(file_paths, label_suffix):
    print(f"Starting calculate_and_plot with label_suffix: {label_suffix}")
    ####### Initialize variables to store cumulative values
    total_data_variance_X = 0
    total_data_variance_Y = 0

    ####### Lists to store individual data variance values for plotting
    data_variance_X_values = []
    data_variance_Y_values = []

    ####### Iterate through each file path
    for idx, file_path in enumerate(file_paths, start=1):
        print(f"Reading file {file_path}")
        if not os.path.exists(file_path):
            print(f"File {file_path} does not exist.")
            continue

        ####### Read the CSV file
        df = pd.read_csv(file_path)

        ####### Concatenate X and Y values
        all_X = pd.concat([df['X_1'], df['X_2']], ignore_index=True)
        all_Y = pd.concat([df['Y_1'], df['Y_2']], ignore_index=True)

        ####### Normalize the data
        data_X = np.divide(all_X, 500)
        data_Y = np.divide(all_Y, 500)

        ####### Compute the covariance matrix
        cov_matrix = np.cov(data_X, data_Y)

        ####### Accumulate the values
        total_data_variance_X += cov_matrix[0][0]
        total_data_variance_Y += cov_matrix[1][1]

        ####### Append individual data variance values for plotting
        data_variance_X_values.append(cov_matrix[0][0])
        data_variance_Y_values.append(cov_matrix[1][1])

    ####### Calculate the average values
    average_data_variance_X = total_data_variance_X / len(file_paths)
    average_data_variance_Y = total_data_variance_Y / len(file_paths)

    print(f"Average Data Variance X{label_suffix}: {average_data_variance_X}")
    print(f"Average Data Variance Y{label_suffix}: {average_data_variance_Y}")

    ####### Calculate relative error between average variances of te_0.0 and other file paths
    return calculate_relative_error(average_data_variance_X, average_data_variance_Y, file_paths[1:], label_suffix)

####### Function to calculate relative error between average variances of reference and other file paths
def calculate_relative_error(reference_variance_X, reference_variance_Y, other_file_paths, label_suffix):
    print(f"Starting calculate_relative_error with label_suffix: {label_suffix}")
    relative_error_results = []

    ####### Iterate through sets of other file paths
    for other_label_suffix, other_file_paths in zip(["1.0e-6", "2.0e-6", "3.0e-6", "4.0e-6",
                                                     "5.0e-6", "6.0e-6", "7.0e-6", "8.0e-6",
                                                     "9.0e-6", "1.0e-5", "2.0e-5", "3.0e-5",
                                                     "4.0e-5", "5.0e-5", "6.0e-5", "7.0e-5",
                                                     "8.0e-5", "9.0e-5"],
                                                    [file_paths_set2, file_paths_set3, file_paths_set4, file_paths_set5,
                                                     file_paths_set6, file_paths_set7, file_paths_set8, file_paths_set9,
                                                     file_paths_set10, file_paths_set11, file_paths_set12,
                                                     file_paths_set13, file_paths_set14, file_paths_set15,
                                                     file_paths_set16, file_paths_set17, file_paths_set18,
                                                     file_paths_set19]):
        print(f"Processing label_suffix: {other_label_suffix}")
        total_relative_error_X = 0
        total_relative_error_Y = 0

        ####### Iterate through each file in the other file paths
        for idx, other_file in enumerate(other_file_paths, start=1):
            print(f"Reading file {other_file}")
            if not os.path.exists(other_file):
                print(f"File {other_file} does not exist.")
                continue

            other_df = pd.read_csv(other_file)

            ####### Concatenate X and Y values
            other_X = pd.concat([other_df['X_1'], other_df['X_2']], ignore_index=True)
            other_Y = pd.concat([other_df['Y_1'], other_df['Y_2']], ignore_index=True)

            ####### Normalize the data
            data_X = np.divide(other_X, 500)
            data_Y = np.divide(other_Y, 500)

            ####### Compute the covariance matrix
            cov_matrix = np.cov(data_X, data_Y)

            ####### Calculate the average variances for the other files
            average_variance_X = np.mean(np.var(data_X))
            average_variance_Y = np.mean(np.var(data_Y))

            ####### Calculate and accumulate the relative error for X and Y
            total_relative_error_X += np.abs((average_variance_X - reference_variance_X) / reference_variance_X)
            total_relative_error_Y += np.abs((average_variance_Y - reference_variance_Y) / reference_variance_Y)

        ####### Calculate the average relative error for X and Y
        average_relative_error_X = total_relative_error_X / len(other_file_paths)
        average_relative_error_Y = total_relative_error_Y / len(other_file_paths)

        print(f"Average Relative Error X{label_suffix} vs {other_label_suffix}: {average_relative_error_X}")
        print(f"Average Relative Error Y{label_suffix} vs {other_label_suffix}: {average_relative_error_Y}")

        ####### Append the results
        relative_error_results.append({
            "silencing rates": other_label_suffix,
            "relative error X": average_relative_error_X,
            "relative error Y": average_relative_error_Y
        })

    ####### Return the results as a DataFrame
    return pd.DataFrame(relative_error_results)

####### Generate file paths for the first set of simulations (te_0.0)
base_path = "/path/"
file_paths_set1 = [f"{base_path}te_0.0_{i}_data.csv" for i in range(1, 101)]
file_paths_set2 = [f"{base_path}te_1.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set3 = [f"{base_path}te_2.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set4 = [f"{base_path}te_3.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set5 = [f"{base_path}te_4.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set6 = [f"{base_path}te_5.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set7 = [f"{base_path}te_6.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set8 = [f"{base_path}te_7.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set9 = [f"{base_path}te_8.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set10 = [f"{base_path}te_9.0e-6_{i}_data.csv" for i in range(1, 101)]
file_paths_set11 = [f"{base_path}te_1.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set12 = [f"{base_path}te_2.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set13 = [f"{base_path}te_3.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set14 = [f"{base_path}te_4.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set15 = [f"{base_path}te_5.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set16 = [f"{base_path}te_6.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set17 = [f"{base_path}te_7.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set18 = [f"{base_path}te_8.0e-5_{i}_data.csv" for i in range(1, 101)]
file_paths_set19 = [f"{base_path}te_9.0e-5_{i}_data.csv" for i in range(1, 101)]

####### Call the function for the first set of file paths
print("Calling calculate_and_plot")
result_df = calculate_and_plot(file_paths_set1, " (te_0.0)")

####### Check the contents of result_df
print("Contents of result_df:")
print(result_df)

####### Check if the directory exists and is writable
print("Checking directory existence and permissions")
directory = "/path/"
file_path = os.path.join(directory, "relative_silencing_error.csv")

if not os.path.exists(directory):
    print(f"Directory {directory} does not exist.")
elif not os.access(directory, os.W_OK):
    print(f"Directory {directory} is not writable.")
else:
    try:
        ####### Attempt to write the DataFrame to CSV
        result_df.to_csv(file_path, index=False)
        print(f"File successfully written to {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

####### Write the result DataFrame to CSV
result_df.to_csv("/path/relative_silencing_error.csv", index=False)
