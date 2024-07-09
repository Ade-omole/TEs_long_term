using Pkg
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("CSV")
using Random
using Distributions
using CSV

####### Function to simulate the dynamics of transposable elements (TEs) in a population
function simulate_TE_dynamics(population_size, generations, alpha, beta1, beta2, delta1, delta2, dt, V, silencing_rate)

    ####### Initialize TE copy numbers for the initial generation
    te_copy_numbers1 = [[200, 200, false] for _ in 1:population_size]
    te_copy_numbers2 = [[200, 200, false] for _ in 1:population_size]

    for i in 1:100

        ####### Simulate the Wright-Fisher process with TE replication and silencing
        for generation in 1:generations

            ####### Create a new generation of TE copy numbers
            new_te_copy_numbers1 = []
            new_te_copy_numbers2 = []

            for individual in 1:population_size

                size = collect(1:population_size)

                index_choice1 = size[rand(1:length(size))]

                remaining_size = [index for index in size if index != index_choice1]

                index_choice2 = remaining_size[rand(1:length(remaining_size))]

                ################ TEs replicate before silencing: copy and paste transposition mechanism

                ####### Randomly select two parent TEs (with replacement) for TE1 with their corresponding tuples
                Ind_choice1 = rand([1,2])
                if Ind_choice1 == 1
                    random_te1 = te_copy_numbers1[index_choice1]
                    corresponding_tuple1 = te_copy_numbers2[index_choice1]
                else
                    random_te1 = te_copy_numbers2[index_choice1]
                    corresponding_tuple1 = te_copy_numbers1[index_choice1]
                end

                ####### Randomly select two parent TEs (with replacement) for TE2 with their corresponding tuples
                Ind_choice2 = rand([1,2])
                if Ind_choice2 == 2
                    random_te2 = te_copy_numbers2[index_choice2]
                    corresponding_tuple2 = te_copy_numbers1[index_choice2]
                else
                    random_te2 = te_copy_numbers1[index_choice2]
                    corresponding_tuple2 = te_copy_numbers2[index_choice2]
                end

                ####### For TE replication per generation for te_copy_numbers1 and te_copy_numbers2
                av_complex = [alpha * random_te1[1]/(beta1 + delta2 + beta2 * random_te1[2]), alpha * random_te2[1]/(beta1 + delta2 + beta2 * random_te2[2])]

                propensity2 = [beta1 * av_complex[1] * dt,  beta1 * av_complex[2] * dt]
                poisson2 = [rand(Poisson(propensity2[1])), rand(Poisson(propensity2[2]))]
                tau2 = [poisson2[1], poisson2[2]]

                propensity3 = [beta2 * av_complex[1] * random_te1[2] * dt, beta2 * av_complex[2] * random_te2[2] * dt]
                poisson3 = [rand(Poisson(propensity3[1])), rand(Poisson(propensity3[2]))]
                tau3 = [poisson3[1], poisson3[2]]

                propensity4 = [delta1 * random_te1[1] * dt, delta1 * random_te2[1] * dt]
                poisson4 = [rand(Poisson(propensity4[1])), rand(Poisson(propensity4[2]))]
                tau4 = [poisson4[1], poisson4[2]]

                propensity5 = [delta1 * random_te1[2] * dt, delta1 * random_te2[2] * dt]
                poisson5 = [rand(Poisson(propensity5[1])), rand(Poisson(propensity5[2]))]
                tau5 = [poisson5[1], poisson5[2]]

                ################################################## Simulate TE replication and silencing

                ############################## First TE tuple
                silencing_status1 = false

                ####### Replication and/or silencing of TE1
                first_replicated1 = false 

                if random_te2[3] == true
                    local first_random_te1 = random_te1[1]
                    silencing_status1 = true
                else
                    ####### use the poisson process for replication
                    first_random_te_rep1 = random_te1[1] + tau2[1] - tau4[1]
                    if true
                        first_replicated1 = true
                        if first_replicated1 && rand() < silencing_rate 
                            first_random_te1 = max(first_random_te_rep1 + 0.0, 0)
                            silencing_status1 = true
                        else
                            first_random_te1 = max(first_random_te_rep1, 0)
                        end
                    else
                        first_random_te1 = random_te1[1]
                    end
                end

                ####### Replication and/or silencing of TE2
                second_replicated1 = false 
                if silencing_status1 == true
                    local second_random_te1 = random_te1[2]
                    silencing_status1 = true
                else
                    ####### use the poisson process for replication
                    second_random_te_rep1 = random_te1[2] + tau3[1] - tau5[1]
                    if true
                        second_replicated1 = true
                        if second_replicated1 && rand() < silencing_rate 
                            second_random_te1 = max(second_random_te_rep1 + 0.0, 0)
                            silencing_status1 = true
                        else
                            second_random_te1 = max(second_random_te_rep1, 0)
                        end
                    else
                        second_random_te1 = random_te1[2]
                    end
                end

                ############################## Second TE tuple
                silencing_status2 = false

                ####### Replication and/or silencing of TE1
                first_replicated2 = false 

                if random_te2[3] == true
                    local first_random_te2 = random_te2[1]
                    silencing_status2 = true
                else
                    ####### use the poisson process for replication
                    first_random_te_rep2 = random_te2[1] + tau2[2] - tau4[2]
                    if true
                        first_replicated2 = true
                        if first_replicated2 && rand() < silencing_rate
                            first_random_te2 = max(first_random_te_rep2 + 0.0, 0)
                            silencing_status2 = true
                        else
                            first_random_te2 = max(first_random_te_rep2, 0)
                        end
                    else
                        first_random_te2 = random_te2[1]
                    end
                end

                ####### Replication and/or silencing of TE2
                second_replicated2 = false # TE2

                if silencing_status2 == true
                    local second_random_te2 = random_te2[2]
                    silencing_status2 = true
                else
                    ####### use the poisson process for replication
                    second_random_te_rep2 = random_te2[2] + tau3[2] - tau5[2]
                    if true
                        second_replicated2 = true
                        if second_replicated2 && rand() < silencing_rate 
                            second_random_te2 = max(second_random_te_rep2 + 0.0, 0)
                            silencing_status2 = true
                        else
                            second_random_te2 = max(second_random_te_rep2, 0)
                        end
                    else
                        second_random_te2 = random_te2[2]
                    end
                end
                
                tuple1 = (first_random_te1, second_random_te1, silencing_status1)
                tuple2 = (first_random_te2, second_random_te2, silencing_status2)

                push!(new_te_copy_numbers1, tuple1)
                push!(new_te_copy_numbers2, tuple2)
            end

            ####### Update TE copy numbers for the next generation
            te_copy_numbers1 = new_te_copy_numbers1
            te_copy_numbers2 = new_te_copy_numbers2
        end

        ####### Specify absolute path for CSV file
        csv_file_path = "/path/te_$(silencing_rate)_$(i)_data.csv"

        ####### Create a dictionary for writing to CSV
        data_dict = Dict("X_1" => [t[1] for t in te_copy_numbers1],
                         "X_2" => [t[1] for t in te_copy_numbers2],
                         "Y_1" => [t[2] for t in te_copy_numbers1],
                         "Y_2" => [t[2] for t in te_copy_numbers2])

        ####### Write data to CSV file
        CSV.write(csv_file_path, data_dict)
    end
end

####### Parameters
alpha = 2/10
beta1 = 1/10
beta2 = 0.001/10
delta1 = 0.5/10
delta2 = 1.5/10
dt = 1.0
V = 500.0
population_size = 1000
generations = 2000

####### Call the function for each silencing rate
for rate in ARGS
    simulate_TE_dynamics(population_size, generations, alpha, beta1, beta2, delta1, delta2, dt, V, parse(Float64, rate))
end

####### silencing_rates = [0, 0.000001, 0.000002, 0.000003, 0.000004, 0.000005, 0.000006, 0.000007, 0.000008, 0.000009, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.0007, 0.00008, 0.0009]
