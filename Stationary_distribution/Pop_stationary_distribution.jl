using Random, Distributions, CSV, DataFrames

######## Constants
chromosome_length = 5000  ######## Length of each haploid chromosome (loci)
num_autonomous_TEs = 200  ######## Initial number of autonomous TEs ("X") per chromosome
num_nonautonomous_TEs = 200  ######## Initial number of nonautonomous TEs ("Y") per chromosome
population_size = 10000  ######## Total population size (number of individuals)
num_diploid_pairs = 1  ######## Each individual has 1 pair of diploid chromosomes (2 haploid chromosomes total)
generations = 10000  ######## Number of generations to simulate

alpha = 2/100
beta1 = 1/100
beta2 = 0.001/100
delta1 = 0.5/100
delta2 = 1.5/100
dt = 1.0
V = 500

######## Function to create a single haploid chromosome with autonomous and nonautonomous TEs
function create_haploid_chromosome()
    chromosome = fill('.', chromosome_length)  ######## Initialize chromosome with empty sites represented by '.'
    autonomous_sites = rand(1:chromosome_length, num_autonomous_TEs)  ######## Randomly select sites for autonomous TEs ('X')
    for site in autonomous_sites
        chromosome[site] = 'X'  ######## Place 'X' at the selected sites
    end
    available_sites = setdiff(1:chromosome_length, autonomous_sites)  ######## Find sites not occupied by 'X'
    nonautonomous_sites = rand(available_sites, num_nonautonomous_TEs)  ######## Randomly select sites for nonautonomous TEs ('Y')
    for site in nonautonomous_sites
        chromosome[site] = 'Y'  ######## Place 'Y' at the selected sites
    end
    return chromosome
end

######## Function to create an individual with 1 pair of diploid chromosomes
function create_individual()
    chromosomes = [create_haploid_chromosome() for _ in 1:(num_diploid_pairs * 2)]  ######## Create two haploid chromosomes for the individual
    return chromosomes
end

######## Function to create the initial population
function create_population(population_size)
    population = [create_individual() for _ in 1:population_size]  ######## Generate a population of individuals
    return population
end

######## Function to perform crossover between homologous haploid chromosomes
function crossover(chromosome1, chromosome2)
    num_crossovers = rand(Poisson(0.1))  ######## Randomly determine the number of crossovers
    if num_crossovers == 0
        return chromosome1, chromosome2  ######## Return chromosomes unchanged if no crossovers occur
    end

    crossover_positions = sort(rand(1:chromosome_length, num_crossovers))  ######## Get crossover positions in sorted order
    for pos in crossover_positions
        if pos < chromosome_length
            ######## Swap segments after each crossover point
            chromosome1[pos:end], chromosome2[pos:end] = chromosome2[pos:end], chromosome1[pos:end]
        end
    end
    return chromosome1, chromosome2
end

######## Function to update a chromosome by duplicating and deleting TEs based on propensity
function update_chromosome(chromosome)
    ######## Step 1: Count the initial number of 'X' and 'Y' elements
    x_old = count(c -> c == 'X', chromosome)
    y_old = count(c -> c == 'Y', chromosome)

    ######## Create a copy of the chromosome
    new_chromosome = copy(chromosome)

    ######## Step 2: Calculate the number of 'X' and 'Y' to be deleted
    num_x_to_delete = rand(Poisson(delta1 * x_old * dt))
    num_y_to_delete = rand(Poisson(delta1 * y_old * dt))

    ######## Randomly delete elements from the new chromosome while preserving positions
    function delete_random_elements!(element, num_to_delete)
        ######## Find all positions of the given element in the chromosome
        positions = findall(c -> c == element, new_chromosome)
        
        ######## Only proceed if there are positions to delete from
        if isempty(positions)
            return
        end

        ######## Randomly select positions to delete
        delete_positions = rand(positions, min(length(positions), num_to_delete))
         
        ######## Replace the selected positions with '.'
        for pos in delete_positions
            new_chromosome[pos] = '.'
        end
    end
    
    ######## Perform random deletion of 'X' and 'Y' on the chromosome
    delete_random_elements!('X', num_x_to_delete)
    delete_random_elements!('Y', num_y_to_delete)

    ######## Step 3: Replication Process
    av_complex = max((alpha * x_old) / (beta1 + delta2 + beta2 * y_old), 0.0)

    ######## Number of elements to add to the chromosomes
    num_x_to_add = rand(Poisson(beta1 * av_complex * dt))
    num_y_to_add = rand(Poisson(beta2 * av_complex * y_old * dt))
    
    ######## Step 4: Find available positions to place new 'X' and 'Y'
    available_positions = findall(c -> c == '.', new_chromosome)
    occupied_positions = falses(length(new_chromosome))  ######## Array to track occupied positions

    ######## Mark initial available positions
    for pos in available_positions
        occupied_positions[pos] = true
    end

    ######## Step 5: Function to place new elements randomly in available positions
    function place_elements!(element, count)
        for _ in 1:count
            ######## Filter out positions that are still marked as available
            valid_positions = filter(p -> occupied_positions[p], available_positions)
            if !isempty(valid_positions)
                pos = rand(valid_positions)  ######## Select a random valid position
                new_chromosome[pos] = element
                occupied_positions[pos] = false  ######## Mark this position as no longer available
            end
        end
    end

    ######## Place 'X' and 'Y' elements into the available positions
    place_elements!('X', num_x_to_add)
    place_elements!('Y', num_y_to_add)

    return new_chromosome
end

######## Function to perform recombination within a diploid pair of chromosomes
function recombine_diploid_pair(diploid_pair)
    chrom1, chrom2 = diploid_pair
    recombined_chrom1, recombined_chrom2 = crossover(copy(chrom1), copy(chrom2))  ######## Perform crossover and return recombined chromosomes
    return recombined_chrom1, recombined_chrom2
end

######## Function to count the number of X and Y elements in each chromosome of the offspring
function count_TEs_per_chromosome(offspring)
    X_counts = []
    Y_counts = []
    for chromosome in offspring
        push!(X_counts, count(c -> c == 'X', chromosome))  ######## Count 'X' elements
        push!(Y_counts, count(c -> c == 'Y', chromosome))  ######## Count 'Y' elements
    end
    return X_counts, Y_counts
end

######## Function to generate the next generation by performing recombination and update
function generate_next_generation(gen)
    global population
    new_population = []

    generation_X_count = 0  ######## Track total X count for the generation
    generation_Y_count = 0  ######## Track total Y count for the generation
    generation_C_count = gen == 1 ? 200 : 0  ######## Track additional characteristics if needed

    ######## Create a DataFrame to store individual chromosome counts for the final generation
    individual_chromosome_counts = DataFrame(Individual_ID = Int[], X_1 = Int[], X_2 = Int[], Y_1 = Int[], Y_2 = Int[])

    for i in 1:population_size
        ######## Randomly select two parents for recombination
        parent_indices = randperm(population_size)[1:2]
        parent1, parent2 = population[parent_indices[1]], population[parent_indices[2]]

        parent1_recombined = []
        parent2_recombined = []

        crossover_occurred = false  ######## Track if crossover occurred

        ######## Recombine chromosomes from both parents
        recombined_chrom1, recombined_chrom2 = recombine_diploid_pair((parent1[1], parent1[2]))
        push!(parent1_recombined, recombined_chrom1)
        push!(parent1_recombined, recombined_chrom2)

        recombined_chrom3, recombined_chrom4 = recombine_diploid_pair((parent2[1], parent2[2]))
        push!(parent2_recombined, recombined_chrom3)
        push!(parent2_recombined, recombined_chrom4)

        ######## Check if crossover occurred
        if recombined_chrom1 != parent1[1] || recombined_chrom2 != parent1[2]
            crossover_occurred = true
        end
        if recombined_chrom3 != parent2[1] || recombined_chrom4 != parent2[2]
            crossover_occurred = true
        end
    
        offspring = []
        ######## Form offspring by selecting one chromosome from each recombined or original parent
        if crossover_occurred
            chrom_from_parent1 = rand(1:2)
            chrom_from_parent2 = rand(1:2)
            push!(offspring, parent1_recombined[chrom_from_parent1])
            push!(offspring, parent2_recombined[chrom_from_parent2])
        else
            chrom_from_parent1 = rand(1:2)
            chrom_from_parent2 = rand(1:2)
            push!(offspring, parent1[chrom_from_parent1])
            push!(offspring, parent2[chrom_from_parent2])
        end

        ######## Update the chromosomes of the offspring
        updated_offspring = [update_chromosome(chromosome) for chromosome in offspring]
        X_counts, Y_counts = count_TEs_per_chromosome(updated_offspring)

        ######## Accumulate generation-level counts
        generation_X_count += sum(X_counts)
        generation_Y_count += sum(Y_counts)

        ######## Add individual counts to the DataFrame
        push!(individual_chromosome_counts, (i, X_counts..., Y_counts...))

        push!(new_population, updated_offspring)
    end

    ######## Update the global population for the next generation
    population = new_population

    ######## Save the individual chromosome counts for the final generation
    if gen == generations
        CSV.write("/paths/pop_distribution_0.1_recombination.csv", individual_chromosome_counts)
        println("Done")
    end

    ######## Create a DataFrame to store generation-level counts
    generation_data = DataFrame(generation = [gen], Total_X = [generation_X_count], Total_Y = [generation_Y_count], Total_C = [generation_C_count])

    ######## Write generation-level counts to a CSV file for TEs dynamics
    csv_file_path = "/paths/total_count_0.1_recombination.csv"
    if gen == 1  ######## Create the file and add the header for the first generation
        CSV.write(csv_file_path, generation_data)
    else  ######## Append subsequent generations without headers
        CSV.write(csv_file_path, generation_data, append=true)
    end

    ######## Print the generation-level counts to the console
    println("Generation $gen: Total X = $generation_X_count, Total Y = $generation_Y_count, Total C = $generation_C_count")
end

######## Time the simulation
@time begin
    ######## Create initial population and run the simulation
    population = create_population(population_size)
    for gen in 1:generations
        generate_next_generation(gen)
    end
end
