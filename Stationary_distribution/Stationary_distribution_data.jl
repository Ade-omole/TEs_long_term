
using Random
using Plots
using Pkg
Pkg.add("DataFrames")
Pkg.add("CSV")
using DataFrames
using CSV

function gillespie_TEs_fixed_step(total_time, rates, time_step)
    time = 0.0
    T, R, H, V = 200, 0, 200, 500
    events = [(time, T, R, H)]

    next_record_time = time_step

    while time < total_time
        propensities = [
            rates[1] * T,
            rates[2] * R * V,
            rates[3] * R * H,
            rates[4] * R * V,
            rates[5] * T,
            rates[6] * H,
        ]
        
        total_propensity = sum(propensities)
        time += -log(rand()) / total_propensity
        if time > total_time
            break
        end
        
        reaction_choice = rand() * total_propensity
        cum_propensity = cumsum(propensities)
        
        if reaction_choice <= cum_propensity[1]
            R += 1
        elseif reaction_choice <= cum_propensity[2]
            R -= 1; T += 1
        elseif reaction_choice <= cum_propensity[3]
            R -= 1; H += 1
        elseif reaction_choice <= cum_propensity[4]
            R -= 1
        elseif reaction_choice <= cum_propensity[5]
            T -= 1
        else
            H -= 1
        end
        
        if time >= next_record_time
            push!(events, (time, T, R, H))
            next_record_time += time_step
        end
    end 
    
    return events
end

# Parameters
total_time = 1100000.0
rates = [2.0, 1.0, 0.5, 1.5, 0.5, 0.5]
# rates = [0.5, 0.1, 0.2, 0.1, 0.2, 0.2]
# rates = [1.1, 0.1, 0.2, 0.1, 0.5, 0.5]
# rates = [0.55, 0.3, 0.4, 0.1, 0.2, 0.2]
time_step = 100.0

# Measure the time taken to run the simulation
@time events = gillespie_TEs_fixed_step(total_time, rates, time_step)

# Extract data for plotting
time_points = [event[1] for event in events]
T_counts = [event[2] for event in events]
H_counts = [event[4] for event in events]

# Write the results to a CSV file
df = DataFrame(time=time_points, x=T_counts/500, y=H_counts/500)
CSV.write("/Users/adekanmiomole/Julia_programming/gillespie_data_10000_distribution.csv", df)