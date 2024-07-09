using Random
using Plots
using Pkg
Pkg.add("DataFrames")
Pkg.add("CSV")
using DataFrames
using CSV

####### Function to simulate the Gillespie algorithm for transposable elements (TEs) with fixed time steps
function gillespie_TEs_fixed_step(total_time, rates, time_step)

    ####### Initialize time and copy numbers for X, C, Y, and V
    time = 0.0
    X, C, Y, global V = 200, 0, 200, 500

    ####### List to store events as tuples of (time, X, C, Y)
    events = [(time, X, C, Y)]

    next_record_time = time_step

    ####### Main simulation loop
    while time < total_time

        ####### Calculate propensities for each reaction
        propensities = [
            rates[1] * X,        ####### Reaction 1: X -> X + C
            rates[2] * C * V,    ####### Reaction 2: C -> X
            rates[3] * C * Y,    ####### Reaction 3: C + Y-> 2Y
            rates[4] * C * V,    ####### Reaction 4: C -> ∅
            rates[5] * X,        ####### Reaction 5: X -> ∅
            rates[6] * Y         ####### Reaction 6: Y -> ∅
        ]
        
        total_propensity = sum(propensities)  ####### Sum of all propensities
        time += -log(rand()) / total_propensity  ####### Time until next reaction
        if time > total_time
            break
        end
        
        ####### Determine which reaction occurs
        reaction_choice = rand() * total_propensity
        cum_propensity = cumsum(propensities)
        
        if reaction_choice <= cum_propensity[1]
            C += 1
        elseif reaction_choice <= cum_propensity[2]
            C -= 1; X += 1
        elseif reaction_choice <= cum_propensity[3]
            C -= 1; Y += 1
        elseif reaction_choice <= cum_propensity[4]
            C -= 1
        elseif reaction_choice <= cum_propensity[5]
            X -= 1
        else
            Y -= 1
        end
        
        ####### Record the event at fixed time steps
        if time >= next_record_time
            push!(events, (time, X, C, Y))
            next_record_time += time_step
        end
    end 
    
    return events
end

####### Parameters
total_time = 1100000.0
rates = [2.0, 1.0, 0.5, 1.5, 0.5, 0.5]
time_step = 100.0

####### Measure the time taken to run the simulation
@time events = gillespie_TEs_fixed_step(total_time, rates, time_step)

####### Extract data for plotting
time_points = [event[1] for event in events]
X_counts = [event[2] for event in events]
Y_counts = [event[4] for event in events]

####### Write the results to a CSV file
df = DataFrame(time=time_points, x=X_counts/V, y=Y_counts/V)
CSV.write("/path/data_distribution.csv", df)
