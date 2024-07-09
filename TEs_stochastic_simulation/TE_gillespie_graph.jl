using Random
using Plots
using DifferentialEquations

####### Function to simulate the Gillespie algorithm for transposable elements (TEs)
function gillespie_TEs(total_time, rates)
    ####### Initialize time and copy numbers for X, C, and Y
    time = 0.0
    X = 200
    C = 200
    Y = 200
    global V = 500  ####### Volume or scaling factor for propensities

    ####### List to store events as tuples of (time, X, Y, C)
    events = [(time, X, Y, C)]

    ####### Main simulation loop
    while time < total_time
        ####### Calculate propensities for each reaction
        propensities = [
            rates[1] * X,        ####### Reaction 1: X -> X + C
            rates[2] * C * V,    ####### Reaction 2: C -> X
            rates[3] * C * Y,    ####### Reaction 3: C + Y -> 2Y
            rates[4] * C * V,    ####### Reaction 4: C -> ∅
            rates[5] * X,        ####### Reaction 5: X -> ∅
            rates[6] * Y         ####### Reaction 6: Y -> ∅
        ]
        
        total_propensity = sum(propensities)  ####### Sum of all propensities
        reaction_time = -log(rand()) / total_propensity  ####### Time until next reaction
        time += reaction_time  ####### Increment the time

        if time > total_time
            break
        end
        
        ####### Determine which reaction occurs
        reaction_choice = rand() * total_propensity
        if reaction_choice <= propensities[1]
            C += 1
        elseif reaction_choice <= propensities[1] + propensities[2]
            C -= 1
            X += 1
        elseif reaction_choice <= propensities[1] + propensities[2] + propensities[3]
            C -= 1
            Y += 1
        elseif reaction_choice <= propensities[1] + propensities[2] + propensities[3] + propensities[4]
            C -= 1
        elseif reaction_choice <= propensities[1] + propensities[2] + propensities[3] + propensities[4] + propensities[5]
            X -= 1
        else
            Y -= 1
        end
        
        ####### Store the event
        push!(events, (time, X, Y, C))
    end 
    
    return events
end

####### Parameters
total_time = 200.0
####### Chosen rates  α,β1,β2,δ2,δ1,δ1
rates = [2.0, 1.0, 0.5, 1.5, 0.5, 0.5] 

####### Run the Gillespie simulation
events = gillespie_TEs(total_time, rates)

####### Extract data for plotting
time_points = [event[1] for event in events]
X_counts = [event[2] for event in events]
Y_counts = [event[3] for event in events]
C_counts = [event[4] for event in events]

####### Function to define the ODE system for TEs
function te_rre(du, u, p, t)
    α, β1, β2, δ1, δ2 = p
    du[1] = u[1] * (β1 * α / (β1 + δ2 + β2 * u[2]) - δ1)  ####### dX/dt
    du[2] = u[2] * (β2 * α * u[1] / (β1 + δ2 + β2 * u[2]) - δ1)  ####### dY/dt
end

tspan = (0.0, total_time)  ####### Time span for the ODE solver
u0 = [0.4, 0.4]  ####### Initial conditions for X and Y

####### Chosen rates
p = (2, 1, 0.5, 0.5, 1.5) 

####### Define the ODE problem
TE_prob = ODEProblem(te_rre, u0, tspan, p)
sol = solve(TE_prob, abstol=1e-8, alg_hints=[:stiff])  ####### Solve the ODE

####### Extract results for ODE
t = sol.t
u = hcat(sol.u...)'  ####### Concatenate solution values horizontally

####### Plot results
plot(time_points, X_counts/V, label=false)  ####### Plot Gillespie X counts
plot!(time_points, C_counts/V, label=false)  ####### Plot Gillespie C counts
plot!(time_points, Y_counts/V, label=false)  ####### Plot Gillespie Y counts
plot!(t, u, lw=2, color=[:blue :green], label=false)  ####### Plot ODE solutions

xlabel!("Time (generation)") 
ylabel!("Copy number (scaled)")
####### title!("Dynamics of TEs")
savefig("/path/TE_gillespie_graph.pdf")  ####### Save the plot to a PDF file
