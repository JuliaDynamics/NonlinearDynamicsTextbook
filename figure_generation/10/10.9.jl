using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using Agents, LightGraphs, DataFrames
# import GLMakie
# using GraphMakie.NetworkLayout


n0 = 100
seed = 42
s = 80 # steps to run
er = n -> erdos_renyi(n, 200; seed)
ba = n -> barabasi_albert(n, 3, 2; complete=true, seed)
ws = n -> watts_strogatz(n, 4, 0.1; seed)

graphs = [g(n0) for g in (er, ba, ws)]

# g = graphs[i](n0)

# %% First version: one agent per node, no moving

# model parameters: 
τ = 0.1  # transmission probability
# d = 0.02 # death probability? Should I?
p = 14   # period to recover

# Create an ABM on a graph `g`
@agent Human GraphAgent begin
    status::Symbol # :S, :I, :R
    days_infected::Int
end
function init_sir_abm(g; τ = 0.1, p = 14, d = 0.02)
    abm = ABM(Human, GraphSpace(g); properties = (τ = τ, p = p, d = d))
    # Add one agent to each node, with infected agent at node 1
    starting_status(i) = i == 1 ? (:I, 1) : (:S, 0)
    fill_space!(abm, starting_status)
    return abm
end
# Define rules of the ABM simulation
function agent_step!(human, abm) # for each agent at each step
    human.status ≠ :I && return # do nothing if not infected
    # Infect all susceptible neighbors with probability τ
    for friend in nearby_agents(human, abm)
        if friend.status == :S && rand() < abm.τ
            friend.status = :I
        end
    end
    # Check if human recovers or dies
    human.days_infected += 1
    if human.days_infected ≥ abm.p
        human.status = :R
    end
    # elseif (rand() < abm.d)
    #     kill_agent!(human, abm)
    # end
end
# run simulation 
g = graphs[1]
abm = init_sir_abm(g)
infected(abm) = count(a.status == :I for a in allagents(abm))
recovered(abm) = count(a.status == :R for a in allagents(abm))
susceptible(abm) = count(a.status == :S for a in allagents(abm))
mdata = [infected, recovered, susceptible]

# infected(x) = count(i == :I for i in x)
# recovered(x) = count(i == :R for i in x)
# adata = [(:status, infected), (:status, recovered)]



names = ("Erdős–Rényi ", "Barabási–Albert", "Watts-Strogatz")
names = ("ER ", "BA", "WS")

fig = figure(figsize = (figx/2, figy))
for (i, g) in enumerate(graphs)
    abm = init_sir_abm(g)
    _, mdf = run!(abm, agent_step!, 100; mdata)
    plot(mdf[!, :infected]; color = "C$(i-1)", label = "$(names[i])")
    plot(mdf[!, :recovered]; color = "C$(i-1)", lw = 1.5, ls = "--")
    plot(mdf[!, :susceptible]; color = "C$(i-1)", lw = 1.5, ls = ":")
end
legend(handlelength=1)
xlabel("time")
fig.tight_layout(pad=0.3)

# %% Ensemble version
using Statistics
fig = figure(figsize = (figx/2, figy))
ax = gca()
for (i, g) in enumerate(graphs)
    _, mdf = ensemblerun!(seed -> init_sir_abm(g), agent_step!, dummystep, 100; mdata, ensemble = 10)
    μdf = combine(groupby(mdf, :step), [:infected, :susceptible, :recovered] .=> mean .=> [:I, :S, :R])
    σdf = combine(groupby(mdf, :step), [:infected, :susceptible, :recovered] .=> std .=> [:I, :S, :R])

    # gdf = groupby(mdf, :ensemble)
    # I = sum(df[!, :infected] for df in gdf) ./ 5
    # R = sum(df[!, :recovered] for df in gdf) ./ 5
    # S = sum(df[!, :susceptible] for df in gdf) ./ 5
    I = μdf[!, :I]
    Iσ = σdf[!, :I]
    ax.plot(I; color = "C$(i-1)", label = "$(names[i])")
    ax.fill_between(0:length(I)-1, I .- Iσ, I .+ Iσ; color = "C$(i-1)", alpha = 0.25)
    # plot(R; color = "C$(i-1)", lw = 1.5, ls = "--")
    # plot(S; color = "C$(i-1)", lw = 1.5, ls = ":")
end
legend(handlelength=1)
xlabel("time (days)")
ylabel("infected agents")
fig.tight_layout(pad=0.3)

wsave(plotsdir("10", "sir_network"), fig)

# %% Second version: pretty much what we have online: many agents per node (city)
# movement between cities and infecting only within node.

nt = 50
seed = 42
s = 80 # steps to run
er = erdos_renyi(nt, 100; seed)
ba = barabasi_albert(nt, 3, 2; complete=true, seed)
ws = watts_strogatz(nt, 4, 0.1; seed)

graphs = (er, ba, ws)

function init_sir_abm_transport(g; τ=0.05, p=14, m = 0.5)
    tprobs = transition_probabilities(g)
    properties = @ntuple tprobs τ p m
    abm = ABM(Human, GraphSpace(g); properties)
    # Add a bunch of agents at each node, set agent0 as infected.
    for i in positions(abm, :random)
        for j in 1:10
            add_agent!(i, abm, :S, 0)
        end
    end
    abm[1].status = :I
    return abm
end

function transition_probabilities(g::AbstractGraph)
    degrees = degree(g)
    distributions = Vector{DiscreteNonParametric}(undef, nv(g))
    for i in 1:nv(g)
        n = neighbors(g, i)
        if isempty(n)
            vals = [i]; probs = [1.0]
        else
            vals = n
            outdtot = sum(degrees[j] for j in n)
            probs = [degrees[j]/outdtot for j in n]
        end
        distributions[i] = DiscreteNonParametric(vals, probs)
    end
    return distributions
end

function agent_step!(human, abm)
    # Travelling
    if rand(abm.rng) < abm.m
        newnode = rand(abm.rng, abm.tprobs[human.pos])
        move_agent!(human, newnode, abm)
    end
    # Infecting
    human.status ≠ :I && return # do nothing if not infected
    for person in agents_in_position(human, abm) # nearby_agents(human, abm, 0)
        if person.status == :S && rand() < abm.τ
            person.status = :I
        end
    end
    # Check if human recovers or dies
    human.days_infected += 1
    if human.days_infected ≥ abm.p
        human.status = :R
    end
end

fig = figure(figsize = (figx/2, figy))
ax = gca()
for (i, g) in enumerate(graphs)
    _, mdf = ensemblerun!(seed -> init_sir_abm_transport(g), agent_step!, dummystep, 100; mdata, ensemble = 10)
    μdf = combine(groupby(mdf, :step), [:infected, :susceptible, :recovered] .=> mean .=> [:I, :S, :R])
    σdf = combine(groupby(mdf, :step), [:infected, :susceptible, :recovered] .=> std .=> [:I, :S, :R])

    # gdf = groupby(mdf, :ensemble)
    # I = sum(df[!, :infected] for df in gdf) ./ 5
    # R = sum(df[!, :recovered] for df in gdf) ./ 5
    # S = sum(df[!, :susceptible] for df in gdf) ./ 5
    I = μdf[!, :I]
    Iσ = σdf[!, :I]
    ax.plot(I; color = "C$(i-1)", label = "$(names[i])")
    ax.fill_between(0:length(I)-1, I .- Iσ, I .+ Iσ; color = "C$(i-1)", alpha = 0.25)
    # plot(R; color = "C$(i-1)", lw = 1.5, ls = "--")
    # plot(S; color = "C$(i-1)", lw = 1.5, ls = ":")
end
legend(handlelength=1)
xlabel("time (days)")
ylabel("infected agents")
# fig.tight_layout(pad=0.3)