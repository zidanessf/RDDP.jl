using CPLEX,Suppressor,PowerModels,CSV,DataFrames,Interpolations
include("../src/RDDP.jl")
using .RDDP
# read input data
const N_stage,K_seg = 24,5
silence()
case = PowerModels.parse_file("datasets/case6ww.m")
pm = instantiate_model(case, ACPPowerModel, PowerModels.build_opf)
case_dict = pm.ref[:it][:pm][:nw][0]
wind_power = CSV.read("datasets/wind_2020.csv",DataFrame)[1:24,:]
load = CSV.read("datasets/load_data.csv",DataFrame)[1:24,:load]
gens = keys(case_dict[:gen])
RU,RD = 1,1
loadmax = sum(max(case_dict[:gen][gen]["pmax"],0) for gen in gens)
load = load * loadmax/maximum(load)
r = 0.6 * loadmax/maximum(wind_power[!,:wf1])
alpha = 0.25
wind_ipt1 = extrapolate(interpolate(wind_power[:,:wf1],BSpline(Linear())),Flat())
Pw_max = [r*(1+alpha)*wind_ipt1(t*24/N_stage) for t in 1:N_stage]
Pw_min = [r*(1-alpha)*wind_ipt1(t*24/N_stage) for t in 1:N_stage]
Pw = [r*wind_ipt1(t*24/N_stage) for t in 1:N_stage]
ess = [1,2]
Emax,Emin,Pmax,E0 = [1,1],[0,0],[0.5,0.5],[0,0]
msro = RDDP.buildMultiStageRobustModel(
    N_stage = N_stage,
    optimizer = CPLEX.Optimizer,
    MaxIteration = 100,
    MaxTime = 60,
    Gap = 0.01,
    use_maxmin_solver = false
) do ro::JuMP.Model,t
    # defining operator's strategy
    @variable(ro,Pg[gen in gens],RDDP.State,initial_value = case_dict[:gen][gen]["pmin"])
    @variable(ro,pw_max,RDDP.Uncertain,lower_bound=Pw_min[t],upper_bound=Pw_max[t])
    @variable(ro,pw,lower_bound=0)
    @variable(ro,Pgcost[gens],lower_bound=0)
    @variable(ro,Ses[e in [1,2]],RDDP.State,initial_value = E0[e])
    @variable(ro,s1,lower_bound=0)
    @variable(ro,s2,lower_bound=0)
    @variable(ro,Pes[e in [1,2]],lower_bound=-Pmax[e],upper_bound=Pmax[e])
    @constraint(ro,pw <= pw_max)
    for gen in gens # constraints of thermal units
        @constraint(ro,Pg[gen].out >= case_dict[:gen][gen]["pmin"]) # power output upper limit 
        @constraint(ro,Pg[gen].out <= case_dict[:gen][gen]["pmax"]) # power output lower limit 
        @constraint(ro,Pg[gen].out - Pg[gen].in <= RU * case_dict[:gen][gen]["pmax"]) # ramp up limit
        @constraint(ro,Pg[gen].in - Pg[gen].out <= RD * case_dict[:gen][gen]["pmax"]) # ramp down limit
        if length(case_dict[:gen][gen]["cost"]) == 3 # piecewise linear cost 
            for k = 1:K_seg
                pk = case_dict[:gen][gen]["pmax"]*k/K_seg
                costk = (case_dict[:gen][gen]["cost"][2])*pk + (case_dict[:gen][gen]["cost"][1])*pk^2
                ratek = 2*case_dict[:gen][gen]["cost"][1]*pk
                @constraint(ro,Pgcost[gen] >= 24/N_stage * (costk + ratek*(Pg[gen].out - pk)))
            end
        elseif length(case_dict[:gen][gen]["cost"]) == 2
            @constraint(ro,Pgcost[gen] == 24/N_stage * case_dict[:gen][gen]["cost"][1]*Pg[gen].out)
        else
            error("cost term must be 2 or 3")
        end
    end
    for e in ess # constraints of storages
        @constraint(ro,Ses[e].out == Ses[e].in - Pes[e]) #state transition
        @constraint(ro,Ses[e].out <= Emax[e])
        @constraint(ro,Ses[e].out >= Emin[e])
    end
    @constraint(ro,sum(Pg[gen].out for gen in gens) + sum(Pes) + pw + s1 - s2 == load[t]) # power balance
    @objective(ro,Min,sum(Pgcost) + 10000*s1 + 10000*s2 + 2000 * (pw_max - pw))
    # defining attacker's strategy
end
RDDP.train(msro)