include("../src/RDDP.jl")
using .RDDP
N_stage = 24
df = CSV.read("datasets/input.csv",DataFrame)
load,others,gas_plan = df[!,:统调负荷],df[!,:除燃气外],df[!,:燃气]
Pmax,Pmin = 400,200
N_gas = 20
is_ready = [x <= N_gas*2/3 for x in 1:N_gas]
msro = RDDP.buildMultiStageRobustModel(
    N_stage = N_stage,
    optimizer = CPLEX.Optimizer,
    MaxIteration = 100,
    MaxTime = 60,
    Gap = 0.01,
    use_maxmin_solver = false
) do ro::JuMP.Model,t
    @variable(ro,ug[gen in 1:N_gas],Bin,RDDP.State,initial_value=1)
end