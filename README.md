# Robust Dual Dynamic Programming (RDDP)
## Usage
```julia
pkg> add RDDP
```
## Example
```julia
# inventory management under uncertain demand
using RDDP,Gurobi,JuMP
N = 10
msro = RDDP.buildMultiStageRobustModel(
    N_stage = N, # number of stages
    optimizer = Gurobi.Optimizer, # choose an optimizer
    MaxIteration = 100, # maximum epoches
    MaxTime = 60, # maximum trainning time
    Gap = 0.01, # optimal gap
    use_maxmin_solver = false # using a vertex enumeration by default
) do ro::JuMP.Model,t::Int
    # defining operator's strategy
    @variable(ro,0 <= x <= 5,RDDP.State,initial_value = 0) # state variable (stock level)
    @variable(ro,0 <= u <= 1) # control variable (supply)
    @variable(ro,0 <= ξ[i = 1:2] <= 1,RDDP.Uncertain) # uncertain parameter (demand)
    @constraint(ro,sum(ξ) <= 1) # constraint abount uncertain parameter
    @constraint(ro, x.out == x.in - u + sum(ξ)) # state transition (temporal dependence of stock level)
    @objective(ro,Min,x) # minimize the stock
end
RDDP.train(msro)
```
## Citing RDDP.jl
If you use RDDP.jl, we ask that you please cite the following:
```
@article{rddp.jl,
  title={Enhancing the Flexibility of Storage Integrated Power System by Multi-stage Robust Dispatch},
  author={Shi, Yunhui and Guo, Chuangxin and Dong, Shufeng and Chen, Zhe and Wang, Luyu},
  journal={IEEE Transactions on Power Systems},
  year={2020},
  publisher={IEEE}
}

@article{georghiou2019robust,
  title={Robust dual dynamic programming},
  author={Georghiou, Angelos and Tsoukalas, Angelos and Wiesemann, Wolfram},
  journal={Operations Research},
  volume={67},
  number={3},
  pages={813--830},
  year={2019},
  publisher={INFORMS}
}
```
