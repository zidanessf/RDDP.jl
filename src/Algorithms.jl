total_solves = 0
number_of_reject = 0
sceanio_not_change_num = 0
FollowerValuePrev = 0
mutable struct MultiStageRobustModel
    upper::Array{JuMP.Model}
    lower::Array{JuMP.Model}
    uncertain::Array{JuMP.Model}
    config::Dict
end

function add_upper_bound(msro,t,n_iter)
    dual_of_obj = @variable(msro.upper[t],lower_bound=0)
    if n_iter == 1
        @constraint(msro.upper[t],sum_y,dual_of_obj == 1)
    end
    v_upper = objective_value(msro.upper[t+1])
    # v_upper = (objective_value(msro.lower[t+1]) + objective_value(msro.upper[t+1]))/2
    #modify objective
    set_normalized_coefficient(msro.upper[t][:upper_bound],dual_of_obj,-v_upper)
    # modify constraint  ∑y_k == 1
    set_normalized_coefficient(msro.upper[t][:sum_y],dual_of_obj,1)
    # add constraint
    # x_i - ∑y_k*x_{ki} + τu_i - τl_i == 0 for x_i in states variables
    for i in 1:length(msro.upper[t][:sum_states])
        set_normalized_coefficient(msro.upper[t][:sum_states][i],dual_of_obj, - value(msro.upper[t+1][:state][i].in))
    end
end

function add_lower_bound(msro,t)
    lower_cut = AffExpr()
    for i in 1:length(msro.lower[t][:state])
        p = dual(FixRef(msro.lower[t+1][:state][i].in))
        msro.lower[t][:penalty] = max(msro.lower[t][:penalty],1.01*abs(p))#dynamic method
        add_to_expression!(lower_cut,p * (msro.lower[t][:state][i].out - value(msro.lower[t+1][:state][i].in)))
    end
    v_lower = objective_value(msro.lower[t+1])
    add_to_expression!(lower_cut,v_lower)
    @constraint(msro.lower[t],msro.lower[t][:cost_to_go] >= lower_cut)
end
function solve_max_min(_max::JuMP.Model,_min::JuMP.Model,binder_max,binder_min;lipshitz_constant=100000,mu_upper_bound=1e7,max_iteration=20)
    try
        @assert _max[:mu_tag]
    catch
        @variable(_max,μ,upper_bound=mu_upper_bound)
        @variable(_max,λ[i=1:length(binder_max)],lower_bound=-lipshitz_constant,upper_bound=lipshitz_constant)
        @objective(_max,Max,μ + sum(λ.*binder_max))
        _max[:mu_tag] = true
    end
    # delete exisiting cuts
    try
        while !isempty(_max[:valid_cut])
            con = pop!(_max[:valid_cut])
            delete(_max,con)
        end
    catch #_max[:valid_cut] not exisit
        _max[:valid_cut] = []
    end
    num = 1
    cut_bound = 0
    while true
        Suppressor.@suppress_out optimize!(_max)
        UB = objective_value(_max)
        binder_value = value.(binder_max)
        for i in 1:length(binder_value)
            fix(binder_min[i],binder_value[i];force=true)
        end
        Suppressor.@suppress_out optimize!(_min)
        LB = objective_value(_min)
        @constraint(_max,_max[:μ] + sum(_max[:λ] .* binder_value) <= LB)
        for i in 1:length(binder_value)
            cut_bound = max(cut_bound,abs(dual(FixRef(binder_min[i]))))
        end
        num += 1
        # @info(UB,LB)
        if abs(UB-LB)/LB <= 0.01
            break
        end
        if num >= max_iteration
            @warn("maxmin algorithm did not converge, the current gap is $(abs(UB-LB)/LB)")
            break
        end
        for i in 1:length(binder_value)
            set_upper_bound(_max[:λ][i],min(lipshitz_constant,1.5 * cut_bound))
            set_lower_bound(_max[:λ][i],- min(lipshitz_constant,1.5 * cut_bound))
        end
    end
end
function ForwardPassPrimal(msro,start,stop)
    additional = Dict()
    additional[:gaphourly] = []
    # additional[:upper] = []
    gap = 0
    differ = 0
    sceanio_not_change = true
    for t in start:stop #前推步骤
        # fix the dayahead decision
        for i in 1:length(msro.lower[t][:state])
            if t > 1
                fix(msro.upper[t][:state][i].in,value(msro.lower[t-1][:state][i].out))
                fix(msro.lower[t][:state][i].in,value(msro.lower[t-1][:state][i].out))
            end
        end
        # if sceanio_not_change_num >= 20
        #     wst_vertex = msro.upper[t][:wst_vertex]
        # else
        #     wst_vertex = vertex_enumeration_primal(msro,t)
        # end
        if msro.config["use_maxmin_solver"]
            solve_max_min(msro.uncertain[t],msro.upper[t],msro.uncertain[t][:uncertain],msro.upper[t][:uncertain],max_iteration=20,lipshitz_constant=2 * msro.config["initial_penalty"])
            msro.upper[t][:wst_case] = value.(msro.uncertain[t][:uncertain])
        else
            wst_obj,wst_case = 0,[]
            for wst in msro.upper[t][:vertex]
                for i in 1:length(wst)
                    fix(msro.upper[t][:uncertain][i],wst[i]; force=true)
                end
                Suppressor.@suppress_out optimize!(msro.upper[t])
                @assert termination_status(msro.upper[t]) == MOI.OPTIMAL
                if objective_value(msro.upper[t]) > wst_obj
                    wst_obj = objective_value(msro.upper[t])
                    wst_case = wst
                end
            end
            msro.upper[t][:wst_case] = wst_case
        end
        for i in 1:length(msro.upper[t][:wst_case])
            fix(msro.lower[t][:uncertain][i],msro.upper[t][:wst_case][i]; force=true)
            fix(msro.upper[t][:uncertain][i],msro.upper[t][:wst_case][i]; force=true)
        end
        Suppressor.@suppress_out optimize!(msro.lower[t])
        Suppressor.@suppress_out optimize!(msro.upper[t])
        @assert termination_status(msro.lower[t]) == MOI.OPTIMAL
        @assert termination_status(msro.upper[t]) == MOI.OPTIMAL
        global total_solves
        total_solves += 2
        if t < length(msro.upper)
            gapt  = value((msro.upper[t][:cost_to_go]) - value(msro.lower[t][:cost_to_go]))/value(msro.upper[t][:cost_to_go])
        else
            gapt = 0
        end
        push!(additional[:gaphourly],round(gapt;digits=3))
    end
    additional[:UpperBound] = objective_value(msro.upper[1])
    additional[:LowerBound] = objective_value(msro.lower[1])
    additional[:Gap] = (additional[:UpperBound] - additional[:LowerBound])/additional[:UpperBound]
    return additional
end
function BackwardPassPrimal(msro,N_ITER,start,stop)
    for t in [stop-x+start for x in start+1:stop]#回代步骤
        # **update the overestimator**
            # **solve the updated upper problem**
            # solve_max_min(msro.upper[t+1][:uncertain_problem],msro.upper[t+1],msro.upper[t+1][:uncertain_problem][:uncertain],msro.upper[t+1][:uncertain])
            # wst_case = value.(msro.upper[t+1][:uncertain_problem][:uncertain])
        if msro.config["use_maxmin_solver"]
            solve_max_min(msro.uncertain[t+1],msro.upper[t+1],msro.uncertain[t+1][:uncertain],msro.upper[t+1][:uncertain],max_iteration=20,lipshitz_constant=2 * msro.config["initial_penalty"])
            msro.upper[t+1][:wst_case] = value.(msro.uncertain[t+1][:uncertain])
        else
            wst_obj,wst_case = 0,[]
            for wst in msro.upper[t+1][:vertex]
                for i in 1:length(wst)
                    fix(msro.upper[t+1][:uncertain][i],wst[i]; force=true)
                end
                Suppressor.@suppress_out optimize!(msro.upper[t+1])
                @assert termination_status(msro.upper[t+1]) == MOI.OPTIMAL
                if objective_value(msro.upper[t+1]) >= wst_obj
                    wst_obj = objective_value(msro.upper[t+1])
                    wst_case = wst
                end
            end
            msro.upper[t+1][:wst_case] = wst_case
        end
        # wst_vertex = msro.upper[t+1][:wst_vertex]
        for i in 1:length(msro.upper[t+1][:wst_case])
            fix(msro.lower[t+1][:uncertain][i],msro.upper[t+1][:wst_case][i]; force=true)
            fix(msro.upper[t+1][:uncertain][i],msro.upper[t+1][:wst_case][i]; force=true)
        end
        # **solve the updated lower problem**
        Suppressor.@suppress_out optimize!(msro.lower[t+1])
        Suppressor.@suppress_out optimize!(msro.upper[t+1])
        @assert termination_status(msro.upper[t+1]) == MOI.OPTIMAL
        @assert termination_status(msro.lower[t+1]) == MOI.OPTIMAL
        global total_solves
        total_solves += 2
        # intradayMaxtobj = value(msro.upper[t][:cost_to_go])
        if objective_value(msro.upper[t+1]) <= value(msro.upper[t][:cost_to_go]) ||  N_ITER == 1
            add_upper_bound(msro,t,N_ITER)
        end
        # add_upper_bound(msro,t,N_ITER)
        # add_upper_bound(msro,t,N_ITER)
        # **update the underestimator**
        if (objective_value(msro.lower[t+1]) >= value(msro.lower[t][:cost_to_go])) ||  N_ITER == 1
            add_lower_bound(msro,t)
        end
        # @info(objective_value(msro.lower[t+1]),value(msro.lower[t][:cost_to_go]))
        # add_lower_bound(msro,t)
    end
    for t in start:stop-1
        for i in 1:length(msro.upper[t][:state]) 
            set_normalized_coefficient(msro.upper[t][:upper_bound],msro.upper[t][:τu][i],-(msro.lower[t][:penalty] + 10000*exp(-N_ITER)))
            set_normalized_coefficient(msro.upper[t][:upper_bound],msro.upper[t][:τl][i],-(msro.lower[t][:penalty] + 10000*exp(-N_ITER)))
        end
    end
end


function mature_stage(msro,additional)
    if minimum(additional[:gaphourly][1:end-1]) < 0.01 && maximum(additional[:gaphourly][1:end-1]) < 0.2 
        stop = findfirst(x->x<=0.01,additional[:gaphourly][1:end-1])
    else 
        stop = length(msro.lower)
    end
    # @info("stop = $stop")
    return stop
end

function train(msro)
    n_iter = 0
    global total_solves,sceanio_not_change_num
    total_solves,sceanio_not_change_num,gap_temp,no_improvement = 0,0,9999,0
    start,stop = 1,length(msro.lower)
    additional = Dict()
    solution_status = DataFrames.DataFrame(Iteration=[],UpperBound=[],LowerBound=[],Gap=[],Time=[],TotalSolves=[]) 
    print_banner(stdout)
    print_iteration_header(stdout)
    t1 = time()
    while true
        n_iter += 1
        additional = ForwardPassPrimal(msro,start,stop)
        if n_iter >= 5
            stop = mature_stage(msro,additional)
            stop = max(stop,4)
        end
        # ______________________________________________________________________________________________________________
        
        additional[:Iteration] = n_iter
        additional[:TotalSolves] = total_solves
        t2 = time()
        additional[:Time] = t2 - t1
        pop!(additional,:gaphourly)
        print_iteration(stdout,additional)
        if n_iter >= 2
            push!(solution_status,additional)
        end
        if n_iter >= msro.config["MaxIteration"] || additional[:Time] >= msro.config["MaxTime"] || no_improvement >= msro.config["no_improvement"] 
            printstyled("Fail to converge. ";color=:red)
            print("$n_iter iterations in $(additional[:Time]) seconds. \n")
            break
        end
        if n_iter >= 10 && (abs(additional[:Gap]) <= msro.config["Gap"] || stop == 1)
            printstyled("Converged. ";color=:green)
            print("$n_iter iterations in $(round(additional[:Time];digits=3)) seconds. \n")
            break
        end
        Suppressor.@suppress_out BackwardPassPrimal(msro,n_iter,start,stop)
    end
    return solution_status
end
function buildMultiStageRobustModel(creator::Function;N_stage::Int,optimizer,MaxIteration=100,MaxTime=3600,Gap=0.01,initial_penalty = 1e6,use_maxmin_solver=false)
    msro = Suppressor.@suppress_out MultiStageRobustModel(
        [JuMP.Model(with_optimizer(optimizer)) for t in 1:N_stage],
        [JuMP.Model(with_optimizer(optimizer)) for t in 1:N_stage],
        [JuMP.Model(with_optimizer(Ipopt.Optimizer)) for t in 1:N_stage],
        Dict([("MaxIteration",MaxIteration),("MaxTime",MaxTime),("Gap",Gap),("no_improvement",100),("initial_penalty",initial_penalty),("use_maxmin_solver",use_maxmin_solver)]))
    for t in 1:N_stage # construct upper bound problem
        m = msro.upper[t]
        m[:state] = []
        m[:initial_value] = []
        m[:uncertain] = []
        # set_parameter(m[:uncertain_problem],"NonConvex",2)
        creator(m,t)
        if t == 1
            for i in 1:length(m[:state])
                fix(m[:state][i].in,m[:initial_value][i])
            end
        end
        if t < N_stage
            @variable(m,cost_to_go,lower_bound=0)
            @variable(m,τu[1:length(m[:state])],lower_bound=0)
            @variable(m,τl[1:length(m[:state])],lower_bound=0)
            @constraint(m,sum_states[i=1:length(m[:state])],m[:state][i].out + τu[i] - τl[i] == 0)
            @constraint(m,upper_bound,cost_to_go >= sum(τu[i]*initial_penalty for i in 1:length(m[:state])) + sum(τl[i]*initial_penalty for i in 1:length(m[:state])))
            @objective(m,Min,objective_function(m) + cost_to_go)
            m[:cost_to_go_now] = cost_to_go
        end
        m[:converging] = false
        m[:penalty] = 0.01
    end
    for t in 1:N_stage # construct lower bound problem
        m = msro.lower[t]
        m[:state] = []
        m[:initial_value] = []
        m[:uncertain] = []
        creator(m,t)
        if t == 1
            for i in 1:length(m[:state])
                fix(m[:state][i].in,m[:initial_value][i])
            end
        end
        if t < N_stage
            @variable(m,cost_to_go,lower_bound=0)
            @objective(m,Min,objective_function(m) + cost_to_go)
            m[:cost_to_go_now] = cost_to_go
        end
        m[:converging] = false
        m[:penalty] = 0.01
    end
    for t in 1:N_stage # construct uncertainty parameter problem
        m = msro.uncertain[t]
        m[:state] = []
        m[:initial_value] = []
        m[:uncertain] = []
        creator(m,t)
        @objective(m,Min,0)
        for v ∈ JuMP.all_variables(m) # delete non-uncertain variables
            if v ∉ m[:uncertain]
                JuMP.delete(m,v)
            end
        end
        for ctype ∈ JuMP.list_of_constraint_types(m)# delete non-uncertain constraints
            for c ∈ JuMP.all_constraints(m,ctype[1],ctype[2])
                try
                    if all(JuMP.normalized_coefficient(c,v) == 0 for v ∈ m[:uncertain])
                        delete(m,c)
                    end
                catch
                end
            end
        end
        m[:valid_cut] = [] # store cuts generated by max-min algorithm
    end
    for t in 1:N_stage # compute vertice
        poly = Polyhedra.polyhedron(msro.uncertain[t], CDDLib.Library())
        msro.upper[t][:vertex] = [v for v in Polyhedra.points(Polyhedra.vrep(poly))]
    end
    # L1_regularization(msro)
    return msro
end
function L1_regularization(msro::MultiStageRobustModel) # not working
    for t in 1:length(msro.lower)
        m = msro.lower[t]
        @variable(m,z_aux[1:length(m[:state])],lower_bound=0)
        @variable(m,z[1:length(m[:state])],lower_bound=0)
        @objective(m,Min,objective_function(m) + 50000*sum(z_aux))
        for i in 1:length(m[:state])
            @constraint(m,z_aux[i] >= z[i] - m[:state][i].in)
            @constraint(m,z_aux[i] >= - z[i] + m[:state][i].in)
        end
    end
end

function evaluate_under_worst_case(msro,worst_case)
    total_cost = 0
    for t in 1:length(worst_case)
        fixprev(msro,t)
        for (idx,gen) in enumerate(keys(msro.case_dict[:windfarm]))
            vertex = [v for v in worst_case]
            windpower = msro.data[t][:wind_power]
            # @assert abs(vertex[idx]) <= 0.5
            fix(msro.lower[t][:Pw_err][gen],vertex[t][idx]*windpower[gen])
        end
        Suppressor.@suppress_out optimize!(msro.lower[t])
        # println(value(m[t][:cost_now]))
        @assert termination_status(msro.lower[t]) == MOI.OPTIMAL #OPTIMAL
        total_cost += value(msro.lower[t][:cost_now])
    end
    return total_cost
    # ForwardPassPrimal(msro,1,length(msro.lower))
    # return sum(value(msro.lower[t][:cost_now]) for t in 1:length(msro.lower))
end
function evaluate_under_worst_case(leader,follower,worst_case)
    total_cost = 0
    for t in 1:length(worst_case)
        fixprev(leader,t)
        for (idx,gen) in enumerate(keys(leader.case_dict[:windfarm]))
            vertex = [v for v in worst_case]
            windpower = leader.data[t][:wind_power]
            # @assert abs(vertex[idx]) <= 0.5
            fix(leader.intraday[t][:Pw_err][gen],vertex[t][idx]*windpower[gen])
        end
        Suppressor.@suppress_out optimize!(leader.intraday[t])
        # println(value(m[t][:cost_now]))
        @assert termination_status(leader.intraday[t]) == MOI.OPTIMAL #OPTIMAL
        total_cost += value(leader.intraday[t][:cost_now])
        fixprev(follower,value.(leader.intraday[t][:binder]),t)
        Suppressor.@suppress_out optimize!(follower.intraday[t])
        total_cost += value(follower.intraday[t][:cost_now])
    end
    return total_cost
end
function evaluate(msro)
    Random.seed!(1235)
    n = 10
    T = length(msro.lower)
    trajectory = []
    vts = msro.vertice
    nk = length(vts[1])
    for i = 1:n
        trajectory_i = [vts[t][randperm(nk)[1]] for t in 1:T]
        push!(trajectory,trajectory_i)
    end
    objectives = []
    Ses = []
    # p = Progress(n)
    worst_index,tmp_worst = 1,0
    for k in 1:n
        # total_cost[k] = value(dayahead[:cost_now])
        m = msro.lower
        total_cost = 0
        for t in 1:T
            fixprev(msro,t)
            for (idx,gen) in enumerate(keys(msro.case_dict[:windfarm]))
                vertex = trajectory[k][t]
                windpower = msro.data[t][:wind_power]
                # @assert abs(vertex[idx])*m[t][:α] <= 0.5
                fix(m[t][:Pw_err][gen],vertex[idx]*m[t][:α]*windpower[gen])
            end
            Suppressor.@suppress_out optimize!(m[t])
            # println(value(m[t][:cost_now]))
            @assert termination_status(m[t]) == MOI.OPTIMAL #OPTIMAL
            total_cost += value(m[t][:cost_now])
        end
        if total_cost > tmp_worst
            tmp_worst = total_cost
            worst_index = k
        end
        push!(objectives,[total_cost,value(m[T][:Ew])/value(m[T][:Ew_max]),value(m[T][:Ses][1]),sum(value(m[t][:load_Cut_cost]) for t in 1:T)])
        push!(Ses,[value(msro.lower[t][:Ses][1]) for t in 1:T])
    end
    trajectory_worst = trajectory[worst_index]
    m = msro.lower
    return [x[1] for x in objectives],[x[2] for x in objectives],[x[3] for x in objectives],[x[4] for x in objectives]
end
function evaluate(leader,follower)
    Random.seed!(1235)
    n = 200
    T = length(leader.intraday)
    trajectory = []
    vts = leader.vertice
    nk = length(vts[1])
    for i = 1:n
        trajectory_i = [vts[t][randperm(nk)[1]] for t in 1:T]
        push!(trajectory,trajectory_i)
    end
    objectives = []
    for k in 1:n
        m = leader.intraday
        total_cost = 0
        for t in 1:T
            fixprev(leader,t)
            for (idx,gen) in enumerate(keys(leader.case_dict[:windfarm]))
                vertex = trajectory[k][t]
                windpower = leader.data[t][:wind_power]
                fix(m[t][:Pw_err][gen],vertex[idx]*m[t][:α]*windpower[gen])
            end
            Suppressor.@suppress_out optimize!(m[t])
            @assert termination_status(m[t]) == MOI.OPTIMAL #OPTIMAL
            fixprev(follower,value.(m[t][:binder]),t)
            Suppressor.@suppress_out optimize!(follower.intraday[t])
            total_cost += value(m[t][:cost_now]) + value(follower.intraday[t][:cost_now])
        end
        push!(objectives,[total_cost,value(m[T][:Ew])/value(m[T][:Ew_max]),value(m[T][:Ses][1]),sum(value(m[t][:load_Cut_cost]) for t in 1:T)])
    end
    return Tuple([x[i] for x in objectives] for i in 1:length(objectives[1]))
end
function peaksOverThresholdEstimator(data::Array,quantile_α)
    u = quantile(data,quantile_α)
    extreme_values = filter(x->x>=u,data)
    nev = (extreme_values .- minimum(extreme_values))/std(extreme_values)
    N = length(nev)
    function f!(F,x)
        F[1] = 1/x[1] - (1/x[2] + 1) * 1/N * sum(nev[i]/(1 + x[1]*nev[i]) for i in 1:N)
        F[2] = 1/N * sum(log(1 + x[1]*nev[i]) for i in 1:N) - x[2]
    end
    x = NLsolve.nlsolve(f!,[0.1,0.1],autodiff=:forward)
    x = x.zero
    dist = Distributions.GeneralizedPareto(x[2]/x[1],x[2])
    return dist,minimum(extreme_values),std(extreme_values)
end
