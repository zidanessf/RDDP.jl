# function print_helper(f, io, args...)
#     f(stdout, args...)
#     f(io, args...)
# end
import Printf
function print_banner(io)
    println(
        io,
        "--------------------------------------------------------------------------------",
    )
    printstyled(io, "                  DEV   fast Robust Dual Dynamic Programming (fRDDP) (c) Yunhui Shi, 2021\n";color=:red)
    println(io)
end
function print_banner_v2(io)
    println(
        io,
        "--------------------------------------------------------------------------------",
    )
    printstyled(io, "                      D_RDDP_PS.jl (c) Shi Yunhui, 2020\n";color=:red)
    println(io)
end
function print_iteration_header(io)
    println(
        io,
        " Iteration    UpperBound       LowerBound        Gap       Time (s)       # Solves",
    )
end
function print_iteration_header_v2(io)
    println(
        io,
        " Iteration    UBLeader       LBLeader       SVFLeader        TVFFollower        Time (s)       # Solves",
    )
end
print_value(x::Real) = lpad(Printf.@sprintf("%1.6e", x), 13)
print_value(x::Int) = Printf.@sprintf("%9d", x)

function print_iteration(io, additional::Dict)
    print(io, print_value(additional[:Iteration]))
    print(io, "   ", print_value(additional[:UpperBound]))
    print(io, "  ", print_value(additional[:LowerBound]))
    print(io, "  ", print_value(additional[:Gap]))
    print(io, "  ", print_value(additional[:Time]))
    print(io, "  ", print_value(additional[:TotalSolves]))
    println(io)
end
function print_iteration_v2(io, additional::Dict)
    print(io, print_value(additional[:Iteration]))
    print(io, "   ", print_value(additional[:UpperBound]))
    print(io, "  ", print_value(additional[:LowerBound]))
    print(io, "  ", print_value(additional[:SVFLeader]))
    print(io, "  ", print_value(additional[:FollowerValue]))
    print(io, "  ", print_value(additional[:Time]))
    print(io, "  ", print_value(additional[:TotalSolves]))
    println(io)
end