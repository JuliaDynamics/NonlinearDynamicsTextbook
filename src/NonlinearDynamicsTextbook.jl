module NonlinearDynamicsTextbook

using Reexport
@reexport using DynamicalSystems, PyPlot, Random, Statistics

function set_plotting_style!()
    include(joinpath(@__DIR__, "style.jl"))
end

export set_plotting_style!

end