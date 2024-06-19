if !any(name -> isdefined(Main, name), [:Makie, :GLMakie, :CairoMakie])
    using CairoMakie
end
import Downloads

# decide theme:
ENV["COLORSCHEME"] = "JuliaDynamics" # or others, see `plottheme.jl`
ENV["BGCOLOR"] = :white       # anything for `backgroundcolor` of Makie
ENV["AXISCOLOR"] = :black           # color of all axis elements (labels, spines, ticks)

try
    Downloads.download(
        "https://raw.githubusercontent.com/Datseris/plottheme/main/plottheme.jl",
        joinpath(@__DIR__, "_plottheme.jl")
    )
    cp(joinpath(@__DIR__, "_plottheme.jl"), joinpath(@__DIR__, "plottheme.jl"); force = true)
    rm(joinpath(@__DIR__, "_plottheme.jl"); force = true)
catch
end

include("plottheme.jl")


"""
    record_interaction(file, figure; framerate = 30, total_time = 10)

Start recording whatever interaction is happening on some `figure` into a video
output in `file` (recommended to end in `".mp4"`).

## Keyword arguments

* `framerate = 30`
* `total_time = 10`: Time to record for, in seconds
* `sleep_time = 1`: Time to call `sleep()` before starting to save.
"""
function record_interaction(file, figure;
        framerate = 30, total_time = 10, sleep_time = 1,
    )
    ispath(dirname(file)) || mkpath(dirname(file))
    sleep(sleep_time)
    framen = framerate*total_time
    record(figure, file; framerate) do io
        for _ in 1:framen
            sleep(1/framerate)
            recordframe!(io)
        end
    end
    return
end
record_interaction(figure::Figure, file; kwargs...) =
record_interaction(file, figure; kwargs...)