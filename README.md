# Nonlinear Dynamics: A concise introduction interlaced with code

This repository holds the material related with the textbook _Nonlinear Dynamics: A concise introduction interlaced with code_, co-authored by George Datseris and Ulrich Parlitz. This textbook will be published by Springer, in the series Undergraduate Lecture Notes in Physics, in 2021.

**Contents**
1. [Sample material](#sample-material)
2. [Tutorials for Julia and related packages](#tutorials-for-julia-and-related-packages)
3. [Reproducing figures](#reproducing-figures)
4. [Exercise datasets](#exercise-datasets)
5. [Multiple choice questions](#multiple-choice-questions)
6. [Interactive applications & videos](#interactive-applications--videos)
7. [Contributing](#contributing)

## Sample material
Sample material *(currently preliminary)* of the book is contained as `.pdf` in the `sample` folder.

## Tutorials for Julia and related packages
Below we provide links to various sources for learning Julia, or the packages that we use in the code snippets in the book.

- https://www.youtube.com/watch?v=Fi7Pf2NveH0 : Short-timed, intensive Julia workshop, that will teach everything necessary to use Julia. Aimed at people already familiar with programming
- https://juliadynamics.github.io/DynamicalSystems.jl/dev/ : Documentation of DynamicalSystems.jl
- https://juliadynamics.github.io/JuliaDynamics/ : Website of the JuliaDynamics organization
- https://www.youtube.com/watch?v=KPEqYtEd-zY : Introduction to solving differential equations in Julia, which is generally useful for nonlinear dynamics
- https://julialang.org/community/ : Resources for asking questions about Julia
- https://discourse.julialang.org/ : Official Julia forum (and also the main platform that newcomers ask questions)
- https://github.com/rveltz/BifurcationKit.jl : Julia software for bifurcation analysis

## Reproducing figures
The accompanying code base used here is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible environment that creates the figures of the book.
The code that creates the figures is in the `figure_generation` folder. Notice however that some figures were made with (or enhanced by) PowerPoint and thus we do not share this here.

To (locally) reproduce a figure, first install Julia and then do the following:

0. Download this repository as-is and export it to some folder.
0. Install `DrWatson` on your general Julia installation by doing:
   ```
   julia> using Pkg; Pkg.add("DrWatson")
   ```
1. Then do:
   ```
   julia> Pkg.activate("path/to/the/downloaded/project/folder")
   julia> Pkg.instantiate() # installs all packages used by the repo
   ```
1. This repository uses the Python library matplotlib (package `PyPlot`) for book plots and Makie.jl for the interactive applications. You could use any plotting package instead, but if you want to use `PyPlot` and replicate exactly the book plots, then run the following commands to ensure a working installation for all operating systems:
   ```
   julia> ENV["PYTHON"] = ""
   julia> Pkg.add("PyCall"); Pkg.build("PyCall")
   julia> Pkg.add("PyPlot"); using PyPlot
   ```

Now all necessary packages are installed and all scripts should run out of the box.
As you will notice, all scripts start with the commands:
```julia
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
```
which ensures that only local directories will be used, as well as the *exact* package versions contained within the repository, leading to full reproducibility.

## Exercise datasets
The datasets that are used in the book exercises are contained in the `exercise_data` folder, all being in the same text-based format. To load the exercise data you only have to do:
```julia
using DelimitedFiles
X = readdlm("path/to/exercise_data/dataset.csv")
```

The same folder contains information of where this data is coming from: `data_explanations.md`. Some data are generated from simulations in the script `generating_code.jl`.

## Multiple choice questions
Multiple choice questions that we use during lecturing to increase student involvement are in the `multiple_choice` folder.

## Interactive applications & videos
In the folder `animations` we provide scripts that launch interactive applications, and also pre-recorded `.mp4` files for convenience.

## Contributing
Do you have suggestions for improving the contents? Perhaps you want to contribute something new, e.g. new exercises or multiple choice questions? Or perhaps you have a suggestion/comment on a book section? Or maybe you found a typo in the book?

Please do open an issue on this GitHub repository, and we'll take it into consideration for the next version!
