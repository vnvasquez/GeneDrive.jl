# GeneDrive.jl

## Overview 

[`GeneDrive.jl`]() is a [Julia](https://julialang.org) package designed for simulating biological dynamics and control. The objectives of the package include: 
* Provide data models that structure inputs to experimental setups and exploit the power of Julia's type system for [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/). 
* Enable the creation of dynamic models that build on the [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) platform.
* Facilitate the formulation of decision models that employ [`JuMP.jl`](https://jump.dev/JuMP.jl/stable/), the domain-specifc modeling language for mathematical optimization embedded in Julia.

## Installation and usage 

`GeneDrive.jl` will work with Julia [version 1.7 and above](https://julialang.org/downloads/). Add the package with:

```julia
julia> ]
(v1.7) pkg> add GeneDrive
```

Begin using the package with: 
```julia
julia> using GeneDrive
```

## Getting started

The [documentation]() features examples as well as more detailed descriptions of package functionalities and applications.  

## License

`GeneDrive.jl` is released under an MIT [License](https://opensource.org/licenses/MIT). This package has been developed in partial fufillment of the requirements for a Master of Science in Electrical Engineering and Computer Science ([EECS](https://eecs.berkeley.edu/research)) at the University of California, Berkeley. 
