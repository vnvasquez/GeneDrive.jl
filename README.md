# GeneDrive.jl

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Main - CI](https://github.com/vnvasquez/GeneDrive.jl/actions/workflows/main-tests.yml/badge.svg)](https://github.com/vnvasquez/GeneDrive.jl/actions/workflows/main-tests.yml)
[![codecov](https://codecov.io/gh/vnvasquez/GeneDrive.jl/branch/main/graph/badge.svg?token=A1C8HACSIP)](https://codecov.io/gh/vnvasquez/GeneDrive.jl)
[![Stable documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://vnvasquez.github.io/GeneDrive.jl/dev/)

## Overview 

[`GeneDrive.jl`](https://vnvasquez.github.io/GeneDrive.jl/dev/) is a [Julia](https://julialang.org) package designed for simulating biological dynamics and control. The objectives of the package include: 
* Provide data models that structure inputs to experimental setups and exploit the power of Julia's type system for [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/). 
* Enable the creation of dynamic models that build on the [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) platform.
* Facilitate the formulation of decision models that employ [`JuMP.jl`](https://jump.dev/JuMP.jl/stable/), the domain-specific modeling language for mathematical optimization embedded in Julia.

## Installation and usage 

`GeneDrive.jl` will work with Julia [version 1.8 and above](https://julialang.org/downloads/). Add the package with:

```julia
julia> 
(v1.9) pkg> add GeneDrive
```

Begin using the package with: 
```julia
julia> using GeneDrive
```

## Getting started

The [documentation](https://vnvasquez.github.io/GeneDrive.jl/dev/) features examples as well as more detailed descriptions of package functionalities.  

## License

`GeneDrive.jl` is released under an MIT [License](https://opensource.org/licenses/MIT). This package has been developed in partial fulfillment of the requirements for a Master of Science in Electrical Engineering and Computer Science ([EECS](https://eecs.berkeley.edu/research)) at the University of California, Berkeley. 
