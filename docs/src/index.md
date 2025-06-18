```@meta
CurrentModule = GeneDrive
```
## Overview

[`GeneDrive.jl`](https://github.com/vnvasquez/GeneDrive.jl) is a [Julia](https://julialang.org) package designed for simulating biological dynamics and control. The current implementation focuses on genetic-based public health interventions that modify populations of disease vectors, including mosquitoes of the *Aedes* and *Anopheles* genera.

The package furnishes a three-part framework for building and analyzing simulations wherein metapopulations are subject to anthropogenic and environmental change:
* Data models that structure inputs to experimental setups and exploit the power of Julia's type system for [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/).
* Dynamic models that build on the [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) platform.
* Decision models that employ [`JuMP.jl`](https://jump.dev/JuMP.jl/stable/), the domain-specific modeling language for mathematical optimization embedded in Julia.

## Installation and usage

`GeneDrive.jl` will work with Julia [version 1.11 and above](https://julialang.org/downloads/). Add the package with:

```julia
julia> Pkg.add(GeneDrive)
```

Begin using the package with:
```julia
julia> using GeneDrive
```

## What's with the name?

Gene drives include both naturally occurring and synthetic genetic elements. They have been harnessed by scientists for potential use in biological control (e.g.: [public health](https://royalsocietypublishing.org/doi/10.1098/rstb.2019.0803), [agriculture](https://royalsocietypublishing.org/doi/10.1098/rspb.2019.1515), [conservation](https://royalsocietypublishing.org/doi/10.1098/rspb.2019.1606)). "Drive" simply refers to these genetic elements being inherited with a high probability. Therefore, their frequency in a population grows quickly - and at the expense of the wildtype population.

While this package is [not exclusively applicable to gene drive research](@ref environmental_dynamics), its name is a nod to this new technological horizon in biological control.

## Citing GeneDrive.jl

[Insert paper citation when available]
