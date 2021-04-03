# Julia Support for SWMF

We provide two Julia packages for processing SWMF data:
  * [Batsrus.jl](https://github.com/henry2004y/Batsrus.jl) for reading and converting simulation data;
  * [VisAnaJulia](https://github.com/henry2004y/VisAnaJulia) for visualizing and analyzing simulation data.

## Prerequisites

Julia 1.6+

## Installation

We provide a Perl script [gitclone](../Scripts/gitclone) for installation:
```
gitclone Batsrus.jl
gitclone VisAnaJulia
```

Alternatively, please follow the instructions in the online documents. Batsrus.jl is a registered Julia package, which can be installed directly in Julia via
```
using Pkg; Pkg.add("Batsrus")
```

## Usage

For more details, please check the [document for Batsrus.jl](https://henry2004y.github.io/Batsrus.jl/dev/) and [document for VisAnaJulia](https://henry2004y.github.io/VisAnaJulia/dev/).

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## Acknowledgments

* All the nice guys who share their codes!