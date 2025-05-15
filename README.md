# VariogramVerify.jl

NOTE: this repo is under active development...

`VariogramVerify.jl` is a Julia package built on `GeoStats.jl` for systematically computing variograms at multiple angles in search for the most appropriate major and minor axes directions for geostatistical interpolation.  This allows user input at the critical interpolation stage of variogram modeling and selection.

One intended use case is in  [Stanford-Mineral-X/MultiGrid.jl](https://github.com/Stanford-Mineral-X/MultiGrid.jl).

## Installation
```
using VariogramVerify
```

## Example
<img src="assets/VariogramVerify_ex.png" alt="Example" width="480px" />
