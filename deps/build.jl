@info "build.jl is running"

using Pkg

Pkg.rm("GroupSlices")

@info "removed"

Pkg.add(PackageSpec(name = "https://github.com/mcabbott/GroupSlices.jl", rev = "julia07fixes"))

@info "added!" 
