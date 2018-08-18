@info "build.jl is running"

using Pkg

Pkg.rm("GroupSlices")

@info "removed GroupSlices"

Pkg.add(PackageSpec(name="GroupSlices", url="https://github.com/mcabbott/GroupSlices.jl", rev="julia07fixes"))

@info "added GroupSlices"
