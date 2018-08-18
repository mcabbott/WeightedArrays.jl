@info "build.jl is running!"

using Pkg

Pkg.rm("GroupSlices")

@info "removed old GroupSlices"

Pkg.add(Pkg.PackageSpec(name="GroupSlices", url="https://github.com/mcabbott/GroupSlices.jl", rev="julia07fixes"))

@info "added new GroupSlices"
