@info "build.jl is running"

using Pkg

try
    Pkg.rm("GroupSlices")
    @info "removed GroupSlices"
catch
    @info "did not remove GroupSlices"
end

Pkg.add(PackageSpec(name="GroupSlices", url="https://github.com/mcabbott/GroupSlices.jl", rev="julia07fixes"))

@info "added GroupSlices"
