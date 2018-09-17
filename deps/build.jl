@info "build.jl is running!"

using Pkg

try
    pkg"rm GroupSlices"
    @info "removed old GroupSlices"
catch
    @info "didn't remov old GroupSlices"
end

# Pkg.add(Pkg.PackageSpec(name="GroupSlices", url="https://github.com/mcabbott/GroupSlices.jl", rev="julia07fixes"))
pkg"add https://github.com/mcabbott/GroupSlices.jl#julia07fixes"

@info "added new GroupSlices"
