module MyTest

import Nonneg
using Test
using Random
using Distributions

const SKIP = false
include("test_queries.jl")
include("test_mech.jl")
include("test_fit.jl")

end # Module
