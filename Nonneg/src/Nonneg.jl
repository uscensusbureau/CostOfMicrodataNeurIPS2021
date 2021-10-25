module Nonneg
""" This is the main entry point into the Nonneg package.

The relevant files are:
- fit.jl, which depends on mechanisms.jl and queries.jl. This is a collection
         of postprocessing routines that take noisy measurements and constructed
         positively weighted datasets with the help of optimizers
- mechanisms.jl, which provides privacy preserving query answers using the Laplace
        or Guassian mechanisms. This is a research-only implementation that is not
        hardened against floating point attacks or attacks on the random number generator.
- queries.jl, which specifies how queries are answered from data.
- data.jl, which loads benchmark datasets
"""

include("fit.jl")
include("data.jl")
end # module
