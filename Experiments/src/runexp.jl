module Experiments
import Pkg
using Serialization
using Random
using Statistics
using Printf

const LOCATION = "$(@__DIR__)"
const MAIN_PACKAGE = "$(LOCATION)/../../Nonneg"
const EXPERIMENTS = "$(LOCATION)/experiments.jl"
const ANALYZE = "$(LOCATION)/analyze.jl"

const RESULTS_DIR = "$(LOCATION)/../results"
const CHARTS_DIR = "$(LOCATION)/../charts"
const SER_EXTENSION = "ser"
const DATA = "data"
const TXT_EXTENSION = "txt"
const DESCRIPTION = "description.txt"

const STAT_COUNT = "Count"
const STAT_MAX = "Max"
const STAT_TOTAL_ERROR = "Total Error"
const STAT_MED = "Median"
const STAT_75 = raw"75\%"
const STAT_90 = raw"95\%"
const STAT_STD_OF_TOTAL = "STD_OF_TOTAL"
const STAT_STD_OF_QUERY = "STD_OF_QUERY"

Pkg.activate(MAIN_PACKAGE)
import Nonneg


include(EXPERIMENTS)
include(ANALYZE)


if abspath(PROGRAM_FILE) == @__FILE__
    #if this file is being executed, instead of imported
    experiments_go() # run the experiments
    analyze_go() # create reports
end


end # MODULE
