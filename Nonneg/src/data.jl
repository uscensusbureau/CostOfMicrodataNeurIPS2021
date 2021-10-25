module Data
import DelimitedFiles
using Random
""" This module contains code for loading benchmark datasets. The
get_datasets() function will retrieve the datasets broken down by category
(1-d synthetic, 2-d synthetic, 2-d PUMS)"""


""" This is the dataset type, which represents a histogram over a collection of records """
struct Dataset
   name # the name of the dataset
   contents # an array of counts or weights that represents the dataset
end

""" Returns the directory containing the PUMS datasets, assumed to be data/ relative to this file's path """
datadir() = "$(pkgdir(Data))/data"


""" Returns the base name (without extension or path) of the files containing PUMS histograms """
realdatafiles() = copy([
        "ST_01_PUMA_01301",
        "ST_17_PUMA_03529",
        "ST_24_PUMA_01004",
        "ST_29_PUMA_01901",
        "ST_36_PUMA_04010",
        "ST_08_PUMA_00803",
        "ST_17_PUMA_03531",
        "ST_26_PUMA_02702",
        "ST_32_PUMA_00405",
        "ST_51_PUMA_01301",
        "ST_13_PUMA_04600",
        "ST_19_PUMA_01700",
        "ST_28_PUMA_01100",
        "ST_36_PUMA_03710",
        "ST_51_PUMA_51255"
    ])



""" Reads PUMS file from the given filename (no path or extension)
The format is a CSV file and represents a 2-d table of counts. The first
row and first column are headers. Returns a Dataset struct.
"""
function readipums(filename)
    the_data = DelimitedFiles.readdlm("$(datadir())/$(filename).csv", ',', Int, skipstart=1, skipblanks=true)[:, 2:end]
    Dataset(replace(filename, "_"=>"-"), the_data)
end

### Toy Data

""" Returns an array of shape `thesize` whose entries are all equa to `val` """
level_data(thesize, val) = fill(val, thesize)

""" Returns an array of shape `thesize` whose entries sequentially increase starting from 1"""
stair_data(thesize) = reshape([x for x in 1:prod(thesize)], thesize)

""" Returns an array of shape `thesize` whose first half contains 0 and last half contains `val`"""
step_data(thesize, val) = reshape([(if x <= prod(thesize)/2  0 else val end) for x in 1:prod(thesize)], thesize)

""" Returns a dataset of shape `thesize` whose first half of elements sequentially increase from 1, while the last half are 0"""
split_stairs(thesize) = reshape([(if x >= prod(thesize)/2 0 else x end) for x in 1:prod(thesize)], thesize)


""" Takes a base dataset `base_data`, sets the first element to `bigval` and optionally
permutes the data and and optionally adds noise (Laplace with scale 1). `seed` can be specified for
the random operations. `name` will be the name of the dataset and `doshuffle` indicates whether
the data array should be randomly permuted. The  result is a dataset struct.
"""
function make_synth_data(base_data, name, bigval; doshuffle=false, addnoise=0, seed=1)
    the_data = copy(base_data)
    the_data[1] = bigval
    the_size = size(the_data)
    rng = MersenneTwister(seed)
    if doshuffle
        shuffle!(rng, the_data)
    end
    if addnoise > 0
        noisy = the_data + addnoise * (log.(rand(rng, the_size...)) - log.(rand(rng, the_size...)))
        the_data = round.(Int, max.(0, noisy))
    end
    Dataset(
          "$(name):$(bigval):$(doshuffle):$(addnoise):$(seed)",
          the_data
        )
end

""" Returns a tuple of 3 arrays. The first is an array of 1-d synthetic datasets,
the 2nd is an array of 2-d synthetic datasets, and the 3rd is an array of pums
histograms """
function get_datasets()
    onesize = 100
    bigval = 10000
    one_d_synth = [
        make_synth_data(level_data(onesize, 0), "Level00_1d", bigval),
        make_synth_data(level_data(onesize, 1), "Level01_1d", bigval),
        make_synth_data(level_data(onesize, 16), "Level16_1d", bigval),
        make_synth_data(level_data(onesize, 32), "Level32_1d", bigval),

        make_synth_data(stair_data(onesize), "Stair_1d", bigval),

        make_synth_data(step_data(onesize, 16), "Step16_1d", bigval),
        make_synth_data(step_data(onesize, 50), "Step50_1d", bigval),

        make_synth_data(split_stairs(onesize), "SplitStairs_1d", bigval)
    ]
    twosize = (10, 10)
    two_d_synth = [
        make_synth_data(level_data(twosize, 0), "Level00_2d", bigval),
        make_synth_data(level_data(twosize, 1), "Level01_2d", bigval),
        make_synth_data(level_data(twosize, 16), "Level16_2d", bigval),
        make_synth_data(level_data(twosize, 32), "Level32_2d", bigval),

        make_synth_data(stair_data(twosize), "Stair_2d", bigval),

        make_synth_data(step_data(twosize, 16), "Step16_2d", bigval),
        make_synth_data(step_data(twosize, 50), "Step50_2d", bigval),

        make_synth_data(split_stairs(twosize), "SplitStairs_2d", bigval)

    ]
    real_ipums = [
        readipums(f) for f in realdatafiles()
    ]
    return (one_d_synth, two_d_synth, real_ipums)
end

end ### MODULE
