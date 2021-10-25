


""" Make sure that the directory exists to store experiments for this experiment and dataset combo.
Returns the name of the directory """
function ensure_dir(expid, dataname)
   expdir = "$(RESULTS_DIR)/$(expid)"
   expdatadir = "$(RESULTS_DIR)/$(expid)/$(dataname)"
   if (! isdir(expdir)) mkdir(expdir) end
   if (! isdir(expdatadir)) mkdir(expdatadir) end
   return expdatadir
end


""" given the name of the experiment, data, the name of the algorithm, and an array of its outputs,
serializes this array in the path for this experiment and dataset.
     """
function save_experiment(expname, dataname, methodname, expresults)
   expdatadir = ensure_dir(expname, dataname)
   filename = "$(expdatadir)/$(methodname).$(SER_EXTENSION)"
   #println(expresults)
   serialize(filename, expresults)
end

""" serializes the dataset in the path for this experiment and dataset """
function save_data(expname, data::Nonneg.Data.Dataset)
   dataname = data.name
   expdatadir = ensure_dir(expname, dataname)
   filename = "$(expdatadir)/$(DATA).$(SER_EXTENSION)"
   if !isfile(filename) serialize(filename, data.contents) end
end

""" Saves the text description of the experiment and dataset """
function save_description(expname, dataname, description)
   expdatadir = ensure_dir(expname, dataname)
   filename = "$(expdatadir)/$(DESCRIPTION)"
   if !isfile(filename)
       open(filename, "w") do io
           println(io, description)
       end
   end
end

############################
# Experiments
###########################
""" Returns query workload for 1-d datasets """
function oned_queries()
   [Nonneg.SumQuery(), Nonneg.IdentityQuery()]
end

""" Returns query workload for 2-d datasets """
function twod_queries()
   numdims = 2
   [ Nonneg.SumQuery(),
     Nonneg.MarginalQuery(numdims=numdims, keep=1),
     Nonneg.MarginalQuery(numdims=numdims, keep=2),
     Nonneg.IdentityQuery()
    ]
end

""" First experiment, oned datasets, total epsilon=1 (Laplace mechanism)"""
function experiment1(oned, twod, pums, algs; iterations)
   expname = "experiment1" # name of this function
   epsilon=1
   description = "1-d Datasets. Sum and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in oned # only look at 1-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = oned_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" Second experiment, oned datasets, total epsilon=0.5 (Laplace mechanism)"""
function experiment2(oned, twod, pums, algs; iterations)
   expname = "experiment2" # name of this function
   epsilon=0.5
   description = "1-d Datasets. Sum and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in oned # only look at 1-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = oned_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" Third experiment, oned datasets, total epsilon=0.1 (Laplace mechanism)"""
function experiment3(oned, twod, pums, algs; iterations)
   expname = "experiment3" # name of this function
   epsilon=0.1
   description = "1-d Datasets. Sum and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in oned # only look at 1-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = oned_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end



""" Fourth experiment, twod datasets, total epsilon=1 (Laplace mechanism)"""
function experiment4(oned, twod, pums, algs; iterations)
   expname = "experiment4" # name of this function
   epsilon=1
   description = "2-d Datasets. Sum, marginal and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in twod # only look at 2-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = twod_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" Fifth experiment, twod datasets, total epsilon=0.5 (Laplace mechanism)"""
function experiment5(oned, twod, pums, algs; iterations)
   expname = "experiment5" # name of this function
   epsilon=0.5
   description = "2-d Datasets. Sum, marginal and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in twod # only look at 2-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = twod_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" Sixth experiment, twod datasets, total epsilon=0.1 (Laplace mechanism)"""
function experiment6(oned, twod, pums, algs; iterations)
   expname = "experiment6" # name of this function
   epsilon=0.1
   description = "2-d Datasets. Sum, marginal and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in twod # only look at 2-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = twod_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" 7th experiment, pums datasets, total epsilon=1 (Laplace mechanism)"""
function experiment7(oned, twod, pums, algs; iterations)
   expname = "experiment7" # name of this function
   epsilon=1
   description = "Pums Datasets. Sum, marginal and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in pums # only look at 2-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = twod_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" 8th experiment, pums datasets, total epsilon=0.5 (Laplace mechanism)"""
function experiment8(oned, twod, pums, algs; iterations)
   expname = "experiment8" # name of this function
   epsilon=0.5
   description = "Pums Datasets. Sum, marginal and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in pums # only look at 2-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = twod_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" 9th experiment, pums datasets, total epsilon=0.1 (Laplace mechanism)"""
function experiment9(oned, twod, pums, algs; iterations)
   expname = "experiment9" # name of this function
   epsilon=0.1
   description = "pums Datasets. Sum, marginal and point queries with Laplace mechanism, epsilon=$epsilon"
   for data in pums # only look at 2-d datasets
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           queries = twod_queries()
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.laplace_mechanism(queries, data=data.contents, epsilon=epsilon, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

function gauss_experiment(data_collection, algs; iterations, expname, description, rho, queries)
   for data in data_collection
       save_data(expname, data)
       save_description(expname, data.name, description)
       for myalg in algs
           results = map(1:iterations) do i
               rng = MersenneTwister(i)
               sanitized = Nonneg.gaussian_mechanism(queries, data=data.contents, rho=rho, rng=rng)
               (solution, status) = myalg(size(data.contents), queries, sanitized)
               if (i % 100 == 0) println("$(data.name), $(myalg), Iteration $i") end
               (solution, Nonneg.status_to_string(status))
           end
           save_experiment(expname, data.name,string(Symbol(myalg)),results)
       end
   end
end

""" 10th experiment, 1-d datasets, total rho=1/2 (Gaussian mechanism)"""
function experiment10(oned, twod, pums, algs; iterations)
    data_collection=oned
    expname="experiment10"
    rho=0.5
    description="1-d Datasets. Sum and point queries with Gaussian mechanism, rho=$rho"
    queries=oned_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end
""" 11th experiment, 1-d datasets, total rho=0.125 (Gaussian mechanism)"""
function experiment11(oned, twod, pums, algs; iterations)
    data_collection=oned
    expname="experiment11"
    rho=0.125
    description="1-d Datasets. Sum and point queries with Gaussian mechanism, rho=$rho"
    queries=oned_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end

""" 12th experiment, 1-d datasets, total rho=0.005 (Gaussian mechanism)"""
function experiment12(oned, twod, pums, algs; iterations)
    data_collection=oned
    expname="experiment12"
    rho=0.005
    description="1-d Datasets. Sum and point queries with Gaussian mechanism, rho=$rho"
    queries=oned_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end

""" 13th experiment,  2-d datasets, total rho=1/2 (Gaussian mechanism)"""
function experiment13(oned, twod, pums, algs; iterations)
    data_collection=twod
    expname="experiment13"
    rho=0.5
    description="2-d Datasets. Sum, marginal and point queries with Gaussian mechanism, rho=$rho"
    queries=twod_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end
""" 14th experiment, 2-d datasets, total rho=0.125 (Gaussian mechanism)"""
function experiment14(oned, twod, pums, algs; iterations)
    data_collection=twod
    expname="experiment14"
    rho=0.125
    description="2-d Datasets. Sum, maginal and point queries with Gaussian mechanism, rho=$rho"
    queries=twod_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end

""" 15th experiment, 2-d datasets, total rho=0.005 (Gaussian mechanism)"""
function experiment15(oned, twod, pums, algs; iterations)
    data_collection=twod
    expname="experiment15"
    rho=0.005
    description="2-d Datasets. Sum, marginal and point queries with Gaussian mechanism, rho=$rho"
    queries=twod_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end

""" 16th experiment,  pums datasets, total rho=1/2 (Gaussian mechanism)"""
function experiment16(oned, twod, pums, algs; iterations)
    data_collection=pums
    expname="experiment16"
    rho=0.5
    description="PUMS Datasets. Sum, marginal and point queries with Gaussian mechanism, rho=$rho"
    queries=twod_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end
""" 17th experiment, pums datasets, total rho=0.125 (Gaussian mechanism)"""
function experiment17(oned, twod, pums, algs; iterations)
    data_collection=pums
    expname="experiment17"
    rho=0.125
    description="PUMS Datasets. Sum, maginal and point queries with Gaussian mechanism, rho=$rho"
    queries=twod_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end

""" 18th experiment, pums datasets, total rho=0.005 (Gaussian mechanism)"""
function experiment18(oned, twod, pums, algs; iterations)
    data_collection=pums
    expname="experiment18"
    rho=0.005
    description="PUMS Datasets. Sum, marginal and point queries with Gaussian mechanism, rho=$rho"
    queries=twod_queries()
    gauss_experiment(data_collection, algs, iterations=iterations, expname=expname, description=description, rho=rho, queries=queries)
end


###########################
# End Experiments
###########################
""" Main function for running all experiments """
function experiments_go(iterations=1_000)
    algs = [
             Nonneg.olsalg,
             Nonneg.nnlsalg,
             Nonneg.maxalg,
             Nonneg.seqalg,
	         Nonneg.weightalg
             ]
    oned, twod, pums = Nonneg.Data.get_datasets()
    println("Exp1")
    experiment1(oned, twod, pums, algs, iterations=iterations)
    println("Exp2")
    experiment2(oned, twod, pums, algs, iterations=iterations)
    println("Exp3")
    experiment3(oned, twod, pums, algs, iterations=iterations)
    println("Exp4")
    experiment4(oned, twod, pums, algs, iterations=iterations)
    println("Exp5")
    experiment5(oned, twod, pums, algs, iterations=iterations)
    println("Exp6")
    experiment6(oned, twod, pums, algs, iterations=iterations)
    println("Exp7")
    experiment7(oned, twod, pums, algs, iterations=iterations)
    println("Exp8")
    experiment8(oned, twod, pums, algs, iterations=iterations)
    println("Exp9")
    experiment9(oned, twod, pums, algs, iterations=iterations)
    println("Exp10")
    experiment10(oned, twod, pums, algs, iterations=iterations)
    println("Exp11")
    experiment11(oned, twod, pums, algs, iterations=iterations)
    println("Exp12")
    experiment12(oned, twod, pums, algs, iterations=iterations)
    println("Exp13")
    experiment13(oned, twod, pums, algs, iterations=iterations)
    println("Exp14")
    experiment14(oned, twod, pums, algs, iterations=iterations)
    println("Exp15")
    experiment15(oned, twod, pums, algs, iterations=iterations)
    println("Exp16")
    experiment16(oned, twod, pums, algs, iterations=iterations)
    println("Exp17")
    experiment17(oned, twod, pums, algs, iterations=iterations)
    println("Exp18")
    experiment18(oned, twod, pums, algs, iterations=iterations)
end
