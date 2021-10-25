const INCLUDE_STD = true # include std error into the tables?


""" Maps experiment name to queries that were used and their readable names"""
const ExpMap = Dict(
    "experiment1" => (oned_queries(), ["Sum", "Id"]),
    "experiment2" => (oned_queries(), ["Sum", "Id"]),
    "experiment3" => (oned_queries(), ["Sum", "Id"]),
    "experiment4" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment5" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment6" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment7" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment8" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment9" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment10" => (oned_queries(), ["Sum", "Id"]),
    "experiment11" => (oned_queries(), ["Sum", "Id"]),
    "experiment12" => (oned_queries(), ["Sum", "Id"]),
    "experiment13" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment14" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment15" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment16" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment17" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"]),
    "experiment18" => (twod_queries(), ["Sum", "Marg1", "Marg2", "Id"])
)

""" maps postprocessing function name to readable name for tables """
const AlgsDict = Dict("olsalg"=>"ols",
                      "nnlsalg" => "nnls",
                      "maxalg" => "maxfit",
                      "seqalg" => "seq",
                      "weightalg" => "reweight"
                      )

""" get list of experiment names"""
function get_experiments()
   unsorted = [x for x in readdir(RESULTS_DIR) if startswith(x, "exp")]
   sort(unsorted, by = x->(filter(!isdigit, x), parse(Int32,filter(isdigit, x))))
end

""" given an experiment name, get the names of datasets it uses"""
function get_datanames(expname)
   readdir("$(RESULTS_DIR)/$(expname)")
end

""" given an experiment name and data name, retrieve that dataset """
function get_data(expname,dataname)
   thedata = deserialize("$(RESULTS_DIR)/$(expname)/$(dataname)/$(DATA).$(SER_EXTENSION)")
end

""" given an experiment name and ata name, get its descrption """
function get_description(expname, dataname)
    strip(read("$(RESULTS_DIR)/$(expname)/$(dataname)/$(DESCRIPTION)", String))
end

""" printable version of data name"""
function pp_dataname(dataname)
    if occursin("PUMA", dataname)
        parts = split(dataname, "-")
        "$(parts[3])$(parts[2])$(parts[4])"
    else
        replace(split(dataname, ":")[1], "_"=>"-")
    end
end

""" given an experiment and data name and an alg name, unserialize its results """
function get_result(expname, dataname, algname)
    filename = "$(RESULTS_DIR)/$(expname)/$(dataname)/$(algname).$(SER_EXTENSION)"
    if isfile(filename)
        deserialize(filename)
    else
        #println("'",filename,"'")
        missing
    end
end

""" Print a summary of model optimization codes for a given experiment"""
function analyze_codes(expname)
    datanames = get_datanames(expname)
    status_counts = Dict([alg=>Dict() for alg in keys(AlgsDict)])
    for d in datanames
        for alg in keys(AlgsDict)
            result = get_result(expname, d, alg)
            if !ismissing(result)
                for (reconstructed, statuses) in result
                    for s in Set(if isa(statuses, Array) statuses else [statuses] end)
                        status_counts[alg][s] = 1 + get(status_counts[alg], s, 0)
                    end
                end
            end
        end
    end
    for alg in sort(collect(keys(AlgsDict)))
        println(alg,":")
        for s in sort(collect(keys(status_counts[alg])))
            println("    ", s, ": ", status_counts[alg][s])
        end
    end
end

""" isok is used to check whether the status codes returned in the optimization merit
its inclusion in the experiments. So if you only want to include results where the optimizer
succeeded with optimal codes, okstats would be ["Optimal"] and isok will be checking if the
actual code belongs in the array okstats. If you provide an array of codes instead of a single
code, it will check that each element of this array is in okstats. If you just want to avoid
complete model failure, set okstats=["OPTIMAL", "ITERATION_LIMIT"]"""
isok(code::String, okstats) = in(code, okstats)
isok(codelist::Array, okstats) = all(isok(code, okstats) for code in codelist)

""" For a given experiment and dataset used in that experiemnt, check all of the runs whose optimization
statuses are in okstats and compute summary statistics like the number, the max expected error of
each query type, the std dev of the expected error, the sum of the expected errors, etc. """
function analyze_experiment_part(expname, dataname, alg; okstats=["OPTIMAL", "ITERATION_LIMIT"])
    (queries, querynames) = ExpMap[expname]
    info = Dict([qn=>Dict() for qn in querynames]) # summary statistics about query errors
    tmperr = Dict() # holds error values for each query aggregated across iterations
    tmpSQerr = Dict() # used to estimate std of worst error query
    answerSizes = Dict() # used to store answer size of each query
    thedata = get_data(expname, dataname)
    thedescription = get_description(expname, dataname)
    result_set = get_result(expname, dataname, alg)
    thecount = 0
    if !ismissing(result_set)
        for (reconstructed, statuses) in result_set
            if isok(statuses, okstats)
                for (q,qn) in zip(queries, querynames)
                    recon_answer = Nonneg.answerquery(q,reconstructed)
                    if !in(qn, keys(answerSizes))
                        answerSizes[qn] = length(recon_answer)
                    end
                    therror = (Nonneg.answerquery(q,thedata) - recon_answer) .^ 2
                    tmperr[qn] = therror .+ get(tmperr, qn, 0)
                    tmpSQerr[qn] = therror .^ 2 .+ get(tmpSQerr, qn, 0)
                end
                thecount += 1
            end
        end
    end
    if (!ismissing(result_set) && thecount > 0)
        for qn in  querynames
            info[qn][STAT_COUNT] = thecount
            info[qn][STAT_MAX] = maximum(tmperr[qn]) / thecount
            info[qn][STAT_TOTAL_ERROR] = sum(tmperr[qn]) / thecount
            (q50, q75, q90) = quantile(vec(tmperr[qn]), [0.5, 0.75, 0.9])
            info[qn][STAT_MED] = q50/thecount
            info[qn][STAT_75] = q75/thecount
            info[qn][STAT_90] = q90/thecount
            focus_index = argmax(tmperr[qn]) #had largest expected value
            #estimate variance of the mean statistic (estimated variance over n)
            allvars =  ((tmpSQerr[qn] ./ thecount) .- (tmperr[qn] ./ thecount).^2) ./thecount
            std_of_query = sqrt(allvars[focus_index])
            std_of_total = sqrt(sum(allvars))
            info[qn][STAT_STD_OF_TOTAL] = std_of_total
            info[qn][STAT_STD_OF_QUERY] = std_of_query
        end
        info
    else
        missing
    end
end


""" For a given experiment and algs you want to compare, create a latex table with `caption` included
in the caption. THe table will be created in a new latex file, and 'input(the latex file)' will be
printed to the strea `ioAll`"""
function analyze_experiment(expname, algs, caption, ioAll; okstats=["OPTIMAL", "ITERATION_LIMIT"])
    SUM_QUERY="Sum"
    datanames = get_datanames(expname)
    thedescription = get_description(expname, datanames[1])
    (queries, querynames) = ExpMap[expname]
    for qn in querynames
        if qn == SUM_QUERY
            continue
        end
        filenameA = "$(CHARTS_DIR)/$(expname)-$(qn)-basic.tex"
        println(ioAll, "\\input{$(expname)-$(qn)-basic.tex}")
        ioA = open(filenameA, "w")
        println(ioA, raw"\begin{table}[H]")
        println(ioA, raw"\begin{tabular}{|l|", repeat("rr|", length(algs)), raw"}")
        #println(ioA, raw"\begin{tabular}{l", repeat("rr", length(algs)), raw"}")
        println(ioA, raw"\multicolumn{1}{c}{}")
        for alg in algs
            println(ioA, " & \\multicolumn{2}{c}{\\textbf{$(alg)}}")
        end
        println(ioA, raw"\\ ") # extra space needed to avoid escaping
        println(ioA, raw"\multicolumn{1}{c|}{\textbf{Dataset}}", repeat("& Total & Max ",length(algs)), raw"\\\hline")
        #println(ioA, raw"\multicolumn{1}{c}{\textbf{Dataset}}", repeat("& Total & Max ",length(algs)), raw"\\\hline")

        filenameB = "$(CHARTS_DIR)/$(expname)-$(qn)-quant.tex"
        #println(ioAll, "\\input{$(expname)-$(qn)-quant.tex}")
        ioB = open(filenameB, "w")
        println(ioB, raw"\begin{table}[H]")
        println(ioB, raw"\begin{tabular}{|l|", repeat("rrr|", length(algs)), raw"}")
        #println(ioB, raw"\begin{tabular}{l", repeat("rrr", length(algs)), raw"}")
        println(ioB, raw"\multicolumn{1}{c}{}")
        for alg in algs
            println(ioB, " & \\multicolumn{3}{c}{\\textbf{$(alg)}}")
        end
        println(ioB, raw"\\ ") # extra space needed to avoid escaping
        println(ioB, raw"\multicolumn{1}{c|}{\textbf{Dataset}}", repeat(raw"& 50\% & 75\% & 90\% ", length(algs)), raw"\\\hline")
        #println(ioB, raw"\multicolumn{1}{c}{\textbf{Dataset}}", repeat(raw"& 50\% & 75\% & 90\% ", length(algs)), raw"\\\hline")

        for d in datanames
            println(ioA, pp_dataname(d))
            stdstring = ""
            println(ioB, pp_dataname(d))
            for alg in algs
                info = analyze_experiment_part(expname, d, alg, okstats=okstats)
                if !ismissing(info)
                    #@printf(ioA, " &\$\\begin{array}{r} %.1f\\\\ \\pm%.1f\\end{array}\$&\$\\begin{array}{r}%.1f\\\\ \\pm%.1f \\end{array}\$ ", info[qn][STAT_TOTAL_ERROR], info[qn][STAT_STD_OF_TOTAL], info[qn][STAT_MAX], info[qn][STAT_STD_OF_QUERY])
                    @printf(ioA, " & %.1f & %.1f ", info[qn][STAT_TOTAL_ERROR], info[qn][STAT_MAX])
                    stdstring = string(stdstring, @sprintf(" &\$\\pm\$%.1f&\$\\pm\$%.1f ", info[qn][STAT_STD_OF_TOTAL], info[qn][STAT_STD_OF_QUERY]))
                    @printf(ioB, " & %.1f & %.1f & %.1f", info[qn][STAT_MED], info[qn][STAT_75], info[qn][STAT_90])
                else
                    print(ioA, " & NA & NA ")
                    stdstring = string(stdstring, "& NA & NA ")
                    print(ioB, " & NA & NA & NA ")
                end
            end
            println(ioA,raw"\\ ")
            if INCLUDE_STD
                println(ioA,"$(stdstring)\\\\")
            end
            println(ioB,raw"\\ ")
        end

        println(ioA, raw"\hline\end{tabular}")
        println(ioA, raw"\caption{Squared Errors (with standard deviations). ", "$(qn) Query. $(caption)",raw"}\label{table:", "$(expname):$(qn):sq",raw"}")
        println(ioA, raw"\end{table}")
        close(ioA)

        println(ioB, raw"\hline\end{tabular}")
        println(ioB, raw"\caption{Squared Error Quantiles. ", "$(qn) Query. $(caption)",raw"}\label{table:", "$(expname):$(qn):quant",raw"}")
        println(ioB, raw"\end{table}")
        close(ioB)
    end
    # now for sum query
    qn = SUM_QUERY
    if in(qn, querynames)
        filenameC = "$(CHARTS_DIR)/$(expname)-$(qn).tex"
        println(ioAll, "\\input{$(expname)-$(qn).tex}")
        ioC = open(filenameC, "w")
        println(ioC, raw"\begin{table}[H]")
        println(ioC, raw"\begin{tabular}{|l|", repeat("r", length(algs)), raw"|}")
        #println(ioC, raw"\begin{tabular}{l", repeat("r", length(algs)), raw"}")
        println(ioC, raw"\multicolumn{1}{c}{\textbf{Dataset}}")
        for alg in algs
            println(ioC, " & \\multicolumn{1}{c}{\\textbf{$(alg)}}")
        end
        println(ioC, raw"\\\hline ") # extra space needed to avoid escaping
        for d in datanames
            println(ioC, pp_dataname(d))
            stdstring=""
            for alg in algs
                info = analyze_experiment_part(expname, d, alg, okstats=okstats)
                if !ismissing(info)
                    #@printf(ioC, " & \$\\begin{array}{r}%.1f\\\\%.1f\\end{array}\$ ", info[qn][STAT_TOTAL_ERROR], info[qn][STAT_STD_OF_TOTAL])
                    @printf(ioC, " & %.1f", info[qn][STAT_TOTAL_ERROR])
                    stdstring = string(stdstring, @sprintf(" & \$\\pm\$%.1f ", info[qn][STAT_STD_OF_TOTAL]))
                else
                    print(ioC, " & NA ")
                    stdstring = string(stdstring, "& NA ")
                end
            end
            println(ioC,raw"\\ ")
            if INCLUDE_STD
                println(ioC,"$(stdstring)\\\\")
            end
        end

        println(ioC, raw"\hline\end{tabular}")
        println(ioC, raw"\caption{Squared Error (with standard deviations). ", "$(qn) Query. $(caption)",raw"}\label{table:", "$(expname):$(qn):sq",raw"}")
        println(ioC, raw"\end{table}")
        close(ioC)

    end
end

""" main function for creating the latex charts"""
function analyze_go()
    ioAll = open("$(CHARTS_DIR)/all.tex", "w")
    println(ioAll, raw"\documentclass{amsart}")
    println(ioAll, raw"\usepackage[margin=1in]{geometry}")
    println(ioAll, raw"\usepackage{float}")
    println(ioAll, raw"\usepackage{multirow}")
    println(ioAll, raw"\begin{document}")
    println(ioAll, raw"\textbf{Complete experimental results. NA indicates that the optimizer was not able to solve the required optimization problems.}")

    analyze_experiment("experiment1",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "1-d datasets. Lap Mechanism (\$\\epsilon=1\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment2",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "1-d datasets. Lap Mechanism (\$\\epsilon=0.5\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment3",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "1-d datasets. Lap Mechanism (\$\\epsilon=0.1\$). ", ioAll, okstats=["OPTIMAL"])

    println(ioAll, raw"\newpage")

    analyze_experiment("experiment4",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "2-d datasets. Lap Mechanism (\$\\epsilon=1\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment5",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "2-d datasets. Lap Mechanism (\$\\epsilon=0.5\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment6",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "2-d datasets. Lap Mechanism (\$\\epsilon=0.1\$). ", ioAll, okstats=["OPTIMAL"])

    println(ioAll, raw"\newpage")

    analyze_experiment("experiment7",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "PUMS datasets. Lap Mechanism (\$\\epsilon=1\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment8",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "PUMS datasets. Lap Mechanism (\$\\epsilon=0.5\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment9",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "PUMS datasets. Lap Mechanism (\$\\epsilon=0.1\$). ", ioAll, okstats=["OPTIMAL"])

    println(ioAll, raw"\newpage")

    analyze_experiment("experiment10",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "1-d datasets. Gauss Mechanism (\$\\rho=0.5\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment11",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "1-d datasets. Gauss Mechanism (\$\\rho=0.125\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment12",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "1-d datasets. Gauss Mechanism (\$\\rho=0.005\$). ", ioAll, okstats=["OPTIMAL"])

    println(ioAll, raw"\newpage")

    analyze_experiment("experiment13",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "2-d datasets. Gauss Mechanism (\$\\rho=0.5\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment14",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "2-d datasets. Gauss Mechanism (\$\\rho=0.125\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment15",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "2-d datasets. Gauss Mechanism (\$\\rho=0.005\$). ", ioAll, okstats=["OPTIMAL"])

    println(ioAll, raw"\newpage")

    analyze_experiment("experiment16",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "PUMS datasets. Gauss Mechanism (\$\\rho=0.5\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment17",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "PUMS datasets. Gauss Mechanism (\$\\rho=0.125\$). ", ioAll, okstats=["OPTIMAL"])
    analyze_experiment("experiment18",  [Nonneg.olsalg, Nonneg.nnlsalg, Nonneg.maxalg, Nonneg.seqalg, Nonneg.weightalg], "PUMS datasets. Gauss Mechanism (\$\\rho=0.005\$). ", ioAll, okstats=["OPTIMAL"])


    println(ioAll, raw"\end{document}")
    close(ioAll)
end
