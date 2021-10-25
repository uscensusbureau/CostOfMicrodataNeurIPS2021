using JuMP
using COSMO
#using OSQP # if you want to use OSSQP
#using Gurobi #if you want to use gurobi

include("mechanisms.jl")

""" converts and optimality status code from JuMP into a string"""
status_to_string(s) = string(s)

""" converts and array of optimality status code from JuMP into an array of strings.
This is useful because different fitting routines will return either a status code or
an array of status codes  """
status_to_string(s::Array) = [status_to_string(x) for x in s]

""" Given a the datashape, returns a tuple of  model and mainvars, where
mainvars has the same shape as datashape and represents the variables to be
optimized over. By default, it uses the COSMO optimizer"""
function getmodel(datashape; nonneg=False)
    #model  = Model(Gurobi.Optimizer) # if you want to use gurobi

    #these are cosmo settings
    model = Model(COSMO.Optimizer)
    set_optimizer_attribute(model, "eps_abs", 1e-7)
    set_optimizer_attribute(model, "eps_rel", 1e-7)
    set_optimizer_attribute(model, "max_iter", 20000)

    # These are OSQP settings
    #model = Model(OSQP.Optimizer)
    #set_optimizer_attribute(model, "eps_abs", 1e-9)
    #set_optimizer_attribute(model, "rel_abs", 1e-9)
    #set_optimizer_attribute(model, "polish", true)
    set_silent(model)
    numvars = reduce(*, datashape)
    mainvars = if nonneg @variable(model,  mainvars[1:numvars], lower_bound=0) else @variable(model, mainvars[1:numvars]) end
    (model, reshape(mainvars, datashape))
end

""" Predicate to check if the model returned an optimal answer """
is_optimal(status) = JuMP.MathOptInterface.OPTIMAL == status

""" Predicate to check if the model stopped early """
stopped_early(status) = JuMP.MathOptInterface.ITERATION_LIMIT == status


""" Adds the constraint query(mainvars) - ans = residuals to the model and returns the new variables residuals """
function addterm(model, query::Query, mainvars, ans)
    numanswers = length(ans)
    tmpvars = @variable(model, [1:numanswers])
    residuals = reshape(tmpvars, size(ans))
    @constraint(model, answerquery(query, mainvars) .- ans .== residuals)
    residuals
end

""" Adds and equality constraint to the model with slack parameters: |query(mainvars) - ans| <= tolerance """
function add_eq_constraint(model, query::Query, mainvars, ans; tolerance=0.001)
    @constraint(model, answerquery(query, mainvars) .- ans .<= tolerance)
    @constraint(model, answerquery(query, mainvars) .- ans .>= -tolerance)
end

""" performs weighted least squares given an array of `Query`s, `NoisyAnswer's, and reweights.
Solves a weighted least squares optimization problem (optionally nonnegative) where the weight is
reweight/variance(query noise)
returns a tuple: the solution and termination status of the model
Note that if a NoisyAnswer has std = 0, then it will be treated as an equality constraint.
The reweights is an array, with one element for each query. This element can either be a constant
(all answers in that query are weighted by that constant) or an array of weights of the same shape
as the query answers would be
"""
function olsfit(datashape, queries::Array{T}, ansarray::Array{S}; reweights=[1 for _ in queries], nonneg=false) where {T <: Query, S <: NoisyAnswer}
    objfunction = 0.0
    model, mainvars = getmodel(datashape, nonneg=nonneg)
    for (q, a, w) in zip(queries, ansarray, reweights)
        if (std(a.dist) == 0)
            add_eq_constraint(model, q, mainvars, a.answer)
        else
            residual = addterm(model, q, mainvars, a.answer)
	        #theweight = w / std(a.dist)^2
            theweight = w ./ std(a.dist)^2
	        #objfunction += sum(r -> theweight*r^2, residual)
            objfunction += sum(theweight .* (residual .^ 2))
        end
    end
    @objective(model, Min, objfunction)
    optimize!(model)
    solutions = value.(mainvars)
    if nonneg
        solutions .= max.(solutions, 0.0)
    end
    (solutions, termination_status(model))
end

""" Finds a closest dataset to the noisy query answers in the weighted Linfintiy metric
we find the smallest x and coresponding histogram D such that std(noise) * x >= |noisy-answer - q(D)|
for all queries. D is not unique. Returns (D, x, termination status)
"""
function nearestfit(datashape, queries::Array{T}, ansarray::Array{S}; nonneg=false) where {T <: Query, S <: NoisyAnswer}
    model, mainvars = getmodel(datashape, nonneg=nonneg)
    ub = @variable(model, ub)
    for (q, a) in zip(queries, ansarray)
        residual = addterm(model, q, mainvars, a.answer)
        @constraint(model, std(a.dist) * ub .>= residual)
        @constraint(model, std(a.dist) * ub .>= -1 .* residual)
    end
    @objective(model, Min, ub^2)
    optimize!(model)
    solutions = value.(mainvars)
    if nonneg #may alter solution and ub value
        solutions .= max.(solutions, 0.0)
    end
    maxresiduals = [maximum(abs.(Nonneg.answerquery(q, solutions) - a.answer))/std(a.dist) for (q,a) in zip(queries, ansarray) if std(a.dist) != 0]
    ubvalue = maximum(maxresiduals)
    #println("mysol ", solutions)
    #println("mysol2 ", value.(mainvars))
    #println("actual ub ", value(ub))
    #println("ub here ", ubvalue)
    if (!is_optimal(termination_status(model))) println("Status: ", termination_status(model)) end
    (solutions, ubvalue, termination_status(model))
end

""" Finds the closest histogram D that minimizes the following error measure
between its answers and the query answers: minimize std weighted Linfinity error
and break ties (secondary objective) using inverse variance weighted least squares.
This is implemented in two stages: findining the minimium distance, and then breaking ties.
Returns the solution and an array of termination statuses. The first status is for finding
the max feasible Linfinity error and the second is for finding the overall solution.
    There is a slack variable added so that the solution can have slightly larger Linfinity distance
    than what was computed in the first stage of the solve. This helps avoid infeasibility errors. """
function maxfit(datashape, queries::Array{T}, ansarray::Array{S}; nonneg=false, slack=0.01) where {T <: Query, S <: NoisyAnswer}
    (_, ub, first_status) = nearestfit(datashape, queries, ansarray, nonneg=nonneg)
    model, mainvars = getmodel(datashape, nonneg=nonneg)
    objfunction = 0.0
    for (q, a) in zip(queries, ansarray)
        residual = addterm(model, q, mainvars, a.answer)
        @constraint(model, std(a.dist) * (ub + slack) .>= residual)
        @constraint(model, std(a.dist) * (ub + slack) .>= -1 .* residual)
        theweight = 1 / std(a.dist)^2
	    objfunction += sum(r -> theweight*r^2, residual)
    end
    @objective(model, Min, objfunction)
    optimize!(model)
    solutions = value.(mainvars)
    if nonneg #may alter solution and ub value
        solutions .= max.(solutions, 0.0)
    end
    #println("ub ", ub)
    if (!is_optimal(termination_status(model))) println("Status2: ", termination_status(model)) end
    #if (termination_status(model) == JuMP.MathOptInterface.NUMERICAL_ERROR) println(objfunction) end

    (solutions, [first_status, termination_status(model)])
end

""" Iteratively solves the ols problem by recursively fitting query j subject to the constraints
that the previous queries must match the answers obtained in the previous solve (fitting query j-1)
returns the solution and list of statuses from intermediate solves.
"""
function sequential_fit(datashape, queries::Array{T}, ansarray::Array{S};slack=0.001) where {T <: Query, S <: NoisyAnswer}
    prev_queries = T[]
    prev_answers = S[]
    statuses = []
    solution = []
    for (q,a) in zip(queries, ansarray)
        current_queries = vcat(prev_queries, [q])
        current_answers = vcat(prev_answers, [a])
        (solution, status) = olsfit(datashape, current_queries, current_answers, nonneg=true)
        if !is_optimal(status) println("Status3: ", status) end
        statuses = vcat(statuses, status)
        new_answer = answerquery(q, solution)
        prev_queries = vcat(prev_queries, [q])
        prev_answers = vcat(prev_answers, [NoisyAnswer(new_answer, Dirac(0))])
    end
    (solution, statuses)
end

""" Returns a mask indicating small noisy query answers, part of the reweight algorithm.
p 1-confidence parameter """
function getmask(ans::NoisyAnswer; p=0.01)
    amount = length(ans.answer)
    dist = ans.dist
    ind_and_x = collect(zip(1:amount, sort(vec(ans.answer))))
    index = findfirst(x -> x[1] * logcdf(dist, x[2]) >= log1p(-p), ind_and_x)
    if index === nothing
        cutoff = ind_and_x[end][2] + 1
    elseif index == 1
        cutoff = ind_and_x[1][2] - 1
    else
        cutoff = ind_and_x[index-1][2]
    end
    ans.answer .<= cutoff

end

""" This is the reweight algorithm from the paper """
function weighted_fit(datashape, queries::Array{T}, ansarray::Array{S}) where {T <: Query, S <: NoisyAnswer}
    #(olsans, status1) = olsfit(datashape, queries, ansarray)
    mod_queries = T[]
    mod_answers = S[]
    mod_weights = []
    for (q, a) in zip(queries, ansarray)
        mask = getmask(a)
        thefilter = 1 .- mask
        if sum(thefilter) != 0 # process the large query answers
            newq = FilterQuery(q, thefilter)
            newa = NoisyAnswer(a.answer .* thefilter, a.dist)
            mod_queries = vcat(mod_queries, newq)
            mod_answers = vcat(mod_answers, newa)
            mod_weights = vcat(mod_weights, [1])
        end
        if sum(mask) != 0 # process the small query answers
            maxmed = max(1, quantile(a.dist, (1/2)^(1/sum(mask))))
            newq1 = FilterQuery(q, mask)
            newa1 = NoisyAnswer(a.answer .* mask, a.dist)
            mod_queries = vcat(mod_queries, newq1)
            mod_answers = vcat(mod_answers, newa1)
            mod_weights = vcat(mod_weights, [1/maxmed^2]/2)

            newq2 = MaskQuery(q, mask)
            newa2 = NoisyAnswer([sum(a.answer .* mask)], Normal(0, sqrt(sum(mask)) * std(a.dist)))
            mod_queries = vcat(mod_queries, newq2)
            mod_answers = vcat(mod_answers, newa2)
            mod_weights = vcat(mod_weights, [1/2])
        end
    end
    olsfit(datashape,mod_queries, mod_answers, reweights=mod_weights, nonneg=true)
end



##### algs ####
# call each one using datahape, queries, ansarray
""" These are the baseline algorithms, all having the same call structure"""
olsalg(datashape, queries, ansarray) = olsfit(datashape, queries, ansarray, nonneg=false)
nnlsalg(datashape, queries, ansarray) = olsfit(datashape, queries, ansarray, nonneg=true)
maxalg(datashape, queries, ansarray) = maxfit(datashape, queries, ansarray, nonneg=true)
seqalg(datashape, queries, ansarray) = sequential_fit(datashape, queries, ansarray)
weightalg(datashape, queries, ansarray) = weighted_fit(datashape, queries, ansarray)
