using Random
using Distributions

""" This file contains code for providing noisy answers to queries and depends
on queries.jl """


include("queries.jl")

""" The struct that represents a noisy answer. It stores the noisy answer along with
the Distribution used to generate the noise. Distribution is provided from the Distributions
package """
struct NoisyAnswer
   answer  # an array representing the noisy answer
   dist::Distribution # the noise distribution
end

""" Given a random number generator `rng`, a distribution `d`, and an array `a`
add independent noise to each element of `a` and return a NoisyAnswer struct
"""
function addnoise(rng::AbstractRNG, d::Distribution, a::Array)
    NoisyAnswer(a + rand(rng, d, size(a)...), d)
end

""" The laplace mechanism for a single query """
function laplace_mechanism(query::Query; data, epsilon, rng::AbstractRNG)
   sens = l1sens(query)
   lap = Laplace(0, sens/epsilon)
   addnoise(rng, lap, answerquery(query, data))
end

""" The Laplace mechanism for an array of Query subtypes. The privacy budget `epsilon` is split evenly across them"""
function laplace_mechanism(queries::Array{T}; data, epsilon, rng::AbstractRNG) where T <: Query
    totalsens = sum(x->l1sens(x), queries)
    lap = Laplace(0, totalsens/epsilon)
    [addnoise(rng, lap, answerquery(q, data)) for q in queries]
end

""" The zCDP gaussian mechanism for a single Query"""
function gaussian_mechanism(query::Query; data, rho, rng::AbstractRNG)
    sen_squared = l2sens(query)^2
    sigma = sqrt(sen_squared/(2*rho))
    gauss = Normal(0, sigma)
    addnoise(rng, gauss, answerquery(query, data))
end

""" The zCDP gaussian mechanism for an array of Query subtypes. The privacy budget `rho` is split evenly across them"""
function gaussian_mechanism(queries::Array{T}; data, rho, rng::AbstractRNG) where T <: Query
    totalsens_squared = sum(x->l2sens(x)^2, queries)
    sigma = sqrt(totalsens_squared/(2*rho))
    gauss = Normal(0, sigma)
    [addnoise(rng, gauss, answerquery(q, data)) for q in queries]
end
