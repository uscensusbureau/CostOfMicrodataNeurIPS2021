using LinearAlgebra

""" This file contains code for manipulating queries. It introduces an
abstract query type Query that should be subtyped for different query classes.
Each query class needs to have 3 functions defined:

1) answerquery(query, data) provides the query answer as an array (it must be an array even
    if the answer has size 1). Data is an array (not a Dataset struct).
    IMPORTANT: this function must work with arrays of numbers along
    with arrays of JuMP variables, since this is used to set constraints in the optimizer using
    the JuMP framework. In particaular, this means that the answerquery function must use simple
    linear algebra (e.g., normal arithmetic operations and indexing, no loops or conditionals)
2)  l1sens(query) returns the L1 sensitivity
3)  l2sens(query) returns the L2 sensitivity
"""


""" The supertype for queries """
abstract type Query end


###############################
""" The Identity query returns the data array, it has l1 and l2 sensitivity of 1"""
struct IdentityQuery <: Query end
answerquery(q::IdentityQuery, data) = data
l1sens(q::IdentityQuery) = 1
l2sens(q::IdentityQuery) = 1

##############################

""" This query represents the sum of all elements in the data array. It has L1 and L2 sensitivity 1"""
struct SumQuery <: Query end
answerquery(q::SumQuery, data::AbstractArray) = [sum(data)]
l1sens(q::SumQuery) = 1
l2sens(q::SumQuery) = 1

##############################

""" The marginal query is constructed as MarginalQuery(numdims, keepdims)
where numdims is the number of dimensions in the dataarray and keepdims are
the dimensions you want to keep. It has L1 and L2 sensitivity of 1.
"""
struct MarginalQuery <: Query
    numdims
    removedims
    MarginalQuery(;numdims, keep) = new(numdims, Tuple(x for x in 1:numdims if !(x in keep)))
end
answerquery(q::MarginalQuery, data::AbstractArray) = dropdims(sum(data, dims=q.removedims), dims=q.removedims)
l1sens(q::MarginalQuery) = 1
l2sens(q::MarginalQuery) = 1

#############################

""" This is a meta query that applies a mask to an existing query.
Given a base query q with answer ans, the mask query takes a mask (array with
same dimensions as ans) and returns the inner product between ans and the mask.
This query is mainly used for the ReWeight function discussed in the paper to
represent the sum of the low queries. Since it is not intended to be used with fresh
noise, the l1sens and l2sens functions provide a loose upper bound on the sensitivity.
 """

struct MaskQuery <: Query
   basequery
   mask
end
answerquery(q::MaskQuery, data::AbstractArray) = [sum(answerquery(q.basequery, data) .* q.mask)]
l1sens(q::MaskQuery) = maximum(abs.(q.mask)) * l1sens(q.basequery)
l2sens(q::MaskQuery) = maximum(abs.(q.mask)) * l2sens(q.basequery)

##########################
""" This is a meta query that zeroes out answers to a base query
Given a base query q with answer a, the filter query takes a mask (bit array)
with same dimension as a and returns a version of array a where elements with 0 mask
are set to 0.
"""
struct FilterQuery <: Query
   basequery
   mask::BitArray
end
answerquery(q::FilterQuery, data::AbstractArray) = answerquery(q.basequery, data) .* q.mask
l1sens(q::FilterQuery) = l1sens(q.basequery)
l2sens(q::FilterQuery) = l2sens(q.basequery)
