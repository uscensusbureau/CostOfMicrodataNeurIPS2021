""" given an list of size 2 arrays/tuples, find correlation between first and second element"""
function correlation(twocols)
    col1 = [x[1] for x in twocols]
    col2 = [x[2] for x in twocols]
    n = length(col1)
    mean1 = sum(col1)/n
    mean2 = sum(col2)/n
    covar1 = sqrt(sum(x -> (x-mean1)^2, col1))
    covar2 = sqrt(sum(x -> (x-mean2)^2, col2))
    top = sum(x -> (x[1]-mean1) * (x[2]-mean2), twocols)
    top/((covar1)*(covar2))
end

@testset "analytical_std" begin
    if SKIP return end
    rng = MersenneTwister() # not secure but repeatable if given seed
    a = 2
    b = 3
    data = zeros(a,b)
    for i in 1:a
        for j in 1:b
            data[i,j]=i+j
        end
    end

    q1 = Nonneg.IdentityQuery()
    q2 = Nonneg.SumQuery()
    q3 = Nonneg.MarginalQuery(numdims=2, keep=1)
    q4 = Nonneg.MarginalQuery(numdims=2, keep=[2])

    queries = [q1, q2, q3, q4]
    epsilons = [1.2, 2.3, 3.4, 5.6]
    rhos = [0.8, 0.9, 1.0, 2.2]

    for i in 1:length(queries)
        lap = Nonneg.laplace_mechanism(queries[1:i], data=data, epsilon=epsilons[i], rng=rng)
        g = Nonneg.gaussian_mechanism(queries[1:i], data=data, rho=rhos[i], rng=rng)
        for noisy in lap
            @test std(noisy.dist) ≈ i/epsilons[i] * sqrt(2) atol=1e-7
        end
        for noisy in g
            @test std(noisy.dist) ≈ sqrt(i/(2*rhos[i])) atol=1e-7
        end
    end
end

@testset "emprical cdf" begin
    if SKIP return end
    rng = MersenneTwister() # not secure but repeatable if given seed
    num_iterations=1_000_000
    a = 2
    b = 3
    data = zeros(a,b)
    for i in 1:a
        for j in 1:b
            data[i,j]=i+j
        end
    end

    q1 = Nonneg.IdentityQuery()
    q2 = Nonneg.SumQuery()
    q3 = Nonneg.MarginalQuery(numdims=2, keep=1)
    q4 = Nonneg.MarginalQuery(numdims=2, keep=[2])

    queries = [q1, q2, q3, q4]
    epsilons = [1.2, 2.3, 3.4, 5.6]
    rhos = [0.8, 0.9, 1.0, 2.2]

    row = 1 # row index of first query result we will check
    col1 = 1 # this and next are the two column indexes of query result we will check
    col2 = 2
    bound = 1.2 # check cdf of noise distributions at this point

    atolstd = 4 # tolerate at most 4 std of error in empirical checks

    for i in 1:length(queries)
        lap = [(Nonneg.laplace_mechanism(queries[1:i], data=data, epsilon=epsilons[i], rng=rng))[1].answer[row, [col1, col2]]  for _ in 1:num_iterations]
        g = [(Nonneg.gaussian_mechanism(queries[1:i], data=data, rho=rhos[i], rng=rng))[1].answer[row, [col1, col2]] for _ in 1:num_iterations]
        mean = Nonneg.answerquery(queries[1], data)[row, [col1, col2]]
        frac_lap = sum(x-> (x-mean)[1] <= bound, lap) / num_iterations
        frac_g = sum(x -> (x-mean)[1] <= bound, g) / num_iterations
        expected_lap = cdf(Laplace(0, i/epsilons[i]), bound)
        # check some cdf values for the actual noise
        @test frac_lap ≈ expected_lap atol=atolstd * sqrt(expected_lap * (1-expected_lap)/num_iterations)
        expected_g = cdf(Normal(0, sqrt(i/(2*rhos[i]))), bound)
        @test frac_g ≈ expected_g atol=atolstd * sqrt(expected_g * (1-expected_g)/num_iterations)
        # check that noise in different queries is not correlated
        corr_lap = correlation(lap)
        @test corr_lap ≈ 0 atol = atolstd * sqrt(0.5/num_iterations)
        corr_g = correlation(g)
        @test corr_g ≈ 0 atol = atolstd * sqrt(0.5/num_iterations)
    end
end


@testset "emprical cdf single query" begin
    if SKIP return end
    #if SKIP return end
    rng = MersenneTwister() # not secure but repeatable if given seed
    num_iterations=1_000_000
    a = 2
    b = 3
    data = zeros(a,b)
    for i in 1:a
        for j in 1:b
            data[i,j]=i+j
        end
    end

    q1 = Nonneg.IdentityQuery()
    q2 = Nonneg.SumQuery()
    q3 = Nonneg.MarginalQuery(numdims=2, keep=1)
    q4 = Nonneg.MarginalQuery(numdims=2, keep=[2])

    queries = [q1, q2, q3, q4]
    epsilons = [1.2, 2.3, 3.4, 5.6]
    rhos = [0.8, 0.9, 1.0, 2.2]

    row = 1 # row index of first query result we will check
    col = 1 # column index of query result we will check
    bound = 1.2 # check cdf of noise distributions at this point
    atolstd = 4 # tolerate at most 4 std of error in empirical checks

    for (q,e,r) in zip(queries, epsilons, rhos)
        lap = [(Nonneg.laplace_mechanism(q, data=data, epsilon=e, rng=rng)).answer[row, col]  for _ in 1:num_iterations]
        g = [(Nonneg.gaussian_mechanism(q, data=data, rho=r, rng=rng)).answer[row, col] for _ in 1:num_iterations]
        mean = Nonneg.answerquery(q, data)[row, col]
        frac_lap = sum(x-> (x-mean) <= bound, lap) / num_iterations
        frac_g = sum(x -> (x-mean) <= bound, g) / num_iterations
        expected_lap = cdf(Laplace(0, 1/e), bound)
        # check some cdf values for the actual noise
        @test frac_lap ≈ expected_lap atol=atolstd * sqrt(expected_lap * (1-expected_lap)/num_iterations)
        expected_g = cdf(Normal(0, sqrt(1/(2*r))), bound)
        @test frac_g ≈ expected_g atol=atolstd * sqrt(expected_g * (1-expected_g)/num_iterations)
    end
end

@testset "correlations" begin
    if SKIP return end
    # test that noise is not correlated, in particular it is not identical
    # for the identity query
    #if SKIP return end
    atolstd = 4 # tolerate at most 4 std of error in empirical checks
    num_rows=1_000_000
    num_cols=2
    epsilon = 1.2
    rho = 0.4
    seed1 = 1
    rng1 = MersenneTwister(seed1)
    q1 = Nonneg.IdentityQuery()
    data = zeros(num_rows, num_cols)
    lap = Nonneg.laplace_mechanism(q1, data=data, epsilon=epsilon, rng=rng1).answer
    g = Nonneg.gaussian_mechanism(q1, data=data, rho=rho, rng=rng1).answer
    cor_lap = cor(lap[:, 1], lap[:, 2])
    cor_g = cor(g[:, 1], g[:, 2])
    @test cor_lap ≈ 0 atol = atolstd * sqrt(0.5/num_rows)
    @test cor_g ≈ 0 atol = atolstd * sqrt(0.5/num_rows)
end

@testset "compare single and multiple versions" begin
    if SKIP return end
    data = [1 2 3; 4 5 6; 7 8 9]
    query = Nonneg.IdentityQuery()
    epsilon = 1.2
    rho = 0.4
    seed = 1
    rng = MersenneTwister(seed)
    lap1 = Nonneg.laplace_mechanism(query, data=data, epsilon=epsilon; rng=rng)
    rng = MersenneTwister(seed)
    lap2 = Nonneg.laplace_mechanism([query], data=data, epsilon=epsilon; rng=rng)[1]
    @test lap1.answer == lap2.answer
    rng = MersenneTwister(seed)
    g1 = Nonneg.gaussian_mechanism(query, data=data, rho=rho; rng=rng)
    rng = MersenneTwister(seed)
    g2 = Nonneg.gaussian_mechanism([query], data=data, rho=rho; rng=rng)[1]
    @test g1.answer == g2.answer
end
