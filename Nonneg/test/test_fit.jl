@testset "no noise fitting" begin
    if SKIP return end
    a = 20
    b = 30
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
    dists = [Normal(0, 1), Laplace(0, 2), Normal(0,3), Laplace(0, 4)]
    answers = [Nonneg.NoisyAnswer(Nonneg.answerquery(q, data), d) for (q,d) in zip(queries, dists)]

    (ols_ans, status1) = Nonneg.olsfit(size(data), queries, answers)
    (nnls_ans, status2) = Nonneg.olsfit(size(data), queries, answers, nonneg=true)

    @test Nonneg.is_optimal(status1)
    @test Nonneg.is_optimal(status2)
    @test sum(abs, ols_ans - data) ≈ 0 atol=0.001
    @test sum(abs, ols_ans - data) ≈ 0 atol=0.001
end

@testset "nonneg" begin
    if SKIP return end
    answer = reshape(collect(-100:99), (40, 5))
    dist = Normal(0,1)
    q1 = Nonneg.IdentityQuery()
    (nnls_ans, status) = Nonneg.olsfit(size(answer), [q1], [Nonneg.NoisyAnswer(answer, dist)], nonneg=true)
    @test Nonneg.is_optimal(status)
    @test sum(abs, max.(answer, 0) - nnls_ans) ≈ 0 atol=0.001
end

@testset "ols" begin
    if SKIP return end
    m = 5
    adder = m * 3
    ans_id = collect(-5:5) .+ 0.0
    ans_sum = sum(ans_id) + adder
    n = length(ans_id)

    expected_sum = (n * ans_sum + sum(ans_id)) / (n+1)
    expected_ans = (expected_sum - sum(ans_id))/n .+ ans_id

    queries = [Nonneg.IdentityQuery(), Nonneg.SumQuery()]
    noisy = [Nonneg.NoisyAnswer(ans_id, Normal(0,1)), Nonneg.NoisyAnswer(ans_sum, Normal(0, 1))]
    (ols_ans, status) = Nonneg.olsfit(size(ans_id), queries, noisy)

    @test Nonneg.is_optimal(status)
    @test sum(abs, ols_ans-expected_ans) ≈ 0 atol=0.001
end

@testset "constrained_ols" begin
    if SKIP return end
    tolerance = 0.01
    n = 3
    queries = [Nonneg.IdentityQuery(), Nonneg.SumQuery()]
    mydist = Normal(0, 4)
    noisyid = rand(mydist, n)
    mysum = 7
    noisy = [Nonneg.NoisyAnswer(noisyid, mydist), Nonneg.NoisyAnswer([mysum], Normal(0,0))]
    expected = (mysum - sum(noisyid)) / n  .+ noisyid
    (ols_ans, status) = Nonneg.olsfit(size(noisyid), queries, noisy)
    @test Nonneg.is_optimal(status)
    @test sum(abs, expected - ols_ans) ≈ 0 atol=tolerance
end

@testset "olsweight" begin
    if SKIP return end
    m = 5
    adder = m * 3
    ans_id = collect(-5:5) .+ 0.0
    ans_sum = sum(ans_id) + adder
    n = length(ans_id)

    weights = [(1,1), (2,3)]

    for (w1, w2) in weights
        expected_sum = (w2 * n * ans_sum + w1 * sum(ans_id)) / (w2 * n + w1)
        expected_ans = (expected_sum - sum(ans_id))/n .+ ans_id

        queries = [Nonneg.IdentityQuery(), Nonneg.SumQuery()]
        noisy = [Nonneg.NoisyAnswer(ans_id, Normal(0,1)), Nonneg.NoisyAnswer(ans_sum, Normal(0, 1))]
        (ols_ans, status) = Nonneg.olsfit(size(ans_id), queries, noisy, reweights=[w1, w2])

        @test Nonneg.is_optimal(status)
        @test sum(abs, ols_ans-expected_ans) ≈ 0 atol=0.001

        noisy = [Nonneg.NoisyAnswer(ans_id, Normal(0,sqrt(2))), Nonneg.NoisyAnswer(ans_sum, Normal(0, sqrt(3)))]
        (ols_ans, status) = Nonneg.olsfit(size(ans_id), queries, noisy, reweights=[2*w1, 3*w2])

        @test Nonneg.is_optimal(status)
        @test sum(abs, ols_ans-expected_ans) ≈ 0 atol=0.001
    end
end

@testset "ols finegrained weight" begin
    if SKIP return end
    side = 3
    amount = side^2
    start = -5
    id_ans = reshape(collect(start:(start + amount -1)), (side,side))
    sum_ans = [sum(id_ans) - side^2]
    sum_weight = [2]
    id_weight = 1 ./ reshape(collect(1:amount), (side, side))
    queries = [Nonneg.IdentityQuery(), Nonneg.SumQuery()]
    noisy = [Nonneg.NoisyAnswer(id_ans, Normal(0,1)), Nonneg.NoisyAnswer(sum_ans, Normal(0, 1))]
    (olsans, status) = Nonneg.olsfit(size(id_ans), queries, noisy, reweights=[id_weight, sum_weight])

    w1 = sum(1 ./ id_weight)
    w2 = 1/sum(sum_weight)
    sum_estimate = (sum_ans[1] ./ w2  + sum(id_ans) ./ w1) / (1/w1 + 1/w2)

    discrepancy = sum_estimate - sum(id_ans)
    expected_ans = id_ans + (discrepancy ./ id_weight) ./ w1
    @test Nonneg.is_optimal(status)
    @test sum(abs, olsans-expected_ans) ≈ 0 atol=0.001
end

@testset "nearest_fit" begin
    if SKIP return end
    tolerance = 0.001
    iterations = 10
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

    dists_choices = [
             [Normal(0, 1), Laplace(0, 2), Normal(0,3), Laplace(0, 4)],
             [Normal(0, 1), Normal(0, 1), Normal(0,1), Normal(0, 1)]
            ]
    for _ in 1:iterations
        for dists in dists_choices
            queries = [q1, q2, q3, q4]
            noises = [rand(d, size(Nonneg.answerquery(q, data))...) for (d,q) in zip(dists, queries)]
            answers = [Nonneg.NoisyAnswer(Nonneg.answerquery(q, data) .+ pert, d) for (q,d,pert) in zip(queries, dists, noises)]
            negated = [Nonneg.NoisyAnswer(Nonneg.answerquery(q, data) .- pert, d) for (q,d,pert) in zip(queries, dists, noises)]

            (ols_ans, ub, status1) = Nonneg.nearestfit(size(data), vcat(queries, queries), vcat(answers, negated))

            @test Nonneg.is_optimal(status1)
            individual_ub = [maximum(abs.(Nonneg.answerquery(q, data) .- a.answer))/std(a.dist) for (q,a) in zip(queries, answers)]
            @test maximum(individual_ub) ≈ ub atol=tolerance
        end
    end
end

@testset "maxfit" begin
    if SKIP return end
    rng = MersenneTwister()
    tolerance=0.1
    datasize=1000
    data = zeros(datasize)
    n = 10
    data[1] = 100
    epsilon = sqrt(2)
    queries = [Nonneg.SumQuery(), Nonneg.IdentityQuery()]
    for _ in 1:n
        answers = Nonneg.laplace_mechanism(queries, data=data, epsilon=epsilon, rng=rng)
        (sol, (status1, status2)) = Nonneg.maxfit(size(data), queries, answers, nonneg=true)

        max_ans = [abs.(Nonneg.answerquery(q, sol) - a.answer) for (q,a) in zip(queries, answers)]
        max_max = maximum([maximum(x) for x in max_ans])

        data_ans = [abs.(Nonneg.answerquery(q, data) - a.answer) for (q,a) in zip(queries, answers)]
        datamax = maximum([maximum(x) for x in data_ans])
        #println("sol", sol)
        #println("datamax ",datamax)
        #println("maxmax ",max_max)
        #println("sol ", sol)
        #println("a1 ", answers[1].answer)
        #println("a2", answers[2].answer)
        #println("-----")
        @test Nonneg.is_optimal(status1)
        @test Nonneg.is_optimal(status2) || Nonneg.stopped_early(status2)
        @test datamax + tolerance >= max_max
    end
end

@testset "sequentialfit" begin
    if SKIP return end
    tol = 2e-3
    datashape = (2,3)
    q2 = Nonneg.IdentityQuery()
    q1 = Nonneg.SumQuery()
    N = Normal(0, 1)
    a2 = [3 2 1; -1 -2 -1]
    a1 = [9]

    expected = [4 3 2; 0 0 0]
    ans2 = Nonneg.NoisyAnswer(a2, N)
    ans1 = Nonneg.NoisyAnswer(a1, N)
    # fit sum query first
    (solution, statuses) = Nonneg.sequential_fit(datashape, [q1, q2], [ans1, ans2])
    @test sum(abs, solution - expected) ≈ 0 atol=tol
end
