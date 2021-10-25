@testset "queries" begin
    if SKIP return end
    data1 = collect(1:10)
    a = 2
    b = 3
    c = 4
    data2 = zeros(a,b,c)
    for i in 1:a
        for j in 1:b
            for k in 1:c
                data2[i,j,k]=i+j+k
            end
        end
    end

    idquery = Nonneg.IdentityQuery()
    sumquery = Nonneg.SumQuery()
    marg12 = Nonneg.MarginalQuery(numdims=3, keep=[1,2])
    marg3 = Nonneg.MarginalQuery(numdims=3, keep=3)

    @test Nonneg.l1sens(idquery) == 1
    @test Nonneg.l2sens(idquery) == 1
    @test Nonneg.l1sens(sumquery) == 1
    @test Nonneg.l2sens(sumquery) == 1
    @test Nonneg.l1sens(marg12) == 1
    @test Nonneg.l2sens(marg12) == 1
    @test Nonneg.l1sens(marg3) == 1
    @test Nonneg.l2sens(marg3) == 1

    @test Nonneg.answerquery(idquery, data1) == data1
    @test Nonneg.answerquery(idquery, data2) == data2

    @test Nonneg.answerquery(sumquery, data1) == [55]
    @test Nonneg.answerquery(sumquery, data2) == [a*b*c*(a+1)/2 + a*b*c*(c+1)/2 + a*b*c*(b+1)/2]

    @test Nonneg.answerquery(marg12, data2) == [c*(i+j)+c*(c+1)/2 for i in 1:a, j in 1:b]
    @test Nonneg.answerquery(marg3, data2) == [a*b*k + b*a*(a+1)/2 + a*b*(b+1)/2 for k in 1:c]

end

@testset "mask and filter query" begin
    if SKIP return end
    tol=1e-7
    r=4
    c=7
    d=3
    dataset = rand(r,c,d)
    mask_rc = rand(r,c)
    mask_d = rand(d)

    marg_rc = Nonneg.MarginalQuery(numdims=3, keep=[1,2])
    marg_d = Nonneg.MarginalQuery(numdims=3, keep=3)

    maskquery_rc = Nonneg.MaskQuery(marg_rc, mask_rc)
    maskquery_d = Nonneg.MaskQuery(marg_d, mask_d)

    filterquery_rc = Nonneg.FilterQuery(marg_rc, mask_rc .>= 0.5)
    filterquery_d = Nonneg.FilterQuery(marg_d, mask_d .>= 0.5)

    expected_rc = sum(dropdims(sum(dataset, dims=3), dims=3) .* mask_rc)
    expected_d = sum(dropdims(sum(dataset, dims=[1,2]), dims=(1, 2)) .* mask_d)

    got_rc = Nonneg.answerquery(maskquery_rc, dataset)
    got_d = Nonneg.answerquery(maskquery_d, dataset)

    @test sum(abs, expected_rc .- got_rc) ≈ 0 atol=tol
    @test sum(abs, expected_d .- got_d) ≈ 0 atol=tol

    expected_filter_rc = dropdims(sum(dataset, dims=3), dims=3) .* (mask_rc .>= 0.5)
    expected_filter_d = dropdims(sum(dataset, dims=[1,2]), dims=(1, 2)) .* (mask_d .>= 0.5)

    got_filter_rc = Nonneg.answerquery(filterquery_rc, dataset)
    got_filter_d = Nonneg.answerquery(filterquery_d, dataset)

    @test sum(abs, expected_filter_rc - got_filter_rc) ≈ 0 atol=tol
    @test sum(abs, expected_filter_d - got_filter_d) ≈ 0 atol=tol

    @test Nonneg.l1sens(maskquery_rc) <= maximum(abs.(mask_rc))
    @test Nonneg.l1sens(maskquery_d) <= maximum(abs.(mask_d))
    @test Nonneg.l2sens(maskquery_rc) <= maximum(abs.(mask_rc))
    @test Nonneg.l2sens(maskquery_d) <= maximum(abs.(mask_d))

    @test Nonneg.l1sens(filterquery_rc) <= 1
    @test Nonneg.l1sens(filterquery_d) <= 1
    @test Nonneg.l2sens(filterquery_rc) <= 1
    @test Nonneg.l2sens(filterquery_d) <= 1

end
