using LazySets: linear_map_inverse

for N in [Float64, Rational{Int}, Float32]
    # random universe
    rand(Universe)

    U = Universe{N}(2)
    B = BallInf(ones(N, 2), N(1))

    # universe is an absorbing element for
    # - convex hull
    @test CH(B, U) == CH(U, B) == CH(U, U) == U
    cha = ConvexHullArray([B, N(2) * B, N(3) * B])
    @test CH(cha, U) == CH(U, cha) == U
    # - union
    @test B ∪ U == U ∪ B == U ∪ U == U
    # - Minkowski sum
    # TODO requires #1099
    # @test B ⊕ U == U ⊕ B == U ⊕ U == U
    # msa = MinkowskiSumArray([B, N(2) * B, N(3) * B])
    # @test msa ⊕ U == U ⊕ msa == U

    # universe is a neutral element for
    # - intersection
    @test B ∩ U == U ∩ B == B
    @test U ∩ U == U
    ia = IntersectionArray([B, N(2) * B, N(3) * B])
    @test ia ∩ U == U ∩ ia == ia

    # dim
    @test dim(U) == 2

    # copy
    U2 = copy(U)
    @test U2 isa Universe{N} && dim(U2) == 2

    # support function and support vector
    @test ρ(N[0, 1], U) == N(Inf)
    @test ρ(N[0, 0], U) == zero(N)
    @test σ(N[0, 1], U) == N[0, Inf]
    @test σ(N[-1, 2], U) == N[-Inf, Inf]

    # boundedness
    @test !isbounded(U)

    # ispolyhedral
    @test ispolyhedral(U)

    # membership
    @test N[0, 0] ∈ U
    @test_throws AssertionError N[0] ∈ U

    # emptiness
    @test !isempty(U)

    # universality
    @test isuniversal(U) && isuniversal(U, true) == (true, N[])

    # an_element
    @test an_element(U) ∈ U

    # constrained_dimensions
    @test constrained_dimensions(U) == Int[]

    # traits
    @test !isboundedtype(Universe)
    @test !isoperationtype(Universe)

    # norm/radius/diameter functions
    @test_throws ErrorException norm(U)
    @test_throws ErrorException radius(U)
    @test_throws ErrorException diameter(U)

    # translation
    @test translate(U, N[1, 2]) == U
    @test translate!(U, N[1, 2]) == U

    # constraints / constraints_list / tosimplehrep
    @test collect(constraints(U)) == constraints_list(U) == Vector{N}()
    A, b = tosimplehrep(U)
    @test A == zeros(N, 0, 2) && b == N[]

    # concrete intersection
    @test intersection(B, U) == intersection(U, B) == B
    @test intersection(U, U) == U

    # intersection emptiness
    res, w = isdisjoint(B, U, true)
    @test !isdisjoint(B, U) && !res && w ∈ B && w ∈ U
    res, w = isdisjoint(U, B, true)
    @test !isdisjoint(U, B) && !res && w ∈ B && w ∈ U
    res, w = isdisjoint(U, U, true)
    @test !isdisjoint(U, U) && !res && w ∈ U
    E = EmptySet{N}(2)
    res, w = isdisjoint(E, U, true)
    @test isdisjoint(E, U) && res && w == N[]
    res, w = isdisjoint(U, E, true)
    @test isdisjoint(U, E) && res && w == N[]

    # subset
    res, w = ⊆(B, U, true)
    @test B ⊆ U && res && w == N[]
    res, w = ⊆(U, B, true)
    @test U ⊈ B && !res && w ∉ B
    res, w = ⊆(U, U, true)
    @test U ⊆ U && res && w == N[]

    # inverse linear map
    M = ones(N, 2, 3)
    @test linear_map_inverse(M, U) == Universe{N}(3)

    # projection
    @test project(Universe{N}(5), [1, 4, 5]) == Universe{N}(3)

    # complement
    C = complement(U)
    @test C == EmptySet{N}(2) && C isa EmptySet{N}

    # permutation
    @test permute(U, [1, 2]) == permute(U, [2, 1]) == U

    # reflect
    @test reflect(U) == U

    # scale/scale!
    @test scale(N(1), U) == scale(N(-1), U) == U
    Z = scale(N(0), U)
    @test Z isa ZeroSet{N} && Z == ZeroSet{N}(dim(U))
    U2 = copy(U)
    scale!(N(1), U2)
    @test U2 == U
    @test_throws ArgumentError scale!(N(0), U2)
end

for N in [Float64, Float32]
    U = Universe{N}(2)

    # rationalize
    U2 = rationalize(U)
    @test U2 isa Universe{Rational{Int}} && dim(U2) == 2
    @test_throws MethodError rationalize(U2)
end

# default Float64 constructor
@test Universe(2) == Universe{Float64}(2)
