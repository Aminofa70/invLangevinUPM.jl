using Test

@testset "inverse Langevin function by Warner" begin
    ILf = invWarner()
    out = ILf(0.9) 
    expected_value = 14.210526315789478
    @test out ≈ expected_value atol=1e-10
end

##########################################################
@testset "inverse Langevin function by Puso" begin
    ILf = invPuso()
    out = ILf(0.9) 
    expected_value = 9.963099630996314
    @test out ≈ expected_value atol=1e-10
end

##########################################################
@testset "inverse Langevin function by Cohen" begin
    ILf = invCohen()
    out = ILf(0.9) 
    expected_value = 10.373684210526319
    @test out ≈ expected_value atol=1e-10
end

##########################################################
@testset "inverse Langevin function by Jedynak" begin
    ILf = invJedynak()
    out = ILf(0.9) 
    expected_value = 10.131192660550456
    @test out ≈ expected_value atol=1e-10
end

##########################################################
@testset "inverse Langevin function by Petrosyan" begin
    ILf = invPetrosyan()
    out = ILf(0.9) 
    expected_value = 9.988638025926525
    @test out ≈ expected_value atol=1e-10
end

##########################################################
@testset "inverse Langevin function by modified Jedynak" begin
    ILf = invJedynak_modified()
    out = ILf(0.9) 
    expected_value = 9.997227885887092
    @test out ≈ expected_value atol=1e-10
end

##########################################################
@testset "inverse Langevin function by Marchi and Arruda" begin
    ILf = invMarchiArruda()
    out = ILf(0.9) 
    expected_value = 9.910674957761998
    @test out ≈ expected_value atol=1e-10
end

##########################################################
@testset "inverse Langevin function by BenitezMontans" begin
    ILf = invBenitezMontans()
    out = ILf(0.9) 
    expected_value = 9.999999587768945
    @test out ≈ expected_value atol=1e-10
end
