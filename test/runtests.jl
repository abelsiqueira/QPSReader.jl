using QPSReader, NLPModels, NLPModelsIpopt, LinearOperators, QuadraticModels

problems = ["GENHS28.SIF"; "HS76.SIF"; "QPTEST.SIF"]
objectives = [9.2717369e-01; -4.6818182e+00; 4.3718750e+00]

@testset "qpsreader" begin
  for (k, p) in enumerate(problems)
    qps = readqps("../$p")
    nlp = QuadraticModel(qps.c, qps.Q, opHermitian(qps.Q), qps.A, qps.lcon, qps.ucon, qps.lvar, qps.uvar, c0=qps.c0)

    output = ipopt(nlp, print_level=0)

    @test output.dual_feas < 1e-6
    @test output.primal_feas < 1e-6
    @test norm(output.objective - objectives[k]) < 1e-6
  end
end
