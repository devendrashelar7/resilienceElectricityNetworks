

# resiliency
begin
    N = 12; setGlobalParameters(); filename = string("N", N);

    global trial = false;
    global printPerformance = true;
    ResWC = 35;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta)); println("Resilience: ", 100(1-objv/LLCmax));
    println("# loads tripped: ",sum(kcv));

end

begin
    m = Model(solver=GurobiSolver(OutputFlag=0))
    @variable(m, x >= 1)
    @objective(m, :Min, x+1)
    solve(m)
    println(getobjectivevalue(m))
end
# resiliency
begin
    N = 12; setGlobalParameters(); filename = string("N", N);

    global trial = false;
    global printPerformance = true;
    ResWC = 35;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getBendersMethodCascadeForSocp(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta)); println("Resilience: ", 100(1-objv/LLCmax));
    println("# loads tripped: ",sum(kcv));

end

begin
    N = 12; setGlobalParameters(); filename = string("N", N);

    global trial = true;
    global printPerformance = true;
    ResWC = 20;
    LLCreq = (1 - ResWC/100) * LLCmax;
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta)); println("Resilience: ", 100(1-objv/LLCmax));

end

socpIneqtbar
begin
    m = Model(solver=GurobiSolver(OutputFlag=0))
    @variable(m, x)
    @constraint(m, con[i=1:2], x >= i)
    println(typeof(con))
    println(typeof(con[1]))
end


begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    getBendersCascadeMinNodesLLC()
end

begin
    N = 69; setGlobalParameters();
    println(getDGUnaffected());
end

begin
    fig = figure("Voltage Vs Distance",figsize=(5,5)); clf();
    ax = gca(); grid("on"); ms = 10; lw = 2; fs = 12;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    dist = copy(nodes); dist[9:12] = collect(5:8);
    plot(dist, vv, label="voltage", color="green", linewidth=lw,
    linestyle="", marker="s", markersize=ms);
    # plot(100oMs/N, 100oobjvs/LLCmax, label="uncontrolled", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    xlabel("Distance from substation",fontdict=font1);
    ylabel("Voltage in p.u.",fontdict=font1);

    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"voltageLTC",filename,".pdf"));
end

begin
    N = 3; delta = zeros(Int64,N); delta[3] = 1;
    # delta[DG] = 1;
    # N = 12; delta = zeros(Int64,N); delta[9:12] = 1;
    setGlobalParameters(); deltaV0 = 0;
    resourceResponse = false;
    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    printResults(pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv)
end

begin
    N = 12; delta = zeros(Int64,N);
    # delta[DG] = 1;
    # N = 12; delta = zeros(Int64,N); delta[9:12] = 1;
    setGlobalParameters(); deltaV0 = 0;
    resourceResponse = false;
    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    printResults(pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv)

    N = 12; delta = ones(Int64,N);
    # delta[DG] = 1;
    # N = 12; delta = zeros(Int64,N); delta[9:12] = 1;
    # setGlobalParameters(); deltaV0 = 0;
    resourceResponse = false;
    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    printResults(pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv)
end

begin
    N = 12; delta = zeros(Int64,N);
    # delta[DG] = 1;
    # N = 12; delta = zeros(Int64,N); delta[9:12] = 1;
    setGlobalParameters(); deltaV0 = 0;
    # resistance=ov; reactance = ov
    # pcmax = ov; qcmax = ov
    resourceResponse = false;
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveSocpModelWithCascade(delta, deltaV0);
    printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)

    N = 12; delta = ones(Int64,N);
    # delta[DG] = 1;
    # N = 12; delta = zeros(Int64,N); delta[9:12] = 1;
    # setGlobalParameters(); deltaV0 = 0;
    resourceResponse = false;
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveSocpModelWithCascade(delta, deltaV0);
    printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)
end

begin
    N = 12; delta = zeros(Int64,N); setGlobalParameters(); deltaV0 = 0;
    for deltaV0 in [0]
        j = deltaV0 + 2;
        for M = 0:0
            i = M+1;
            delta = zeros(Int64, N); delta[end-M+1:end] = 1;
            println("nodes attacked : ", nodes[delta.==1]);
            pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, oobjv,
            lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltaV0);

            # oobjvs[i,j] = oobjv;
        end
    end
end

# bruteForceComputationStorageSocp
begin
    N = 12; setGlobalParameters(); filename = "N12";
    bfobjvs = zeros(Float64,0); bfkcvs = zeros(Int64,N,0); bfkgvs = zeros(Int64,N,0);  bfsumkcvs = zeros(Int64,0); bfsumkgvs = zeros(Int64,0); bfDeltas = zeros(Int64,N,0); bfMs = zeros(Int64,0);

    for M = 3:3
        println("M : ", M);
        bestDelta, bestDeltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv =  getBruteForceCascadeSocpAttack(M-1);

        bfobjvs = vcat(bfobjvs,objv); bfkcvs = hcat(bfkcvs, kcv); bfkgvs = hcat(bfkgvs, kgv); bfsumkcvs = vcat(bfsumkcvs, sum(kcv)); bfsumkgvs = vcat(bfsumkgvs, sum(kgv)); bfDeltas = hcat(bfDeltas, bestDelta); bfMs = vcat(bfMs, M-1);
        # bfdeltaV0[M] = bestDeltaV0;
    end

    # jldopen(string(mypath,"bruteForceCascadeSocp",filename,".jld"), "w") do file
    #     write(file, "bfMs", bfMs);
    #     write(file, "bfobjvs", bfobjvs);
    #     # write(file, "bfsumkmvs", bfsumkmvs);
    #     write(file, "bfkcvs", bfkcvs);
    #     write(file, "bfkgvs", bfkgvs);
    #     # write(file, "bfkmvs", bfkmvs);
    #     write(file, "bfDeltas", bfDeltas);
    #     # write(file, "bfdeltaV0", bfdeltaV0);
    # end

end



# bruteForceComputationStorageSocp
begin
    N = 12; setGlobalParameters(); filename = "N12";
    bfobjvs = zeros(Float64,0); bfkcvs = zeros(Int64,N,0); bfkgvs = zeros(Int64,N,0);  bfsumkcvs = zeros(Int64,0); bfsumkgvs = zeros(Int64,0); bfDeltas = zeros(Int64,N,0); bfMs = zeros(Int64,0);

    for M = 4:4
        println("M : ", M);
        tic();
        bestDelta, bestDeltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv =  getBruteForceCascadeSocpAttack(M-1);
        println("Time taken = ", toq());
        withCuts = false;
        tic();
        bestDelta, bestDeltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv1, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv =  getBruteForceCascadeSocpAttack(M-1,withCuts);
        println("Time taken = ", toq());
        println("objvs = ", objv, " ", objv1)
    end
end

# bruteVsBenders
begin
    N = 12; setGlobalParameters(); #v0dis = 0.03;
    printPerformance = false;  filename = string("N", N);
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs, bfdeltaV0;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    afn = string(mypath,"bruteForceCascadeSocp",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"]; bfDeltas = load(afn)["bfDeltas"];

    for i = 1:length(bfobjvs)
        LLCreq = bfobjvs[i];

        delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv, ellv = getBendersMethodCascadeForSocp(LLCreq);

        println("Wanted cardinality : ", bfMs[i]);
        println("Got cardinality : ", sum(delta));
        println("LLCreq : ", bfobjvs[i]);
        println("Benders attack : ", nodes[delta.==1]);
        println("Optimal attack : ", nodes[bfDeltas[:,i].==1]);
        println("optimal objective: ", bfobjvs[i]);
    end

end

println(combinations(collect(1:6),2))
@everywhere begin
  using Combinatorics
  combs = collect(combinations(collect(1:6),2));
end
y = @parallel hcat for i in combs
  println(i)
  i
end
using DataStructures

println(subsets(collect(1:6),2))

x = combinations([1,2,3],2) |> collect
println(combs[2])
println(length(combs))

begin
  @everywhere function myCompare(tupleA, tupleB)
    (objvA, deltaA) = tupleA;
    (objvB, deltaB) = tupleB;
    if objvA > objvB
      return tupleA
    else
      return tupleB
    end
  end

  myCompare([2,4],[4,5])

  @time y = @parallel myCompare for i in combs
    i
  end

end

begin
    m = Model(solver=MosekSolver());
    p1 = 0.3; q1 = 0.1; r1 = 0.01; x1= 0.01
    @variable(m, nu0);
    @variable(m, P1);@variable(m, Q1);@variable(m, nu1); @variable(m, ell1);

    @constraint(m, substation, nu0 == 1)
    @constraint(m, realCons, P1 == p1 + r1*ell1); @constraint(m, reacCons, Q1 == q1 + x1*ell1);
    @constraint(m, voltDrop, nu1 == nu0 - 2*(r1*P1+x1*Q1) + (r1^2+x1^2)*ell1);

    @constraint(m, ohm, norm([nu0; ell1; P1; Q1]) <= nu0+ell1);


    @objective(m, :Min, ell1);
    println(m)
    solve(m);

    println(getvalue([nu1 ell1 P1 Q1]));


    println(getdual([realCons reacCons voltDrop]))
    println(getdual(ohm));

    println("obj = ", getobjectivevalue(m));
    println(getdual(voltDrop)*0+p1*getdual(realCons)+q1*getdual(reacCons) + getdual(substation)*1)
end
