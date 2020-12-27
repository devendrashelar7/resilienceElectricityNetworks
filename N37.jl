begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    global printPerformance = true;
    ResWC = 65;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltav0, pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    # delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getLinearAssistBendersCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta));
    println("% DGs attacked : ", 100sum(delta)/N)
    println("Resilience: ", 100(1-objv/LLCmax));
end

begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    println(DG)
    getDGUnaffected()
end
begin
    N = 36; setGlobalParameters(); filename = string("N", N);

    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv,
    lacv, lovrv, lcontv, lsdv = getEasyBenchmarkValues();

    println("voltage : ", vv);
    println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);

    println("pgv : ", pgv)

    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv,
    lacv, lovrv, lcontv, lsdv = getpgecontingencyValues();

    println("pge-contingency values")
    println("voltage : ", vv)
    println("objv : ", objv)
    println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);
    println("pgv : ", pgv);
    println("qgv : ", qgv);

    delta = zeros(Int64,N); deltaV0 = 0;
    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv,
    lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltaV0);

    println("Online cascade - no attack")
    println("voltage : ", vv)
    println("objv : ", objv)
    println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);
    println("pgv : ", pgv);
    println("qgv : ", qgv);

    delta[end-5+1:end] = 1;
    v0dis = 0.03; deltaV0 = 1;
    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv,
    lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltaV0);

    println("Online cascade - all attack")
    println("voltage : ", vv)
    println("objv : ", objv)
    println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);
    println("pgv : ", pgv);
    println("qgv : ", qgv);
end

println(LLCmax)
string("y",1)
println("y" + 1)

# bendersVsOnline
begin
    v0dis = 0.001;
    bnMs, bnminDeltav, bnobjvs, bnsumkcvs, bnsumkgvs,
    deltaV0s = getBendersCascadeMinNodesLLC();

    clf();  ms = 8; lw = 3;lw1 = 2;
    println(maximum(bnMs));

    oMsIndex = getIndexUniqueLargestIndex(bnMs);
    oMs = bnMs[oMsIndex]; odeltaV0s = deltaV0s[oMsIndex];
    ooDeltas = bnminDeltav[:,oMsIndex]; oobjvs = [];
    for i = 1:length(oMs)
        delta = ooDeltas[:,i]; deltaV0 = odeltaV0s[i];
        pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, oobjv, prv, qrv,
        lacv, lovrv, lcontv, olsdv = getOnlineCascade(delta, deltaV0);
        oobjvs = [oobjvs; olsdv];
    end

    fig = figure(L"$objv$ vs $M$",figsize=(5,5));
    ax = gca(); grid("on"); ms = 10; lw = 2; fs = 12;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bnMs/N, 100bnobjvs, label="BD", color="green", linewidth=lw,
    linestyle="--", marker="s", markersize=ms);
    plot(100oMs/N, 100oobjvs/LLCmax, label="uncontrolled", color="red", linewidth=lw,
    linestyle="--", marker="o", markersize=ms);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% Post-contingency cost",fontdict=font1);

    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"bendersVsOnlineConvex",filename,".pdf"));

end

#Investigate the master cuts
begin
    # println(bfobjvs); LLCreq = bfobjvs[1];
    LLCreq = 13000
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    # bnobjvs = [bnobjvs; objv]; bnMs = [bnMs; sum(delta)];
    println("LLCreq : ", LLCreq);
    println("benders : objv, deltaV0, sum(delta), ", objv, ", ", deltaV0, ", ",
    sum(delta));
    println("nodes attacked : ", nodes[delta.==1]);
end

# bruteVsBenders
begin
    N = 36; setGlobalParameters(); filename = "N36";
    bfobjvs = zeros(Float64,0); bfkcvs = zeros(Int64,N,0); bfkgvs = zeros(Int64,N,0);  bfsumkcvs = zeros(Int64,0); bfsumkgvs = zeros(Int64,0); bfDeltas = zeros(Int64,N,0); bfMs = zeros(Int64,0);

    # for M = 1:8
    for M = 3:3
        println("M : ", M);
        bestDelta, bestDeltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv =  getBruteForceCascadeAttack(M-1);

        bfobjvs = vcat(bfobjvs,objv); bfkcvs = hcat(bfkcvs, kcv); bfkgvs = hcat(bfkgvs, kgv); bfsumkcvs = vcat(bfsumkcvs, sum(kcv)); bfsumkgvs = vcat(bfsumkgvs, sum(kgv)); bfDeltas = hcat(bfDeltas, bestDelta); bfMs = vcat(bfMs, M-1);
        # bfdeltaV0[M] = bestDeltaV0;
    end

    jldopen(string(mypath,"bruteForceCascade",filename,".jld"), "w") do file
        write(file, "bfMs", bfMs);
        write(file, "bfobjvs", bfobjvs);
        write(file, "bfsumkcvs", bfsumkcvs);
        write(file, "bfsumkgvs", bfsumkgvs);
        # write(file, "bfsumkmvs", bfsumkmvs);
        write(file, "bfkcvs", bfkcvs);
        write(file, "bfkgvs", bfkgvs);
        # write(file, "bfkmvs", bfkmvs);
        write(file, "bfDeltas", bfDeltas);
        # write(file, "bfdeltaV0", bfdeltaV0);
    end

end

# bruteVsBenders
begin
    N = 36; setGlobalParameters(); #v0dis = 0.03;
    printPerformance = false;  filename = string("N", N);
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs, bfdeltaV0;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    afn = string(mypath,"bruteForceCascade",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"]; bfDeltas = load(afn)["bfDeltas"];

    i = 6; LLCreq = bfobjvs[i]; println(LLCreq);
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);

    println("Wanted cardinality : ", bfMs[i]);
    println("Got cardinality : ", sum(delta));
    println("LLCreq : ", bfobjvs[i]);
    println("Benders attack : ", nodes[delta.==1]);
    println("Optimal attack : ", nodes[bfDeltas[:,i].==1]);
    println("optimal objective: ", bfobjvs[i]);

end



begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    # v0dis = 0.01;
    delta = zeros(Int64,N);
    attack = [28, 32, 34, 35, 36];
    delta[attack] = 1; deltaV0 = 0;

    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);

    println("nodes attacked : ", nodes[delta.==1]);
    println("% objv: ", 100objv/LLCmax)
    println("objv: ", objv)
    # println(kcv);
    # println(vv);
    # println(pgv/srmax);
    # println(qgv/srmax);
end

# bruteVsBrute
begin
    N = 36; setGlobalParameters(); #v0dis = 0.03;
    printPerformance = false;  filename = string("N", N);
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs, bfdeltaV0;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    afn = string(mypath,"bruteForceCascade",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"]; bfDeltas = load(afn)["bfDeltas"]; bfDeltaV0 = load(afn)["bfdeltaV0"];
    bfobjvs1 = zeros(Float64, length(bfMs));

    for i = 1:length(bfMs)
        delta = bfDeltas[:,i]; deltaV0 = bfDeltaV0[i];
        println(deltaV0);
        pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
        bfobjvs1[i] = objv;

    end
    clf();
    plot(bfobjvs, bfobjvs1);
    savefig(string(mypath,"bruteVsBrute",filename,".pdf"));
end

# spatial structure of attacks
# bruteVsBrute
begin
    N = 36; setGlobalParameters(); #v0dis = 0.03;
    printPerformance = false;  filename = string("N", N);
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs, bfdeltaV0;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    afn = string(mypath,"bruteForceCascade",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"]; bfDeltas = load(afn)["bfDeltas"]; bfDeltaV0 = load(afn)["bfdeltaV0"];
    bfobjvs1 = zeros(Float64, length(bfMs));

    for i = 1:length(bfMs)
        delta = bfDeltas[:,i]; deltaV0 = bfDeltaV0[i];
        # println(deltaV0);
        println("nodes attacked: ", nodes[delta.==1]);
        # pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);

    end
end


begin
    N = 36;
    nodes = collect(1:N); srand(716);
    # Y = (nodes[randperm(N)]);
    println((nodes[randperm(N)])[1:Int64(N/2)])
end

#microgrid
begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    global printPerformance = true;
    ResWC = 80; global trial = false;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getBendersMethodMicrogrid(LLCreq);
    # delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getLinearAssistBendersCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta));
    println("% DGs attacked : ", 100sum(delta)/N)
    println("Resilience: ", 100(1-objv/LLCmax));
end

# bruteVsBenders
begin
    N = 36; setGlobalParameters(); #v0dis = 0.03;
    printPerformance = false;  filename = string("N", N);
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs, bfdeltaV0;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    afn = string(mypath,"bruteForceCascade",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"]; bfDeltas = load(afn)["bfDeltas"];

    TT = 11
    delta = zeros(Int64,N);
    delta[DG[bfDeltas[DG,TT] .== 1]] = 1
    deltav0 = 0
    println(delta)
    println(DG[bfDeltas[DG,TT] .== 1])
    println(DG)

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveModelWithCascade(delta, deltav0);

    println(kgv)
    println(nodes[kgv .== 1])

end
