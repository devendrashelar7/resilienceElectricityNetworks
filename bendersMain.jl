begin
    N = 12; setGlobalParameters(); filename = string("N", N);
    # v0dis = 0.01;
end



begin
    N = 12; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.03; M = 2; delta = zeros(Int64, N); delta[end-M+1:end] = 1;
    deltav0 = -1;

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltav0);

    clf(); myperm = [1,2,3,4,5,9,6,10,7,11,8,12]; ax = gca();

    myX = [1,2,3,4,5,6,7,8,5,6,7,8];
    plot(nodes[1:8], vv[1:8], linewidth=2, marker="s",markersize=15);
    plot(nodes[1:8], 0.95ones(nodes[1:8]), color="red",linewidth=3);
    plot(nodes[1:8], 1.05ones(nodes[1:8]), color="red",linewidth=3);
    plot(nodes[1:8], 1ones(nodes[1:8]), color="black",linewidth=3, linestyle="--");

    fs = 25;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);
    xlabel("Distance from substation",fontdict=font1);
    ylabel("Voltage",fontdict=font1); ylim(0.85,1.06);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    savefig(string(mypath,"tempVoltage",filename,".pdf"));

    printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, pgv, qgv);
    srmaxp = srmax; pc2sr = 2; pg2sr = 1; sd2sr = 0; se2sr = 1;
    modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr);
    println(round(vv,3))
    println("loads shed : ", nodes[kcv.==1])
    println("DG shed : ", nodes[kgv.==1])
end

begin
    LLCreq = 190.52
    delta, deltav0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    println(delta)
end


begin
    N = 12; setGlobalParameters(); filename = string("N", N);
    # v0dis = 0.01;
    delta = zeros(Int64,N);
    attack = [9,10,11,12];
    # attack = []
    delta[attack] = 1; deltav0 = 0;

    pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);

    println("nodes attacked : ", nodes[delta.==1]);
    println("% objv: ", 100objv/LLCmax)

    println(kcv);
    println(vv);
    println(prv/srmax);
    println(qrv/srmax);
end




# Benders all M
begin
    N = 24; setGlobalParameters(); filename = string("N",N);
    bnMs, bnminDeltav, bnobjvs, bnsumkcvs, bnsumkgvs = getBendersCascadeMinNodesLLC();
    println(bnMs);
    # bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);
    # bnsumkcvs = bnsumkcvs[bnind];
    # bnsumkgvs = bnsumkgvs[bnind];

    # jldopen(string(mypath,"bendersCascade",filename,".jld"), "w") do file
    #     write(file, "bnMs", bnMs);
    #     write(file, "bnobjvs", bnobjvs);
    #     write(file, "bnsumkcvs", bnsumkcvs);
    #     write(file, "bnsumkgvs", bnsumkgvs);
    #     # write(file, "bnsumkmvs", bnsumkmvs);
    # end

end

begin
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    jldopen(string(mypath,"bruteForceCascade",filename,".jld"), "r") do file
        bfMs = read(file, "bfMs");
        bfobjvs = read(file, "bfobjvs");
        bfsumkcvs = read(file, "bfsumkcvs");
        bfsumkgvs = read(file, "bfsumkgvs");
        # bfsumkmvs = read(file, "bfsumkmvs");
    end
    jldopen(string(mypath,"bendersCascade",filename,".jld"), "r") do file
        bnMs = read(file, "bnMs");
        bnobjvs = read(file, "bnobjvs");
        bnsumkcvs = read(file, "bnsumkcvs");
        bnsumkgvs = read(file, "bnsumkgvs");
        # bnsumkmvs = read(file, "bnsumkmvs");
    end
    # jldopen(string(mypath,"onlineCascade",filename,".jld"), "r") do file
    #   oobjvs = read(file, "oobjvs");
    #   sumokcvs = read(file, "sumokcvs");
    #   sumokgvs = read(file, "sumokgvs");
    #   okcvs = read(file, "okcvs");
    #   okgvs = read(file, "okgvs");
    # end
    println("read files");
    clf(); nd = length(bfobjvs); ms = 8; lw = 3;lw1 = 3;
    println(size(bfMs), " ", size(bfobjvs));
    println(bfMs);
    println(bnMs);
    plot(bfMs, bfobjvs/LLCmax, label=L"$brute$", color="red", linewidth=lw, linestyle="--", marker="x", markersize=ms)
    # plot(Ms, oobjvs/LLCmax, label=L"$cascade$", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    plot(bnMs, bnobjvs, label=L"benders", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Cost"); grid("on");
    # ylim(-50,maximum(oobjvs)+100);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"objCascadeMethods",filename,".pdf"));

    clf();
    plot(bfMs, bfsumkcvs, label="brute sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    # plot(bfMs, bfsumkgvs, label="brute sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    # plot(Ms, bfsumkmvs, label="brute sum(kmv)", color="red", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    # plot(Ms, sumokcvs, label="cascade sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="+", markersize=ms)
    # plot(Ms, sumokgvs, label="cascade sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="o", markersize=ms)
    plot(bnMs, bnsumkcvs, label="benders sum(kcv)", color="red", linewidth=lw1, linestyle="-", marker="x", markersize=ms)
    # plot(bnMs, bnsumkgvs, label="benders sum(kgv)", color="blue", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    # plot(bnminNodesv, bnsumkmvs, label="benders sum(kmv)", color="green", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Number"); grid("on");
    ylim(0,N+4);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"kckgCascademethods",filename,".pdf"));

end

begin
    N = 6; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.0;
    levp1 = 10;
    levp1 = length(SDI)+1;
    bfMs = zeros(Int64, levp1);
    bfobjvs = zeros(Float64, levp1); bfsumkcvs = zeros(Int64, levp1);
    bfsumkgvs = zeros(Int64, levp1); bfsumkmvs = zeros(Int64, levp1);
    bfkcvs = zeros(Int64, N, levp1); bfkmvs = zeros(Int64, N, levp1);
    bfkgvs = zeros(Int64, N, levp1); bfDeltas = zeros(Int64, N, levp1);

    for M = 1:levp1
        bestDelta, bestDeltaV0, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv,
        qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv,
        lmgv = getBruteForceMicrogridAttack(M-1);
        bfobjvs[M] = objv; bfsumkcvs[M] = sum(kcv); bfsumkgvs[M] = sum(kgv);
        bfsumkmvs[M] = sum(kmv);
        bfkcvs[:,M] = kcv; bfkgvs[:,M] = kgv; bfkmvs[:,M] = kmv;
        bfDeltas[:,M] = bestDelta; bfMs[M] = M-1;
    end

    bnMs, bnminDeltav, bnobjvs, bnsumkcvs, bnsumkgvs,
    bnsumkmvs = getBendersMicrogridMinNodesLLC();
    println(bnMs);
    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);
    bnsumkcvs = bnsumkcvs[bnind];
    bnsumkgvs = bnsumkgvs[bnind];
    bnsumkmvs = bnsumkmvs[bnind];

    clf();
    plot(bfMs, bfobjvs/LLCmax, label=L"$brute$", color="red", linewidth=lw, linestyle="--", marker="x", markersize=ms)
    # plot(Ms, oobjvs/LLCmax, label=L"$cascade$", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    plot(bnMs, bnobjvs, label=L"benders", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Cost"); grid("on");
    # ylim(-50,maximum(oobjvs)+100);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"test1",filename,".pdf"));

    clf();
    plot(bfMs, bfsumkcvs, label="brute sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    # plot(bfMs, bfsumkgvs, label="brute sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    plot(bfMs, bfsumkmvs, label="brute sum(kmv)", color="red", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    # plot(Ms, sumokcvs, label="cascade sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="+", markersize=ms)
    # plot(Ms, sumokgvs, label="cascade sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="o", markersize=ms)
    plot(bnMs, bnsumkcvs, label="benders sum(kcv)", color="red", linewidth=lw1, linestyle=":", marker="x", markersize=ms)
    # plot(bnMs, bnsumkgvs, label="benders sum(kgv)", color="blue", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    plot(bnMs, bnsumkmvs, label="benders sum(kmv)", color="green", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Number"); grid("on");
    bfl = length(bfMs); bnl = length(bnMs);
    ymax = max(maximum(bfsumkcvs), maximum(bfsumkgvs), maximum(bfsumkmvs),
    maximum(bnsumkcvs), maximum(bnsumkgvs), maximum(bnsumkmvs));
    ylim(0,ymax+4);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"test2",filename,".pdf"));

end

begin
    N = 6; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.0;

    bnMs, bnminDeltav, bnobjvs, bnsumkcvs, bnsumkgvs,
    bnsumkmvs = getBendersMicrogridMinNodesLLC();
    println(bnMs);
    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);
    bnsumkcvs = bnsumkcvs[bnind];
    bnsumkgvs = bnsumkgvs[bnind];
    bnsumkmvs = bnsumkmvs[bnind];

    bnMs1, bnminDeltav1, bnobjvs1, bnsumkcvs1, bnsumkgvs1,
    bnsumkmvs1 = getBendersMicrogridMinNodesLLC();
    println(bnMs);
    bnMs1, bnobjvs1, bnind1 = getUniqueFirstSortSecondArrays(bnMs1, bnobjvs1);
    bnsumkcvs1 = bnsumkcvs1[bnind1];
    bnsumkgvs1 = bnsumkgvs1[bnind1];
    bnsumkmvs1 = bnsumkmvs1[bnind1];

    clf();
    plot(bfMs, bfobjvs/LLCmax, label=L"$brute$", color="red", linewidth=lw, linestyle="--", marker="x", markersize=ms)
    # plot(Ms, oobjvs/LLCmax, label=L"$cascade$", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    plot(bnMs, bnobjvs, label=L"benders", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Cost"); grid("on");
    # ylim(-50,maximum(oobjvs)+100);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"test1",filename,".pdf"));

    clf();
    plot(bfMs, bfsumkcvs, label="brute sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    # plot(bfMs, bfsumkgvs, label="brute sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    plot(bfMs, bfsumkmvs, label="brute sum(kmv)", color="red", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    # plot(Ms, sumokcvs, label="cascade sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="+", markersize=ms)
    # plot(Ms, sumokgvs, label="cascade sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="o", markersize=ms)
    plot(bnMs, bnsumkcvs, label="benders sum(kcv)", color="red", linewidth=lw1, linestyle=":", marker="x", markersize=ms)
    # plot(bnMs, bnsumkgvs, label="benders sum(kgv)", color="blue", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    plot(bnMs, bnsumkmvs, label="benders sum(kmv)", color="green", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Number"); grid("on");
    bfl = length(bfMs); bnl = length(bnMs);
    ymax = max(maximum(bfsumkcvs), maximum(bfsumkgvs), maximum(bfsumkmvs),
    maximum(bnsumkcvs), maximum(bnsumkgvs), maximum(bnsumkmvs));
    ylim(0,ymax+4);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"test2",filename,".pdf"));

end

# Brute force all M
begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    levp1 = 10;
    levp1 = length(SDI)+1;
    levp1 = length(DG)+1;
    bfMs = zeros(Int64, levp1);
    bfobjvs = zeros(Float64, levp1); bfsumkcvs = zeros(Int64, levp1);
    bfsumkgvs = zeros(Int64, levp1); bfsumkmvs = zeros(Int64, levp1);
    bfkcvs = zeros(Int64, N, levp1);
    bfkgvs = zeros(Int64, N, levp1); bfDeltas = zeros(Int64, N, levp1);
    bfdeltav0 = zeros(Int64, levp1);

    for M = 1:levp1
        println("M : ", M);
        bestDelta, bestDeltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv,
        pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBruteForceCascadeAttack(M-1);
        # bestDeltaV0 = 0; bestDelta = zeros(Int64,N); bestDelta[DG[end-M+2:end]] = 1;
        # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
        # lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(bestDelta, bestDeltaV0);


        bfobjvs[M] = objv; bfkcvs[:,M] = kcv; bfkgvs[:,M] = kgv;
        bfsumkcvs[M] = sum(kcv); bfsumkgvs[M] = sum(kgv);
        bfDeltas[:,M] = bestDelta; bfMs[M] = M-1; bfdeltav0[M] = bestDeltaV0;
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
        write(file, "bfdeltav0", bfdeltav0);
    end

    # println(bfMs);
end

# bruteVsBenders loss/LLCmax vs cardinality
begin
    N = 24; setGlobalParameters(); #v0dis = 0.03;
    global epsilon = 0.1;
    printPerformance = false;  filename = string("N", N);
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs, bfdeltav0;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    jldopen(string(mypath,"bruteForceCascade",filename,".jld"), "r") do file
        bfMs = read(file, "bfMs");
        bfobjvs = read(file, "bfobjvs");
        bfdeltav0 = read(file, "bfdeltav0");
        println("Reached here.")

        bnDeltas = zeros(Int64,N,0);

        # v0dis = 0.0;
        bnobjvs = []; bnMs = [];
        for i = 1:length(bfobjvs)
            LLCreq = bfobjvs[i];
            println(LLCreq);
            # delta, deltav0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
            # objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getLinearAssistBendersCascade(LLCreq);
            delta, deltav0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
            # if objv >= LLCreq
            bnobjvs = [bnobjvs; objv]; bnMs = [bnMs; sum(delta)];
            bnDeltas = hcat(bnDeltas, delta);
            # end
        end

        clf(); nd = length(bfobjvs);
        println(size(bfMs), " ", size(bfobjvs));
        println(bfdeltav0);
        println(bfobjvs);
        println(bfMs);
        println(bnMs);
        ax = gca(); grid("on"); ms = 12; lw = 3; fs = 20;
        font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

        plot(100bfMs/N, 100bfobjvs/LLCmax, label="Optimal", color="red", linewidth=lw,
        linestyle="-", marker="o", markersize=ms)
        plot(100bnMs/N, 100bnobjvs/LLCmax, label="Benders", color="green", linewidth=lw,
        linestyle="--", marker="s", markersize=ms/2);
        xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% Post-contingency cost",fontdict=font1); # ylim(0,100);
        setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
        legend(loc = "upper left", fontsize=15); ## legend position
        new_position = [0.15,0.15,0.8,0.8] # Position Method 2
        ylim(0,50)
        ax[:set_position](new_position)
        savefig(string(mypath,"bruteVsBenders",filename,".pdf"),bbox_inches="tight");
        jldopen(string(mypath,"bendersCascade",filename,".jld"), "w") do file
            write(file, "bnMs", bnMs);
            write(file, "bnobjvs", bnobjvs);
            write(file, "bnDeltas", bnDeltas);
            # write(file, "bnsumkcvs", bnsumkcvs);
            # write(file, "bnsumkgvs", bnsumkgvs);
        end

    end
end

# print benders delta
begin
    N = 36; filename = string("N",N);
    afn = string(mypath,"bruteForceCascade",filename,".jld");
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];
    bfDeltas = load(afn)["bfDeltas"];

    afn = string(mypath,"bendersCascade",filename,".jld");
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"];
    bnDeltas = load(afn)["bnDeltas"];
    for i = 1:19
        println("i ",i)
        println(bfDeltas[DG,i]);
        println(bnDeltas[DG,i]);

        println("******************");
    end
end

# bruteVsBenders1 resilience vs cardinality
begin
    N = 36; filename = string("N",N); setGlobalParameters();
    afn = string(mypath,"bruteForceCascade",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    afn = string(mypath,"bendersCascade",filename,".jld"); # bendersCascadeFilename
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"];

    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bfMs/N, 100-100bfobjvs/LLCmax, label="Optimal", color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms)
    plot(100bnMs/N, 100-100bnobjvs/LLCmax, label="Benders", color="green", linewidth=lw,
    linestyle="--", marker="s", markersize=ms/2);
    xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    # legend(loc = "upper right", fontsize=15); ## legend position
    legend(loc = "lower left", fontsize=fs); ## legend position
    # new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    ylim(50,100)
    ylim(60,100)
    # ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersPres",filename,".pdf"),bbox_inches="tight");
end



# bruteVsBenders2
begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    afn = string(mypath,"bruteForceCascade",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    afn = string(mypath,"bendersCascade",filename,".jld"); # bendersCascadeFilename
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"];

    fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5));
    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 15;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bfobjvs/LLCmax, 100bfMs/N, label="Optimal", color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms)
    plot(100bfobjvs/LLCmax, 100bnMs/N, label="Benders", color="green", linewidth=lw,
    linestyle="--", marker="s", markersize=ms/2);
    ylabel("% Number of nodes attacked",fontdict=font1); xlabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    legend(loc = "upper left", fontsize=15); ## legend position
    new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    # ylim(0,50)
    ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBenders",filename,".pdf"));
end

# bruteVsBendersSocp2
begin
    N = 12; setGlobalParameters(); filename = string("N",N);
    afn = string(mypath,"bruteForceCascadeSocp",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    afn = string(mypath,"bendersCascadeMinNodesSocp",filename,".jld"); # bendersCascadeFilename
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"];

    fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5));
    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 15;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bfobjvs/LLCmax, 100bfMs/N, label="Optimal", color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms)
    plot(100bnobjvs/LLCmax, 100bnMs/N, label="Benders", color="green", linewidth=lw,
    linestyle="--", marker="s", markersize=ms/2);
    ylabel("% Number of nodes attacked",fontdict=font1); xlabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    legend(loc = "upper left", fontsize=15); ## legend position
    new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    # ylim(0,50)
    ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersSocp",filename,".pdf"));
end

using PyCall
using PyPlot
# bruteVsBendersEpsilons1
begin
    N = 36; setGlobalParameters();  filename = string("N", N);
    global epsilon = 20; printPerformance = false;
    afn = string(mypath,"bruteForceCascade",filename,".jld");
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 20;
    font1 = Dict("family"=>"serif", "color"=>"black", "weight"=>"normal", "size"=>fs);
    plot(100bfMs/N, 100-100bfobjvs/LLCmax, label="Optimal", color="red", linewidth=lw, linestyle="-");
    # , marker="o", markersize=ms/2);

    mycolors = ["black","green","blue","magenta"];
    mylinestyles = ["-","--","-.",":"];
    mymarkers = ["s","+","x","o"];
    epsilons = [10,1,0.1,0.05];
    for j = 1:length(epsilons)
        global epsilon = epsilons[j];
        bnDeltas = zeros(Int64,N,0); bnobjvs = []; bnMs = [];
        for i = 1:length(bfobjvs)
            LLCreq = bfobjvs[i];
            delta, deltav0, pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
            bnobjvs = [bnobjvs; objv]; bnMs = [bnMs; sum(delta)];
            bnDeltas = hcat(bnDeltas, delta);
        end
        mylabel = "Benders " * L"\epsilon = " * string(epsilon);
        plot(100bnMs/N, 100-100bnobjvs/LLCmax, label= mylabel, color=mycolors[j], linewidth=lw, linestyle=mylinestyles[j]);
        # , marker=mymarkers[j], markersize=ms);
    end

    println("Got here")
    xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    legend(loc = "lower left", fontsize=15); ## legend position
    new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    ylim(50,110)
    ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersEpsilons",filename,".pdf"), bbox_inches="tight");

end

begin
    using PyPlot
    plot([1 2 3 4]',[1 2 3 4]');
    fs = 18; ax = gca()
    font1 = Dict("family"=>"serif", "color"=>"black", "weight"=>"normal", "size"=>fs);
    xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    savefig(string(mypath,"random.pdf"), bbox_inches="tight");
end

# bruteVsBendersEpsilons2
begin
    N = 36; setGlobalParameters();  filename = string("N", N);
    global epsilon = 20; printPerformance = false;
    afn = string(mypath,"bruteForceCascade",filename,".jld");
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 20;
    font1 = Dict("family"=>"serif", "color"=>"black", "weight"=>"normal", "size"=>fs);
    plot(100-100bfobjvs/LLCmax, 100bfMs/N, label="Optimal", color="red", linewidth=lw, linestyle="-");
    # , marker="o", markersize=ms/2);

    mycolors = ["green","blue","magenta"];
    mylinestyles = ["--","-.",":"];
    mymarkers = ["s","+","x"];
    epsilons = [10,1,0.1];
    for j = 1:length(epsilons)
        global epsilon = epsilons[j];
        bnDeltas = zeros(Int64,N,0); bnobjvs = []; bnMs = [];
        for i = 1:length(bfobjvs)
            LLCreq = bfobjvs[i];
            delta, deltav0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
            bnobjvs = [bnobjvs; objv]; bnMs = [bnMs; sum(delta)];
            bnDeltas = hcat(bnDeltas, delta);
        end
        mylabel = "Benders" * L"\epsilon = " * string(epsilon);
        plot( 100-100bnobjvs/LLCmax, 100bnMs/N,label= mylabel, color=mycolors[j], linewidth=lw, linestyle=mylinestyles[j]);
        # , marker=mymarkers[j], markersize=ms);
    end

    ylabel("% Number of nodes attacked",fontdict=font1); xlabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    legend(loc = "upper left", fontsize=15); ## legend position
    new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    ylim(-5,50)
    ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersEpsilons",filename,".pdf"), bbox_inches="tight");

end

pwd()

# bruteVsBenders write, read and plot
begin
    global bnMs, bnobjvs
    jldopen(string(mypath,"bendersCascade",filename,".jld"), "w") do file
        write(file, "bnMs", bnMs);
        write(file, "bnobjvs", bnobjvs);
        # write(file, "bnsumkcvs", bnsumkcvs);
        # write(file, "bnsumkgvs", bnsumkgvs);
    end


end

begin
    abc = load(string(mypath,"bendersCascade",filename,".jld"))["bnMs"];
    println(abc)

end





begin
    afn = string(mypath,"bendersCascade",filename,".jld"); # bendersCascadeFilename
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"]; bnDeltas = load(afn)["bnDeltas"];
    println(bnMs);
    println(bnDeltas[:,6]);
    delta1 = bnDeltas[:,6];
    delta0 = zeros(delta1); delta0[8:12] = 1; deltav0 = 0;
    pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta0, deltav0);
    objv0 = objv;
    pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta1, deltav0);
    objv1 = objv;
    println(objv0," ", objv1)
end

begin
    println(bfobjvs); LLCreq = bfobjvs[1];
    delta, deltav0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    # bnobjvs = [bnobjvs; objv]; bnMs = [bnMs; sum(delta)];
    println("LLCreq : ", LLCreq);
    println("benders : objv, deltav0, sum(delta), ", objv, ", ", deltav0, ", ",
    sum(delta));
end

begin
    println("What");
    N = 12; attack = [7,8,11,12];
    # N = 6; attack = [2,4,6];

    setGlobalParameters(); M = 3;  filename = string("N", N);
    # bestDelta, bestDeltaV0, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv,
    # qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv,
    # lmgv = getBruteForceMicrogridAttack(M);


    bestDelta = zeros(Int64,N); bestDelta[end-M+1:end] = 1;
    # bestDelta[attack] = 1;

    # println(se2sr);
    deltav0 = -1; se2sr = semax[EG[1]]/srmax; println(se2sr);

    srmaxp = srmax; pc2sr = 4; pg2sr = 1; sd2sr = 0; se2sr = 1;
    modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr);
    sd2srs = linspace(0,5,20); m1 = length(sd2srs);
    v0diss = linspace(0,0.03,2); n1 = length(v0diss);
    println(sd2srs);
    objvs = Array(Float64,m1,n1); sumkmvs = Array(Int64,m1,n1);
    lacvs = Array(Float64,m1,n1); lovrvs = Array(Float64,m1,n1);
    lcontvs = Array(Float64,m1,n1); lsdvs = Array(Float64,m1,n1);
    legvs = Array(Float64,m1,n1); lmgvs = Array(Float64,m1,n1);
    sumkcvs = Array(Int64,m1,n1); sumkevs = Array(Int64,m1,n1);
    oobjvs = Array(Float64,m1,n1); osumkmvs = Array(Int64,m1,n1);
    osumkcvs = Array(Int64,m1,n1);

    for j = 1:length(v0diss)
        v0dis = v0diss[j];
        for i = 1:length(sd2srs)
            sd2sr = sd2srs[i]; pdmax[1:N] = sd2sr * srmax; qdmax = pdmax/3;
            println(i)
            # pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv,
            # qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv,
            # legv = evaluateSlaveModelWithMicrogrid(bestDelta, deltav0);
            pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(bestDelta, deltav0);
            objvs[i,j], sumkcvs[i,j] = objv, sum(kcv);
            lacvs[i,j], lovrvs[i,j], lcontvs[i,j], lsdvs[i,j] =
            lacv, lovrv, lcontv, lsdv;
            # sumkevs[i,j] = sum(kev); sumkmvs[i,j] = sum(kmv); lmgvs[i,j], legvs[i,j] = lmgv, legv;

            # pv, prv, qv, qrv, betav, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lovrv,
            # lcontv, lsdv = getOnlineCascade(bestDelta, deltav0);
            # oobjvs[i,j], osumkcvs[i,j], osumkmvs[i,j] = objv, sum(kcv), sum(kmv);
        end
    end

    lw = 2; ms = 10;
    clf(); plot(sd2srs, objvs/LLCmax, label="Loss (seq)", color="green", linewidth=lw,
    linestyle="-", marker="x", markersize=ms);
    # plot(sd2srs, oobjvs/LLCmax, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    savefig(string(mypath,"objVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, sumkmvs, label="Loss (seq)", color="green", linewidth=lw,
    linestyle="-", marker="x", markersize=ms);
    # plot(sd2srs, osumkmvs, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    savefig(string(mypath,"sumkmvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, sumkcvs, label="Loss", color="green", linewidth=lw,
    linestyle="-", marker="x", markersize=ms);
    # plot(sd2srs, osumkcvs, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    savefig(string(mypath,"sumkcvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, (lacvs + lovrvs + lcontvs)/LLCmax, label="Loss", color="green", linewidth=lw,
    linestyle="-", marker="x", markersize=ms);
    plot(sd2srs, lsdvs/LLCmax, label="Loss (onine)", color="red", linewidth=lw,
    linestyle="--", marker="o", markersize=ms);
    savefig(string(mypath,"lreg1Vssd2sr",filename,".pdf"));

end

begin
    clf(); plot(sd2srs, objvs/LLCmax, linewidth=lw); fs = 12;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);
    mylegend = [L"\Delta \nu = 0.00",L"\Delta \nu = 0.03",L"\Delta \nu = ",L"\Delta \nu = "];
    myxlabel = "Disturbance to reserves ratio";
    # plot(sd2srs, oobjvs, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    legend(mylegend, loc = "upper left", fontsize=15)
    xlabel(myxlabel,fontdict=font1); ylabel(L"L_{SD}",fontdict=font1)
    ax = gca(); grid("on");
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    savefig(string(mypath,"objVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, sumkmvs, linewidth=lw,
    linestyle="-", markersize=ms);ylim(0, length(MGL)+1);
    # plot(sd2srs, osumkmvs, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    # ax[:xaxis][:set_major_locator](Mx)
    legend(mylegend, loc = "upper left", fontsize=15)
    savefig(string(mypath,"sumkmvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, sumkcvs, linewidth=lw, markersize=ms); ylim(0, N+1);
    # plot(sd2srs, osumkcvs, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    legend(mylegend, loc = "upper left", fontsize=15)
    savefig(string(mypath,"sumkcvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, sumkevs, linewidth=lw); ylim(-1, length(EG) + 1); ylim(0, length(EG)+1);
    legend(mylegend, loc = "upper left", fontsize=15)
    # plot(sd2srs, osumkcvs, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    savefig(string(mypath,"sumkevVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, lcontvs, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15)
    savefig(string(mypath,"lcontvVssd2sr",filename,".pdf"));
    clf(); plot(sd2srs, legvs, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15)
    savefig(string(mypath,"legvVssd2sr",filename,".pdf"));

    # clf(); plot(sd2srs, lmgvs, linewidth=lw);
    # legend(mylegend, loc = "upper left", fontsize=15)
    # savefig(string(mypath,"lmgvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, lsdvs/LLCmax, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15)
    xlabel(myxlabel); ylabel(L"L_{SD}")
    savefig(string(mypath,"lsdvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, lcontvs + lsdvs, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15)
    savefig(string(mypath,"llcVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, (lacvs + lovrvs + lcontvs)/LLCmax, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15);
    # xlabel(L"\frac{pd^a}{\overline{sr}}");
    xlabel(myxlabel,fontdict=font1);
    ylabel(L"L^{GC}",fontdict=font1);
    ax = gca(); grid("on"); llfs = 8;
    setp(ax[:get_yticklabels](),fontsize=llfs); setp(ax[:get_xticklabels](),fontsize=llfs);
    legend(mylegend, loc = "upper left", fontsize=fs);
    savefig(string(mypath,"lreg1Vssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, (lacvs + lovrvs + lcontvs + lsdvs)/LLCmax, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15);
    plot(sd2srs, objvs/LLCmax, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15);
    savefig(string(mypath,"allVssd2sr",filename,".pdf"));


end

# seqVsOnline
begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    # srmaxp = srmax; pc2sr = 2.2; pg2sr = 1; sd2sr = 0; se2sr = 1;
    # modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr);

    v0dis = 0.02;
    delta = zeros(Int64, N);
    ludg = length(DG);
    deltav0 = -1;
    oobjvs = zeros(ludg+1,2);
    sobjvs = zeros(ludg+1,2);
    mobjvs = zeros(ludg+1,2);

    i = 0; j = 0;
    for deltav0 in [0,-1]
        j = deltav0 + 2;
        for M = 0:ludg
            i = M+1;
            delta = zeros(Int64, N); delta[DG[end-M+1:end]] = 1;
            # println("nodes attacked : ", nodes[delta.==1]);
            pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, sobjv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);
            # pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv,
            # qcv, mobjv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv,
            # legv = evaluateSlaveModelWithMicrogrid(delta,deltav0);

            sobjvs[i,j] = sobjv;
             # mobjvs[i,j] = mobjv;

        end
    end
    i = 0; j = 0;
    for deltav0 in [0,-1]
        j = deltav0 + 2;
        for M = 0:ludg
            i = M+1;
            delta = zeros(Int64, N); delta[DG[end-M+1:end]] = 1;
            println("nodes attacked : ", nodes[delta.==1]);
            pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, oobjv, pgv, qgv,
            lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltav0);

            oobjvs[i,j] = oobjv;
        end
    end

    clf();
    fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    ax = gca(); grid("on"); ms = 6; lw = 2; fs = 15;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    Ms = collect(0:length(DG));
    plot(100Ms/N, 100oobjvs[:,2]/LLCmax, color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms);
    plot(100Ms/N, 100sobjvs[:,2]/LLCmax, color="green",  linewidth=lw,
    linestyle="--", marker="s", markersize=ms);
    mylegend = [
    "Uncontrolled " * L"\Delta \mathrm{v}_0 = 0",
    "BD " * L"\Delta \mathrm{v}_0 = 0"];
    ylim(0,80);
     grid("on");

    legend(mylegend, loc = "upper left", fontsize=15);
    # println(oobjvs);
    # println(sobjvs);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% Post-contingency cost",fontdict=font1);
    savefig(string(mypath, "seqVsOnline_1_", filename,".pdf"));
    clf(); fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    # Ms = collect(0:N);
    plot(100Ms/N, 100oobjvs[:,1]/LLCmax, color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms);
    plot(100Ms/N, 100sobjvs[:,1]/LLCmax, color="green",  linewidth=lw,
    linestyle="--", marker="s", markersize=ms);
    mylegend = [
    "Uncontrolled " * L"\Delta \mathrm{v}_0 = " * string(v0dis),
    "Benders " * L"\Delta \mathrm{v}_0 = " * string(v0dis)];
    ylim(0,80);
    grid("on");

    legend(mylegend, loc = "upper left", fontsize=15);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% Post-contingency cost",fontdict=font1);
    savefig(string(mypath, "seqVsOnline_2_", filename,".pdf"));

    # printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    # objv, pgv, qgv);
    # srmaxp = srmax; pc2sr = 2; pg2sr = 1; sd2sr = 0; se2sr = 1;
    # modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr);
    # println(round(vv,3))
    # println("loads shed : ", nodes[kcv.==1])
    # println("DG shed : ", nodes[kgv.==1])
end
a = collect(0:N); println(a)

begin
    N = 12; setGlobalParameters(); filename = string("N",N);
    global pcmax, qcmax, v0dis
    v0diss = [0,0.01,0.02,0.03];
    scmax = (pcmax.^2 + qcmax.^2).^0.5; sgmax = (pgmax.^2 + qgmax.^2).^0.5;
    pfn = 10; pf = pi/2*collect(0:pfn)/pfn;
    deltav0 = -1; delta = zeros(Int64,N);
    N == 12? delta[DG[end-2:end]] = 1 : delta[DG[8:12] = 1];

    ppcl = zeros(Float64, length(pf), length(v0diss));
    for i = 1:pfn + 1
        pcmax = cos(pf[i]) * scmax; qcmax = sin(pf[i]) * scmax;
        for j = 1:length(v0diss)
            v0dis = v0diss[j];
            pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);
            ppcl[i,j] = objv;
        end
    end
    clf(); fig = figure("PPCL vs powerFactor",figsize=(5,5)); ax = gca(); grid("on");
    mycolors = ["green","red","blue","black"];
    linestyles = ["-","--","-.",":"]
    mylegends = L"\Delta \mathrm{v}_0 = " .* string.(v0diss); lw = 2; ms = 6;
    # plot(0,0);
    for j = 1:length(v0diss)
        plot(100cos.(pf), 100ppcl[:,j]/LLCmax, color=mycolors[j], linewidth=lw, linestyle=linestyles[j], marker="o", markersize=ms);
        # mylegends[j] = L"\Delta \mathrm{v}_0 = " * string(v0diss[j]);
    end
    # mylegend = [
    # "Uncontrolled " * L"\Delta \mathrm{v}_0 = 0.015",
    # "Benders " * L"\Delta \mathrm{v}_0 = 0.015"];
    # mylegends
    N == 12 ? ylim(0,120): ylim(0,40);
    grid("on");
    fs = 15;
    legend(mylegends, loc = "upper left", fontsize=fs);
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);
    xlabel("% Power factor",fontdict=font1);
    ylabel("% Post-contingency cost",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    new_position = [0.15,0.12,0.80,0.80] # Position Method 2
    ax[:set_position](new_position)
    ff = string(mypath, "ppclVsPowerFactor", filename,".pdf"); println(ff);
    savefig(ff);

end
string.([1, 2, 3])


begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.02;
    # M = 2; delta = zeros(Int64, N); delta[end-M+1:end] = 1;
    deltav0 = -1;
    bfDeltas = load(string(mypath,"bruteForceCascade",filename,".jld"))["bfDeltas"];
    delta = copy(bfDeltas[:,16]);

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltav0);

    println(sum(delta))
    println(sum(kgv))
end


# seqVsOnline
begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    # srmaxp = srmax; pc2sr = 2.2; pg2sr = 1; sd2sr = 0; se2sr = 1;
    # modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr);

    v0dis = 0.02;
    delta = zeros(Int64, N);
    ludg = length(UDG);
    deltav0 = -1;
    oobjvs = zeros(ludg+1,2); odeltas = zeros(Int64, N, ludg+1, 2);

    i = 0; j = 0;
    for M = 0:ludg
        for deltav0 in [-1,0]
            bestObjv = 0; bestDelta = zeros(Int64,N);
            for attackSet in collect(combinations(UDG,M))
                delta = zeros(Int64,N); delta[attackSet] = 1;

                pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltav0);

                if objv >= bestObjv
                    bestObjv = objv; bestDelta = copy(delta);
                end
            end
            oobjvs[M+1,deltav0+2] = bestObjv;
            odeltas[:,M+1,deltav0+2] = bestDelta;
        end
    end

    jldopen(string(mypath,"onlineCascade",filename,".jld"), "w") do file
        write(file, "oobjvs", oobjvs);
        write(file, "odeltas", odeltas);
    end
end

# Benders all M
begin
    N = 24; setGlobalParameters(); filename = string("N",N); v0dis = 0.02;

    deltav0 = -1; deltav0min = -1;
    bnMs, bnminDeltav, bnobjvs, bnsumkcvs, bnsumkgvs = getBendersCascadeMinNodesLLC();
    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);

    deltav0 = 0; deltav0min = 0;
    bnMs1, bnminDeltav, bnobjvs1, bnsumkcvs, bnsumkgvs = getBendersCascadeMinNodesLLC();
    bnMs1, bnobjvs1, bnind = getUniqueFirstSortSecondArrays(bnMs1, bnobjvs1);

    jldopen(string(mypath,"bendersOnlineCascade",filename,".jld"), "w") do file
        write(file, "bnMs", bnMs);
        write(file, "bnobjvs", bnobjvs);
        write(file, "bnMs1", bnMs1);
        write(file, "bnobjvs1", bnobjvs1);
    end
end

# Random attacks
begin
    N = 36; setGlobalParameters(); filename = string("N", N);
    # srmaxp = srmax; pc2sr = 2.2; pg2sr = 1; sd2sr = 0; se2sr = 1;
    # modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr);

    v0dis = 0.02;
    delta = zeros(Int64, N);
    ludg = length(UDG);
    deltav0 = -1;     nr = 50;
    robjvs = zeros(ludg+1,2,nr);
    oobjvs = zeros(ludg+1,2);
    i = 0; j = 0;
    for i = 1:nr
        for deltav0 in [-1,0]
            rUDG = UDG[randperm(ludg)];
            delta = zeros(Int64,N);

            for M = 0:ludg

                M > 0? delta[rUDG[M]] = 1 : nothing;
                # println(nodes[delta.==1])

                pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltav0);
                # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);
                robjvs[M+1,deltav0+2,i] = objv;
            end
        end
    end

    for i = 1:ludg+1
        for j = 1:2
            oobjvs[i,j] = maximum(robjvs[i,j,:]);
        end
    end

    jldopen(string(mypath,"onlineCascade",filename,".jld"), "w") do file
        write(file, "oobjvs", oobjvs);
    end
    jldopen(string(mypath,"randomCascade",filename,".jld"), "w") do file
        write(file, "robjvs", robjvs);
    end
end
begin
    afn = string(mypath,"randomCascade",filename,".jld");
    robjvs = load(afn)["robjvs"];
    size(robjvs)
end

xobjvs = rand(3,3); println(robjvs); robjvs = sort(robjvs, 2); println(robjvs); println(xobjvs); maximum(xobjvs[1,:])

# plot benders vs online
begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    afn = string(mypath,"bendersOnlineCascade",filename,".jld");
    bnMs = load(afn)["bnMs"]; bnMs1 = load(afn)["bnMs1"];
    bnobjvs = load(afn)["bnobjvs"]; bnobjvs1 = load(afn)["bnobjvs1"];
    afn = string(mypath,"onlineCascade",filename,".jld");
    oobjvs = load(afn)["oobjvs"]; #odeltas = load(afn)["odeltas"];
    v0dis = 0.02;
    afn = string(mypath,"randomCascade",filename,".jld");
    robjvs = load(afn)["robjvs"];

    clf();
    fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    ax = gca(); grid("on"); ms = 6; lw = 2; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    Ms = collect(0:length(DG));
    plot(100Ms/N, 100-100oobjvs[:,2]/LLCmax, color="red", linewidth=lw, linestyle="-", marker="o", markersize=ms);
    plot(100bnMs1/N, 100-100bnobjvs1/LLCmax, color="green",  linewidth=lw, linestyle="--", marker="s", markersize=ms);
    mylegend = [
    "Worst attack, response (b)",
    "Worst attack, response (c)",
    "Random attacks, response (b)"];
    ylim(20,100);
     grid("on");

    nr = 20;
    println(size(robjvs));
    robjvs = sort(robjvs, 3);
    alpha1 = 0.3; slw = 1; sms = 1;
    for i = 1:nr
        plot(100Ms/N, 100-100robjvs[:,2,i]/LLCmax, color="black", linewidth=slw, linestyle="--",alpha=alpha1);
    end

    legend(mylegend, loc = "lower left", fontsize=12);
    # println(oobjvs);
    # println(sobjvs);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% System resilience",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    savefig(string(mypath, "seqVsOnline_1_", filename,".pdf"),bbox_inches="tight");

    clf(); fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    # Ms = collect(0:N);
    plot(100Ms/N, 100-100oobjvs[:,1]/LLCmax, color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms);
    plot(100bnMs/N, 100-100bnobjvs/LLCmax, color="green",  linewidth=lw,
    linestyle="--", marker="s", markersize=ms);

    for i = 1:nr
        plot(100Ms/N, 100-100robjvs[:,1,i]/LLCmax, color="black", linewidth=slw, linestyle=":",alpha=alpha1);
    end

    mylegend = [
    "Worst attack, response (b)",
    "Worst attack, response (c)",
    "Random attacks, response (b)"];
    ylim(20,100);
    grid("on");

    legend(mylegend, loc = "lower left", fontsize=12);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% System resilience",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    # annotate("Value of !",
	# xy=[x;y,# Arrow tip
	# xytext=[x+dx;y+dy], # Text offset from tip
	# xycoords="data", # Coordinates in in "data" units
	# arrowprops=["facecolor"=>"black"])
    savefig(string(mypath, "seqVsOnline_2_", filename,".pdf"),bbox_inches="tight");
end


# Random attacks Socp
begin
    N = 24; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.02;
    delta = zeros(Int64, N);
    ludg = length(UDG);
    deltaV0 = -1;     nr = 50;
    robjvs = zeros(ludg+1,2,nr);
    oobjvs = zeros(ludg+1,2);

    # Setup
    withCascade=true;
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell = getBasicSlaveSocpModel(withCascade);

    @variable(sm, kcvar[1:length(LOADS)], Bin); # kc
    @variable(sm, kgvar[1:length(DG)], Bin); # kg
    kc = zeros(AffExpr,N); kg = zeros(AffExpr,N);
    kc[LOADS] = kcvar; kg[DG] = kgvar;

    addObjectiveDisconnectSocpModel(sm, t, beta, kc, ell);
    sm, deltaineq, deltaineqb, dgAttack = addAttackInequalities(sm, pg, qg, kg, delta, deltaV0);
    sm, deltaeq, deltaeqb, subvolt = addAttackEqualities(sm, v0, delta, deltaV0);
    sm,  basicIneq, basicIneqb = addBasicSocpInequalities(sm, pg, qg, beta, v, t, kc, kg, vpar, P, Q, ell);
    sm, basicSocpEq, basicSocpEqb = addBasicSocpEqualities(sm, v, vpar, p, q, pg, qg, beta, P, Q, kg, ell);

    ludg = length(UDG);
    lloads = length(LOADS);
    if ludg > 0
        @constraint(sm, dgLow, kg[UDG] .>= (vgmin - v)[UDG]);
        @constraint(sm, dgUp, kg[UDG] .>= (v - vgmax)[UDG]);
    end
    @constraint(sm, betaeq, beta[LOADS] .== (1-kc)[LOADS]);

    # println("Before cascade\n",sm);

    @constraint(sm, kgintermediate[i=1:ludg], kg[UDG[i]] >= 0);
    @constraint(sm, loadVoltInterMin[i=1:lloads], kc[LOADS[i]] >= 0);
    @constraint(sm, loadVoltInterMax[i=1:lloads], kc[LOADS[i]] >= 0);
    ### Setup complete

    i = 0; j = 0;
    for i = 1:nr
        println("iter = ",i)
        for deltaV0 in [-1,0]
            rUDG = UDG[randperm(ludg)];
            delta = zeros(Int64,N);

            for M = 0:ludg

                M > 0? delta[rUDG[M]] = 1 : nothing;
                # println(nodes[delta.==1])
                JuMP.setRHS(subvolt, v0nom + v0dis * deltaV0);
                for k = 1:ludg
                    JuMP.setRHS(kgintermediate[k], 0);
                    JuMP.setRHS(deltaineq[k], delta[DG[k]]);
                    JuMP.setRHS(deltaineq[ludg+k], (1-delta[DG[k]])*pgmax[DG[k]]);
                    JuMP.setRHS(deltaineq[2ludg+k], (1-delta[DG[k]])*qgmax[DG[k]]);
                end
                for k = 1:lloads
                    JuMP.setRHS(loadVoltInterMin[k], 0);
                    JuMP.setRHS(loadVoltInterMax[k], 0);
                end

                sstatus = solve(sm);
                if sstatus != :Optimal
                    println("Before cascade\n",sm);
                    println("Model solve status : ", sstatus);
                    return;
                end
                kginv = getvalue(kg);
                vinv = getvalue(v);

                for k = 1:ludg
                    JuMP.setRHS(kgintermediate[k], kginv[DG[k]]);
                end

                for k = 1:lloads
                    JuMP.setRHS(loadVoltInterMin[k], vcmin[LOADS[k]] - vinv[LOADS[k]]);
                    JuMP.setRHS(loadVoltInterMax[k], vinv[LOADS[k]] - vcmax[LOADS[k]]);
                end

                objv = getobjectivevalue(sm);
                # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getOnlineCascadeSocp(delta, deltaV0);
                # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);
                robjvs[M+1,deltaV0+2,i] = objv;
            end
        end
    end

    for i = 1:ludg+1
        for j = 1:2
            oobjvs[i,j] = maximum(robjvs[i,j,:]);
        end
    end

    jldopen(string(mypath,"onlineCascadeSocp",filename,".jld"), "w") do file
        write(file, "oobjvs", oobjvs);
    end
    jldopen(string(mypath,"randomCascadeSocp",filename,".jld"), "w") do file
        write(file, "robjvs", robjvs);
    end
end

# plot benders vs online
begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    afn = string(mypath,"bendersOnlineCascade",filename,".jld");
    bnMs = load(afn)["bnMs"]; bnMs1 = load(afn)["bnMs1"];
    bnobjvs = load(afn)["bnobjvs"]; bnobjvs1 = load(afn)["bnobjvs1"];
    afn = string(mypath,"onlineCascade",filename,".jld");
    oobjvs = load(afn)["oobjvs"]; #odeltas = load(afn)["odeltas"];
    v0dis = 0.02;
    afn = string(mypath,"randomCascade",filename,".jld");
    robjvs = load(afn)["robjvs"];

    clf();
    fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    ax = gca(); grid("on"); ms = 6; lw = 2; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    Ms = collect(0:length(DG));
    plot(100Ms/N, 100-100oobjvs[:,2]/LLCmax, color="red", linewidth=lw, linestyle="-", marker="o", markersize=ms);
    plot(100bnMs1/N, 100-100bnobjvs1/LLCmax, color="green",  linewidth=lw, linestyle="--", marker="s", markersize=ms);
    mylegend = [
    "Worst attack, response (b)",
    "Worst attack, response (c)",
    "Random attacks, response (b)"];
    ylim(20,100);
     grid("on");

    nr = 50;
    println(size(robjvs));
    robjvs = sort(robjvs, 3);
    alpha1 = 0.3; slw = 1; sms = 1;
    for i = 1:nr
        plot(100Ms/N, 100-100robjvs[:,2,i]/LLCmax, color="black", linewidth=slw, linestyle="--",alpha=alpha1);
    end

    legend(mylegend, loc = "lower left", fontsize=12);
    # println(oobjvs);
    # println(sobjvs);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% System resilience",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    savefig(string(mypath, "seqVsOnlineSocp_1_", filename,".pdf"),bbox_inches="tight");

    clf(); fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid(true);
    # Ms = collect(0:N);
    plot(100Ms/N, 100-100oobjvs[:,1]/LLCmax, color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms);
    plot(100bnMs/N, 100-100bnobjvs/LLCmax, color="green",  linewidth=lw,
    linestyle="--", marker="s", markersize=ms);

    for i = 1:nr
        plot(100Ms/N, 100-100robjvs[:,1,i]/LLCmax, color="black", linewidth=slw, linestyle=":",alpha=alpha1);
    end

    mylegend = [
    "Worst attack, response (b)",
    "Worst attack, response (c)",
    "Random attacks, response (b)"];
    ylim(20,100);
    grid(true);

    legend(mylegend, loc = "lower left", fontsize=12);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% System resilience",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    # annotate("Value of !",
	# xy=[x;y,# Arrow tip
	# xytext=[x+dx;y+dy], # Text offset from tip
	# xycoords="data", # Coordinates in in "data" units
	# arrowprops=["facecolor"=>"black"])
    savefig(string(mypath, "seqVsOnlineSocp_2_", filename,".pdf"),bbox_inches="tight");
end

begin
    abc = ones(3,3,3);
    println(abc)
    jldopen(string(mypath,"onlineCascade",filename,".jld"), "w") do file
        write(file, "abc", abc);
    end
    afn = string(mypath,"onlineCascadeTrial",filename,".jld"); # bendersCascadeFilename
    abc1 = load(afn)["abc"];
    println(abc1)
end

# Benders all M Socp
begin
    N = 20; setGlobalParameters(); filename = string("N",N); v0dis = 0.02;

    deltav0 = 0; deltav0min = 0;
    bnMs, bnminDeltav, bnobjvs = getBendersSocpCascadeMinNodesLLC();

    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);

    # jldopen(string(mypath,"bendersCascadeMinNodesSocp",filename,".jld"), "w") do file
    #     write(file, "bnMs", bnMs);
    #     write(file, "bnobjvs", bnobjvs);
    # end
end

# bruteVsBendersSocp1 resilience vs cardinality
begin
    N = 36; filename = string("N",N); setGlobalParameters();
    afn = string(mypath,"bruteForceCascadeSocp",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    afn = string(mypath,"bendersCascadeMinNodesSocp",filename,".jld"); # bendersCascadeFilename
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"];

    fig = figure(L"$\mathcal{R}$ vs $M$",figsize=(8,6));
    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bfMs/N, 100-100bfobjvs/LLCmax, label="Optimal", color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms)
    plot(100bnMs/N, 100-100bnobjvs/LLCmax, label="Benders", color="green", linewidth=lw,
    linestyle="--", marker="s", markersize=ms/2);
    xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    # legend(loc = "upper right", fontsize=15); ## legend position
    legend(loc = "lower left", fontsize=fs); ## legend position
    # new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    ylim(50,100)
    ylim(40,100)
    # ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersSocpPres",filename,".pdf"),bbox_inches="tight");
end


# begin
#     bnMs = 0;
#     bnobjvs = 0;
#     jldopen(string(mypath,"abc.jld"), "w") do file
#         write(file, "bnMs", bnMs);
#         write(file, "bnobjvs", bnobjvs);
#     end
# end

begin
    N = 24; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.02;
    delta = zeros(Int64, N);
    ludg = length(UDG);
    deltaV0 = 0;     nr = 10;

    smmip, pmip, qmip, pgmip, qgmip, betamip, Pmip, Qmip, vmip, tmip, v0mip, kcmip, kgmip, xvarmip, vparmip, ineqmip, ineqbmip, eqmip, eqbmip, deltaeqmip, deltaeqbmip, subvoltmip, deltaineqmip, deltaineqbmip, loadAncestorSuccessorCutsmip, dgAncestorSuccessorCutsmip, ellmip  = getSlaveSocpModelWithCascade(delta, deltaV0, true);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts, ell = getSlaveSocpModelWithCascade(delta, deltaV0, false);

    cumTime1 = 0
    cumTime2 = 0
    modifyinglp = false;
    (startSeq, endSeq) = N == 118 ? (35, 45) : N == 36? (8,18) : (2, 12);

    for i = 1:nr
        println("iter = ",i)
        rUDG = UDG[randperm(ludg)];
        delta = zeros(Int64,N);

        for M = startSeq:endSeq #0:ludg

            delta[rUDG[1:M]] = 1;

            deltaineqmip = modifyBinaryAttackConstraints(delta, deltaV0, deltaineqmip, modifyinglp);

            dgAncestorSuccessorCutsmip = modifyAncestorSuccessorDGCuts(dgAncestorSuccessorCutsmip, delta)
            tic()
            solve(smmip)
            # M == 2? println(smmip) : nothing

            cumTime1 += toq()

            # println(nodes[delta.==1])

            deltaineq = modifyBinaryAttackConstraints(delta, deltaV0, deltaineq, modifyinglp);
            # M == 2? println(sm) : nothing
            tic()
            solve(sm)
            # objv = getobjectivevalue(sm);
            cumTime2 += toq()
            R1 = round(100(1-getobjectivevalue(smmip)/LLCmax),2);
            R2 = round(100(1-getobjectivevalue(sm)/LLCmax),2);
            println("i, M, R1, R2 = ",i, ", ", M,", ", R1 ,", ",R2)
        end
    end
    println("cumTime1 = ", round(cumTime1,2))
    println("cumTime2 = ", round(cumTime2,2))
end

begin
  function getArray(i)
    b = i * ones(10,1)
    # println(i)
    b
  end
  a = SharedArray{Float64}(10,10)

  @parallel for i = 1:10
    a[i,:] = getArray(i)
  end

  println(a)
end
println(CPU_CORES)
begin
  function fme(x)
    # println("x is ",x)
    return (x^2, x^3);
  end

  (ysq, ycu) = @parallel vcat for i= 1:100
    (ysqi,ycui) = fme(i);
  end;

  println(ysq," ", ycu)
end

@time begin
  function fme(x)
    # println("x is ",x)
    return x^2, x^3;
  end
  N = 10^6;
  ysq = Array{Float64}(N); ycu = Array{Float64}(N);
  for i= 1:100
    ysq[i], ycu[i] = fme(i);
  end
end

@time begin
  # @everywhere function fme(x)
  #   # println("x is ",x)
  #   return x^2, x^3;
  # end

  N = 10^6;
  ysq = @parallel vcat for i= 1:N
    ysqi, ycui = fme(i);
  end;

  # println(ysq)
end

begin
  function chut(delta)
    m = Model(solver=GurobiSolver(OutputFlag = 0));
    @variable(m, x[1:N]);
    @variable(m, y[1:N]);
    @constraint(m, xconstr[i=1:N], x[i] >= i * delta[i]);
    @constraint(m, yconstr[i=1:N], x[i] >= i * delta[i]);
    @objective(m, :Min, sum(x[i] + y[i] for i =1:N));
    solve(m);
    xv = getvalue(x); yv = getvalue(y)
    (xv, yv')
  end

  res = chut([1,1,0,0])
end

@time begin

  @everywhere using JuMP, Gurobi
  @everywhere N = 3;
  # @everywhere function chut(delta)
  #   m = Model(solver=GurobiSolver(OutputFlag = 0));
  #   @variable(m, x[1:N]);
  #   @variable(m, y[1:N,1:N]);
  #   @constraint(m, xconstr[i=1:N], x[i] >= i * delta[i]);
  #   @constraint(m, yconstr[i=1:N,j=1:N], y[i,j] >= i * delta[i]);
  #   @objective(m, :Min, sum(x[i] + y[i,j] for i =1:N, j = i ));
  #   solve(m);
  #   xv = getvalue(x); yv = getvalue(y)
  #   println("delta ", delta, " yv ", yv)
  #   (xv, yv')
  # end
  # chut(zeros(N,1))

  for M = 0:N
    attackSets = collect(combinations(collect(1:N),M));
    # res = SharedArray{Tuple{Any, Any}}
    res = @parallel vcat for i = 1:length(attackSets)
      attackSet = attackSets[i];
      delta = zeros(Int64,N); delta[attackSet] = 1;
      chut(delta);
    end
    println(res);
  end
  println("*****")
end

# Brute force all M
begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    levp1 = 10;
    levp1 = length(SDI)+1;
    levp1 = length(DG)+1;
    bfMs = zeros(Int64, levp1);
    bfobjvs = zeros(Float64, levp1); bfsumkcvs = zeros(Int64, levp1);
    bfsumkgvs = zeros(Int64, levp1); bfsumkmvs = zeros(Int64, levp1);
    bfkcvs = zeros(Int64, N, levp1);
    bfkgvs = zeros(Int64, N, levp1); bfDeltas = zeros(Int64, N, levp1);
    bfdeltav0 = zeros(Int64, levp1);

    for M = 1:levp1
        println("M : ", M);
        bestDelta, bestDeltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv,
        pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBruteForceCascadeAttack(M-1);
        # bestDeltaV0 = 0; bestDelta = zeros(Int64,N); bestDelta[DG[end-M+2:end]] = 1;
        # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
        # lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(bestDelta, bestDeltaV0);


        bfobjvs[M] = objv; bfkcvs[:,M] = kcv; bfkgvs[:,M] = kgv;
        bfsumkcvs[M] = sum(kcv); bfsumkgvs[M] = sum(kgv);
        bfDeltas[:,M] = bestDelta; bfMs[M] = M-1; bfdeltav0[M] = bestDeltaV0;
    end

    jldopen(string(mypath,"bruteForceRecovery",filename,".jld"), "w") do file
        write(file, "bfMs", bfMs);
        write(file, "bfobjvs", bfobjvs);
        # write(file, "bfkmvs", bfkmvs);
        write(file, "bfDeltas", bfDeltas);
        write(file, "bfdeltav0", bfdeltav0);
    end

    # println(bfMs);
end

# SOCP Benders with different epsilons
begin
    N = 118; setGlobalParameters(); filename = string("N",N); withCuts = true; v0dis = 0.02;

    deltav0 = 0; deltav0min = 0; bnMs, bnminDeltav, bnobjvs = 0, 0, 0;
    linear = false;
    epsilon = 50; withVariableEpsilon = true; nhighest = 0;

    tic();
    if linear
        bnMs, bnminDeltav, bnobjvs = getBendersCascadeMinNodesLLC(epsilon, withVariableEpsilon, nhighest);

        delta = bnminDeltav[:,1];
        sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts, ell = getSlaveSocpModelWithCascade(delta, deltav0);
        solve(sm);
        bnobjvs[1] = getobjectivevalue(sm);

        for i = 2:length(bnMs)
            delta = bnminDeltav[:,i];
            deltaineq = modifyBinaryAttackConstraints(delta, deltav0, deltaineq);
            solve(sm);
            bnobjvs[i] = getobjectivevalue(sm);
        end
    else
        bnMs, bnminDeltav, bnobjvs = getBendersSocpCascadeMinNodesLLC(epsilon, withVariableEpsilon, nhighest);
    end

    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);
    time = toq();
    println("time = ", time)

    baseFileName = "bendersCascadeMinNodes";
    baseFileName = !linear ? string(baseFileName, "Socp") : baseFileName;
    appendage = withVariableEpsilon ? nhighest : epsilon;

    jldopen(string(mypath,baseFileName,appendage,filename,".jld"), "w") do file
        write(file, string("bnMs",appendage), bnMs);
        write(file, string("bnObjvs",appendage), bnobjvs);
        write(file, "time", time)
    end

    afn = string(mypath,"bruteForceCascadeSocp",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    maxOptimalityGap, avgOptimalityGap = computeOptimalityGap(bfMs, bfobjvs, bnMs, bnobjvs);
    println("maxOptimalityGap = ", round(maxOptimalityGap,2));
    println("avgOptimalityGap = ", round(avgOptimalityGap,2));
end

string("bnMs",epsilon)
# bruteVsBendersSocp1 resilience vs cardinality
begin
    N = 24; filename = string("N",N); setGlobalParameters();
    afn = string(mypath,"bruteForceCascadeSocp",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    afn = string(mypath,"bendersCascadeMinNodesSocp0",filename,".jld");
    bnMs0 = load(afn)["bnMs0"]; bnObjvs0 = load(afn)["bnObjvs0"];
    afn = string(mypath,"bendersCascadeMinNodesSocp1",filename,".jld");
    bnMs1 = load(afn)["bnMs1"]; bnObjvs1 = load(afn)["bnObjvs1"];
    afn = string(mypath,"bendersCascadeMinNodesSocp2",filename,".jld");
    bnMs2 = load(afn)["bnMs2"]; bnObjvs2 = load(afn)["bnObjvs2"];
    afn = string(mypath,"bendersCascadeMinNodesSocp20",filename,".jld");
    bnMs20 = load(afn)["bnMs20"]; bnObjvs20 = load(afn)["bnObjvs20"];
    afn = string(mypath,"bendersCascadeMinNodesSocp50",filename,".jld");
    bnMs50 = load(afn)["bnMs50"]; bnObjvs50 = load(afn)["bnObjvs50"];

    fig = figure(L"$\mathcal{R}$ vs $M$",figsize=(8,6));
    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 1; lw = 1; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bfMs/N, 100-100bfobjvs/LLCmax, label="Optimal", color="red", linewidth=3lw, linestyle="-", marker="o", markersize=12ms)
    plot(100bnMs0/N, 100-100bnObjvs0/LLCmax, label="m = 0", color="blue", linewidth=lw, linestyle=":", marker="s", markersize=6ms);
    plot(100bnMs1/N, 100-100bnObjvs1/LLCmax, label="m = 1", color="black", linewidth=2lw, linestyle="--", marker="o", markersize=6ms);
    plot(100bnMs2/N, 100-100bnObjvs2/LLCmax, label="m = 2", color="blue", linewidth=lw, linestyle="--", marker="s", markersize=6ms);
    plot(100bnMs20/N, 100-100bnObjvs20/LLCmax, label=string("",L"\epsilon", " = 20"), color="cyan", linewidth=lw/2, linestyle="--", marker="s", markersize=6ms);
    if N == 24
        plot(100bnMs50/N, 100-100bnObjvs50/LLCmax, label=string("",L"\epsilon", " = 50"), color="magenta", linewidth=2lw, linestyle=":", marker="s", markersize=6ms);
    end

    afn = string(mypath,"bendersCascadeMinNodes1",filename,".jld");
    bnMs1 = load(afn)["bnMs1"]; bnObjvs1 = load(afn)["bnObjvs1"];
    plot(100bnMs1/N, 100-100bnObjvs1/LLCmax, label="linear, m = 1", color="green", linewidth=2lw, linestyle=":", marker="o", markersize=6ms);
    # afn = string(mypath,"bendersCascadeMinNodes0",filename,".jld");
    # bnMs0 = load(afn)["bnMs0"]; bnObjvs0 = load(afn)["bnObjvs0"];
    # plot(100bnMs0/N, 100-100bnObjvs0/LLCmax, label="linear, m = 0", color="blue", linewidth=2lw, linestyle=":", marker="d", markersize=6ms);

    xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% System resilience",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    # legend(loc = "upper right", fontsize=15); ## legend position
    lfs = N == 24 ? 0.8fs : fs;
    legend(loc = "lower left", fontsize=lfs); ## legend position
    # new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    ylim(50,100)
    ylim(40,100)
    # ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersSocpPres",filename,".pdf"),bbox_inches="tight");
end

# SOCP Benders with different epsilons
begin
    N = 36; setGlobalParameters(); filename = string("N",N); v0dis = 0.02;

    deltav0 = 0; deltav0min = 0;
    epsilon = 50; withVariableEpsilon = true; nhighest = 2;

    tic();
    bnMs, bnminDeltav, bnobjvs = getBendersSocpCascadeMinNodesLLC(epsilon, withVariableEpsilon, nhighest);

    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);
    println("time = ", toq())

    if withVariableEpsilon
        jldopen(string(mypath,"bendersCascadeMinNodesSocp",nhighest,filename,".jld"), "w") do file
            write(file, string("bnMs",nhighest), bnMs);
            write(file, string("bnObjvs",nhighest), bnobjvs);
        end
    else
        jldopen(string(mypath,"bendersCascadeMinNodesSocp",epsilon,filename,".jld"), "w") do file
            write(file, string("bnMs",epsilon), bnMs);
            write(file, string("bnObjvs",epsilon), bnobjvs);
        end
    end
end


begin
    N = 24; filename = string("N",N);
    afn = string(mypath,"bendersCascadeMinNodesSocp50",filename,".jld");
    bnMs50 = load(afn)["bnMs50"]; bnObjvs50 = load(afn)["bnObjvs50"];

    println(bnMs50)
    println(bnObjvs50)
end

# computeOptimalityGap
begin
    N = 24; filename = string("N",N);
    afn = string(mypath,"bruteForceCascadeSocp",filename,".jld"); bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    epsilon = 50; nhighest = 0; linear = false; withVariableEpsilon = true;

    baseFileName = "bendersCascadeMinNodes";
    if !linear
        baseFileName = string(baseFileName, "Socp")
    end
    appendage = withVariableEpsilon ? nhighest : epsilon;
    afn = string(mypath,baseFileName,appendage,filename,".jld");

    bnMs = load(afn)[string("bnMs",appendage)]; bnobjvs = load(afn)[string("bnObjvs",appendage)];
    # time = load(afn)["time"];

    maxOptimalityGap, avgOptimalityGap = computeOptimalityGap(bfMs, bfobjvs, bnMs, bnobjvs);
    println("maxOptimalityGap = ", round(maxOptimalityGap,2));
    println("avgOptimalityGap = ", round(avgOptimalityGap,2));
    println("length(bnMs) = ", length(bnMs))
    # println("time = ", time)
end
