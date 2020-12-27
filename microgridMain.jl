begin
    diversification = false; resourceResponse = true;
    N = 12; setGlobalParameters(); filename = string("N", N);
    # v0dis = 0.03;
    delta = zeros(Int64,N); deltav0 = 0; deltaf0 = 0;
    delta[9:12] = 1
    sm, p, q, pg, qg, psyn, qsyn, pvsi, qvsi, ppq, qpq, beta, P, Q, v, v0, f, f0, tvolt, tfreq, xvar, vpar, fpar, kc, kg, km, krsyn, krvsi, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt, subFreq, deltaineq, deltaineqb, dgAttack = getSlaveModelWithMicrogrid(delta, deltav0, deltaf0);
    # println(sm)
    # sstatus = solve(sm);
    #
    # pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, Pv, Qv, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, lovrv, lcontv, lsdv, lmgv = getMicrogridValues(sm, sstatus, p, q, pg, qg, psyn, qsyn, pvsi, qvsi, ppq, qpq, beta, P, Q, v, v0, f, f0, tvolt, tfreq, kc, kg, km, krsyn, krvsi);

    for i = 1:6
      println(dgAttack[i])
      JuMP.setRHS(dgAttack[i], delta[DG[i]])
      println("after : ", dgAttack[i])
    end
    println(typeof(dgAttack))
    println(typeof(sm))

    # printResults(pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv);
    # println("pev: ", round.(pev,3));

    # println(kmv)

end

begin
    N = 24; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.02;
    delta = zeros(Int64,N); deltav0 = -1; deltaf0 = -1;
    delta[:] = 1

    pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, Pv, Qv, lovrv, lcontv, lsdv, lmgv = evaluateSlaveModelWithMicrogrid(delta, deltav0, deltaf0);

    printMGResults(pv, qv, pgv, qgv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, Pv, Qv);
    # println("pev: ", round.(pev,3));
    # println(kmv)
end

begin # Choosing mq and possibly mp
    N = 12; setGlobalParameters(); filename = string("N", N);
    v0dis = 0.00;
    delta = zeros(Int64,N); deltav0 = -1;
    delta[5:8] = 1

    mq = 0.1;
    semax[RES] = 1srmax;
    pv, pgv, qv, qgv, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv,   qcv, objv, prv, qrv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv = evaluateSlaveModelWithMicrogrid(delta, deltav0);
    printMGResults(pv, pgv, qv, qgv, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, prv, qrv, Pv, Qv);
    # println("pev: ", round.(pev,3));
    # println(kmv)
end

begin
  @everywhere N = 36;
  @everywhere setGlobalParameters();
  @everywhere filename = string("N", N);
  bestDelta, bestdeltav0, pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, Pv, Qv, lovrv, lcontv, lsdv, lmgv = getBruteForceMicrogridAttack(1)
end

begin
  N = 24; setGlobalParameters(); filename = string("N", N);
  v0dis = 0.2; deltav0 = -1; deltaf0 = 0;

  bestDelta, pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, Pv, Qv, lovrv, lcontv, lsdv, lmgv = getBruteForceMicrogridAttackEfficient(11, deltav0, deltaf0);
  println(kmv)
end


begin

  nprocs()
  @everywhere N = 8;
  @everywhere m = Model(solver=GurobiSolver(OutputFlag=0));
  @variable(m, x[1:N]);
  @constraint(m, myAttack, x .>= 0);
  @objective(m, :Min, sum(x));
  @everywhere begin
    using JuMP, Gurobi
    # solve(m)
    # println(myid()," ", getobjectivevalue(m));
  end
  @everywhere function gettemp()

    m = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(m, x[1:N]);
    @constraint(m, myAttack, x .>= 0);
    @objective(m, :Min, sum(x));
    m, myAttack
  end
  ms = Array{JuMP.Model}(nprocs());
  myAttacks = Array{Array{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}},1}}(nprocs());
  for i = 1:nprocs()
    ms[i], myAttacks[i] = gettemp();
  end


  bestobj = @parallel max for i = 1:N

      delta = zeros(Int64,N);
      delta[:] = 0; delta[i] = i;

      pid = myid();
      for j = 1:N
        if j == i
          JuMP.setRHS((myAttacks[pid])[i], i);
        else
          JuMP.setRHS((myAttacks[pid])[j], 0);
        end
      end
      solve(ms[pid]);
      objv = getobjectivevalue(ms[pid])
      println("pid = ", pid, "\ndelta = ",delta,"\nobjv = ",objv);
      objv
  end
  println("best = ",bestobj);
end
# Brute force all M
begin
    N = 24; setGlobalParameters(); filename = string("N", N);
    levp1 = length(DG)+1;
    # levp1 = 2;
    v0dis = 0.2;

    bfMs = zeros(Int64, levp1);
    bfobjvs = zeros(Float64, levp1); bfsumkcvs = zeros(Int64, levp1); bfsumkgvs = zeros(Int64, levp1); bfsumkmvs = zeros(Int64, levp1); bfsumkevs = zeros(Int64, levp1);
    bfkcvs = zeros(Int64, N, levp1); bfkmvs = zeros(Int64, N, levp1); bfkgvs = zeros(Int64, N, levp1); bfkevs = zeros(Int64, N, levp1); bfDeltas = zeros(Int64, N, levp1); bfdeltav0 = zeros(Int64, levp1);

    deltav0 = -1; deltaf0 = 0;

    for M = 1:levp1
      println(M-1);
        bestDelta, pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, Pv, Qv, lovrv, lcontv, lsdv, lmgv = getBruteForceMicrogridAttackEfficient(M-1, deltav0, deltaf0);
        bfobjvs[M] = objv;
        # bfdeltav0[M] = bestDeltaV0;
        bfDeltas[:,M] = bestDelta;
        bfMs[M] = M-1;
    end

    jldopen(string(mypath,"bruteForceMicrogrid",filename,".jld"), "w") do file
        write(file, "bfMs", bfMs);
        write(file, "bfobjvs", bfobjvs);

        write(file, "bfDeltas", bfDeltas);
        # write(file, "bfdeltav0s", bfdeltav0s);
        # write(file, "bfdeltaf0s", bfdeltaf0s);
    end

    # println(bfMs);
end

begin
    N = 12; setGlobalParameters(); #v0dis = 0.03;
    LLCreq = 1025.75
    delta, pv, pgv, qv, qgv, betav, vv, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv, lmgv, legv = getBendersMethodMicrogrid(LLCreq);
end

# bruteVsBenders
begin
  N = 36; setGlobalParameters(); #v0dis = 0.03;
  printPerformance = false;  filename = string("N", N);
  afn = string(mypath,"bruteForceMicrogrid",filename,".jld");
  bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];
  v0dis = 0.2

  bnDeltas = zeros(Int64,N,0);

  # v0dis = 0.0;
  bnobjvs = []; bnMs = [];
  for i = 1:length(bfobjvs)
      LLCreq = bfobjvs[i];
      println(LLCreq);
      # delta, deltav0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
      # objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getLinearAssistBendersCascade(LLCreq);
      delta, pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, Pv, Qv, lovrv, lcontv, lsdv, lmgv = getBendersMethodMicrogrid(LLCreq);
      # if objv >= LLCreq
      bnobjvs = [bnobjvs; objv]; bnMs = [bnMs; sum(delta)];
      bnDeltas = hcat(bnDeltas, delta);
      # end
  end

  clf(); nd = length(bfobjvs); ms = 10; lw = 3;lw1 = 2;
  println(size(bfMs), " ", size(bfobjvs));
  # println(bfdeltav0);
  println(bfobjvs);
  println(bfMs);
  println(bnMs);
  ax = gca(); grid("on"); ms = 12; lw = 3; fs = 15;
  font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

  plot(100bfMs/N, 100bfobjvs/LLCmax, label="Brute Force", color="red", linewidth=lw,
  linestyle="-", marker="o", markersize=ms)
  plot(100bnMs/N, 100bnobjvs/LLCmax, label="Benders", color="green", linewidth=lw,
  linestyle="--", marker="s", markersize=ms/2);
  xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% Post-contingency cost",fontdict=font1); # ylim(0,100);
  setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
  legend(loc = "upper left", fontsize=15); ## legend position
  new_position = [0.15,0.15,0.8,0.8] # Position Method 2
  # ylim(0,50)
  ax[:set_position](new_position)
  savefig(string(mypath,"bruteVsBendersMicrogrid",filename,".pdf"),bbox_inches="tight");
  jldopen(string(mypath,"bendersMicrogrideps",filename,".jld"), "w") do file
      write(file, "bnMs", bnMs);
      write(file, "bnobjvs", bnobjvs);
      write(file, "bnDeltas", bnDeltas);
  end

end

# bruteVsBenders1
begin
    N = 36; filename = string("N",N);
    afn = string(mypath,"bruteForceMicrogrid",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    afn = string(mypath,"bendersMicrogrid",filename,".jld"); # bendersCascadeFilename
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"];

    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bfMs/N, 100-100bfobjvs/LLCmax, label="Optimal", color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms)
    plot(100bnMs/N, 100-100bfobjvs/LLCmax, label="Benders", color="green", linewidth=lw,
    linestyle="--", marker="s", markersize=ms/2);
    xlabel("% Number of nodes attacked",fontdict=font1); ylabel("% Post-contingency cost",fontdict=font1); # ylim(0,100);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    legend(loc = "upper left", fontsize=15); ## legend position
    # new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    # # ylim(0,80)
    # ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersMicrogrid",filename,".pdf"),bbox_inches="tight");
end

# bruteVsBenders2
begin
    N = 24; setGlobalParameters(); filename = string("N",N);
    afn = string(mypath,"bruteForceMicrogrid",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];

    afn = string(mypath,"bendersMicrogrid",filename,".jld"); # bendersCascadeFilename
    bnMs = load(afn)["bnMs"]; bnobjvs = load(afn)["bnobjvs"];

    fig = figure("System resilience vs M",figsize=(5,5));
    clf(); nd = length(bfobjvs);
    ax = gca(); grid("on"); ms = 12; lw = 3; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    plot(100bfMs/N, 100-100bfobjvs/LLCmax, label="Optimal", color="red", linewidth=lw, linestyle="-", marker="o", markersize=ms)
    plot(100bnMs/N, 100-100bfobjvs/LLCmax, label="Benders", color="green", linewidth=lw, linestyle="--", marker="s", markersize=0.7ms);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% System resilience",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);

    legend(loc = "lower left", fontsize=fs); ## legend position
    # new_position = [0.15,0.15,0.8,0.8] # Position Method 2
    N == 24 ? ylim(40,100) : ylim(70,100);
    ax[:set_position](new_position)
    savefig(string(mypath,"bruteVsBendersMicrogrid",filename,".pdf"), bbox_inches="tight");
end

# Benders all M
begin
    N = 24; setGlobalParameters(); filename = string("N",N);


    v0dis = 0.03;
    global deltav0min = -1;
    deltav0 = -1;
    bnMs, bnminDeltav, bnobjvs = getBendersMicrogridMinNodesLLC();
    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);

    global deltav0min = 0;
    deltav0 = 0;
    bnMs1, bnminDeltav, bnobjvs1 = getBendersMicrogridMinNodesLLC();
    bnMs1, bnobjvs1, bnind = getUniqueFirstSortSecondArrays(bnMs1, bnobjvs1);

    jldopen(string(mypath,"bendersOnlineMicrogrid",filename,".jld"), "w") do file
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
    ludg = length(DG);
    deltav0 = -1;
    nr = 20;
    # nr = 5;
    robjvs = zeros(ludg+1,2,nr);
    oobjvs = zeros(ludg+1,2);
    bcobjvs = zeros(ludg+1,2);
    deltav0s = [-1,0];
    i = 0; j = 0;

    # for j = 1:2
    #   deltav0 = deltav0s[j];
    #   for M = 0:ludg
    #     pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);
    #     bcobjvs[M+1,j] = objv;
    #   end
    # end
    for i = 1:nr

      rDG = DG[randperm(ludg)];
      delta = zeros(Int64,N);

      for M = 0:ludg
        M > 0? delta[rDG[M]] = 1 : nothing;

        for j = 1:2
          deltav0 = deltav0s[j];

          if i == 1
            pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);
            bcobjvs[M+1,j] = objv;
          end

          # println(nodes[delta.==1])

          pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, objv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltav0);
          deltav0 == -1 ? println(objv) : nothing;

          if M == 0 #|| i == 2
            objv = WVR * (-deltav0) * v0dis;
          end
          robjvs[M+1,j,i] = objv;
        end
      end
    end

    for i = 1:ludg+1
        for j = 1:2
            oobjvs[i,j] = maximum(robjvs[i,j,:]);
            println(oobjvs[i,:]);
        end
    end

    jldopen(string(mypath,"bendersCascadeMic",filename,".jld"), "w") do file
        write(file, "bcobjvs", bcobjvs);
    end
    jldopen(string(mypath,"onlineCascadeMic",filename,".jld"), "w") do file
        write(file, "oobjvs", oobjvs);
    end
    jldopen(string(mypath,"randomCascadeMic",filename,".jld"), "w") do file
        write(file, "robjvs", robjvs);
    end
end

# plot benders vs online
begin
    N = 24; setGlobalParameters(); filename = string("N",N);
    afn = string(mypath,"bendersOnlineMicrogrid",filename,".jld");
    bnMs = load(afn)["bnMs"];
    bnMs1 = load(afn)["bnMs1"];
    bnobjvs = load(afn)["bnobjvs"]; bnobjvs1 = load(afn)["bnobjvs1"];
    afn = string(mypath,"onlineCascadeMic",filename,".jld");
    oobjvs = load(afn)["oobjvs"]; #odeltas = load(afn)["odeltas"];
    v0dis = 0.02;
    afn = string(mypath,"randomCascadeMic",filename,".jld");
    robjvs = load(afn)["robjvs"];

    afn = string(mypath,"bendersCascadeMic",filename,".jld");
    bcobjvs = load(afn)["bcobjvs"];

    clf();
    fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    ax = gca(); grid("on"); ms = 8; lw = 2; fs = 20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    Ms = collect(0:length(DG));
    plot(100Ms/N, 100-100oobjvs[:,2]/LLCmax, color="red", linewidth=lw, linestyle="-", marker="o", markersize=ms);
    plot(100Ms/N, 100-100bcobjvs[:,2]/LLCmax, color="blue",  linewidth=lw, linestyle=":", marker="x", markersize=ms);
    plot(100bnMs1/N, 100-100bnobjvs1/LLCmax, color="green",  linewidth=lw, linestyle="--", marker="s", markersize=ms);

    # mylegend = [
    # L"\mathcal{R}_{NR},"*" worst attack",
    # L"\mathcal{R}_{Mm},"*" worst attack",
    # L"\mathcal{R}_{MG},"*" worst attack",
    # L"\mathcal{R}_{NR},"*" random attacks"];
    # mylegend = [
    # "Worst attack, response (b)",
    # "Worst attack, response (c)",
    # "Worst attack, response (d)",
    # "Random attacks, response (b)"];
    mylegend = [
    "Worst attack, response (b)",
    "Worst attack, response (c)",
    "Worst attack, response (d)",
    "Random attacks, response (b)"];
    ylim(0,100);
     grid("on");

    nr = 20;
    println(size(robjvs));
    robjvs = sort(robjvs, 3);
    alpha1 = 0.3; slw = 1; sms = 1;
    for i = 1:nr
        plot(100Ms/N, 100-100robjvs[:,2,i]/LLCmax, color="black", linewidth=slw, linestyle=":",alpha=alpha1);
    end

    lfs = 13;
    legend(mylegend, loc = "lower left", fontsize=lfs);
    # println(oobjvs);
    # println(sobjvs);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% System resilience",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    savefig(string(mypath, "seqVsOnlineVsMic_1_", filename,".pdf"),bbox_inches="tight");

    clf(); fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    # Ms = collect(0:N);
    plot(100Ms/N, 100-100oobjvs[:,1]/LLCmax, color="red", linewidth=lw, linestyle="-", marker="o", markersize=ms);
    plot(100Ms/N, 100-100bcobjvs[:,1]/LLCmax, color="blue",  linewidth=lw, linestyle=":", marker="x", markersize=ms);
    plot(100bnMs/N, 100-100bnobjvs/LLCmax, color="green",  linewidth=lw, linestyle="--", marker="s", markersize=ms);


    for i = 1:nr
        plot(100Ms/N, 100-100robjvs[:,1,i]/LLCmax, color="black", linewidth=slw, linestyle=":",alpha=alpha1);
    end


    ylim(0,100);
    grid("on");

    legend(mylegend, loc = "lower left", fontsize=lfs);
    xlabel("% Number of nodes attacked",fontdict=font1);
    ylabel("% System resilience",fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    # annotate("Value of !",
	# xy=[x;y,# Arrow tip
	# xytext=[x+dx;y+dy], # Text offset from tip
	# xycoords="data", # Coordinates in in "data" units
	# arrowprops=["facecolor"=>"black"])
    savefig(string(mypath, "seqVsOnlineVsMic_2_", filename,".pdf"),bbox_inches="tight");
end

# Benders all M
begin

    bnMs, bnminDeltav, bnobjvs, bnsumkcvs, bnsumkgvs,
    bnsumkmvs = getBendersMicrogridMinNodesLLC();
    println(bnMs);
    bnMs, bnobjvs, bnind = getUniqueFirstSortSecondArrays(bnMs, bnobjvs);
    bnsumkcvs = bnsumkcvs[bnind];
    bnsumkgvs = bnsumkgvs[bnind];
    bnsumkmvs = bnsumkmvs[bnind];

    jldopen(string(mypath,"bendersMicrogrid",filename,".jld"), "w") do file
        write(file, "bnMs", bnMs);
        write(file, "bnobjvs", bnobjvs);
        write(file, "bnsumkcvs", bnsumkcvs);
        write(file, "bnsumkgvs", bnsumkgvs);
        write(file, "bnsumkmvs", bnsumkmvs);
    end
end

begin
    global bfMs, bfobjvs, bfsumkcvs, bfsumkgvs, bfsumkcvs;
    global bnMs, bnobjvs, bnsumkcvs, bnsumkgvs, bnsumkcvs;
    jldopen(string(mypath,"bruteForceMicrogrid",filename,".jld"), "r") do file
        bfMs = read(file, "bfMs");
        bfobjvs = read(file, "bfobjvs");
        bfsumkcvs = read(file, "bfsumkcvs");
        bfsumkgvs = read(file, "bfsumkgvs");
        bfsumkmvs = read(file, "bfsumkmvs");
    end
    jldopen(string(mypath,"bendersMicrogrid",filename,".jld"), "r") do file
        bnMs = read(file, "bnMs");
        bnobjvs = read(file, "bnobjvs");
        bnsumkcvs = read(file, "bnsumkcvs");
        bnsumkgvs = read(file, "bnsumkgvs");
        bnsumkmvs = read(file, "bnsumkmvs");
    end
    # jldopen(string(mypath,"onlineCascade",filename,".jld"), "r") do file
    #   oobjvs = read(file, "oobjvs");
    #   sumokcvs = read(file, "sumokcvs");
    #   sumokgvs = read(file, "sumokgvs");
    #   okcvs = read(file, "okcvs");
    #   okgvs = read(file, "okgvs");
    # end
    println("read files");
    clf(); nd = length(bfobjvs);
    println(size(bfMs), " ", size(bfobjvs));
    println(bfMs);
    println(bnMs);
    plot(bfMs, bfobjvs/LLCmax, label=L"$brute$", color="red", linewidth=lw, linestyle="--", marker="x", markersize=ms)
    # plot(Ms, oobjvs/LLCmax, label=L"$cascade$", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    plot(bnMs, bnobjvs, label=L"benders", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Cost"); grid("on");
    # ylim(-50,maximum(oobjvs)+100);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"objMethods",filename,".pdf"));

    clf();
    plot(bfMs, bfsumkcvs, label="brute sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="o", markersize=ms)
    plot(bfMs, bfsumkgvs, label="brute sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    # plot(Ms, bfsumkmvs, label="brute sum(kmv)", color="red", linewidth=lw, linestyle="-", marker="+", markersize=ms)
    # plot(Ms, sumokcvs, label="cascade sum(kcv)", color="blue", linewidth=lw1, linestyle=":", marker="+", markersize=ms)
    # plot(Ms, sumokgvs, label="cascade sum(kgv)", color="green", linewidth=lw, linestyle="-", marker="o", markersize=ms)
    plot(bnMs, bnsumkcvs, label="benders sum(kcv)", color="red", linewidth=lw1, linestyle=":", marker="x", markersize=ms)
    plot(bnMs, bnsumkgvs, label="benders sum(kgv)", color="blue", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    # plot(bnminNodesv, bnsumkmvs, label="benders sum(kmv)", color="green", linewidth=lw, linestyle="-", marker="x", markersize=ms)
    xlabel(L"$M$"); ylabel(L"Number"); grid("on");
    ylim(0,N+4);
    legend(loc = "upper left", fontsize=15); ## legend position
    savefig(string(mypath,"kckgMethods",filename,".pdf"));
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
        bestDelta, bestDeltaV0, pv, pgv, qv, qgv, betav, vv, v0v, tv, kcv, kgv, pcv,
        qcv, objv, prv, qrv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv,
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

begin
    a = [1,2,3,4]; a[end-2:end] = 0; println(a)
end

println(["h "," 2"] .* string(L"y_x"))

begin
    println("What");
    N = 12; attack = [7,8,11,12];
    # N = 6; attack = [2,4,6];

    setGlobalParameters(); M = 4;  filename = string("N", N);
    # bestDelta, bestDeltaV0, pv, pgv, qv, qgv, betav, vv, v0v, tv, kcv, kgv, pcv,
    # qcv, objv, prv, qrv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv,
    # lmgv = getBruteForceMicrogridAttack(M);


    bestDelta = zeros(Int64,N); bestDelta[end-M+1:end] = 1;
    # bestDelta[attack] = 1;

    # println(se2sr);
    deltav0 = -1; se2sr = semax[EG[1]]/srmax; println(se2sr);

    srmaxp = srmax; pc2sr = 10; pr2sr = 1; sd2sr = 0; se2sr = 1;
    modifyLoadAndDGParameters(srmaxp, pc2sr, pr2sr, sd2sr, se2sr);
    sd2srs = linspace(0,5,20); m1 = length(sd2srs);
    v0diss = linspace(0,0.06,2); n1 = length(v0diss);
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
            println(i, " ", sd2sr)
            pv, pgv, qv, qgv, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv,
            qcv, objv, prv, qrv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv,
            legv = evaluateSlaveModelWithMicrogrid(bestDelta, deltav0);
            # pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv,
            # lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(bestDelta, deltav0);
            objvs[i,j], sumkcvs[i,j] = objv, sum(kcv);
            lacvs[i,j], lovrvs[i,j], lcontvs[i,j], lsdvs[i,j] =
            lacv, lovrv, lcontv, lsdv;
            sumkevs[i,j] = sum(kev); sumkmvs[i,j] = sum(kmv); lmgvs[i,j], legvs[i,j] = lmgv, legv;

            # pv, pgv, qv, qgv, betav, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lovrv,
            # lcontv, lsdv = getOnlineCascade(bestDelta, deltav0);
            # oobjvs[i,j], osumkcvs[i,j], osumkmvs[i,j] = objv, sum(kcv), sum(kmv);
        end
    end

    ax = gca(); grid("on"); ms = 10; lw = 2; fs = 12;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    clf();
    fig = figure("pyplot_subplot_mixed",figsize=(10,10))
    plot(sd2srs, objvs/LLCmax, label="Loss (seq)", color="green", linewidth=lw,
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

    clf();
    fig = figure("pyplot_subplot_mixed",figsize=(10,10))
    plot(sd2srs, (lacvs + lovrvs + lcontvs)/LLCmax, label="Loss", color="green", linewidth=lw,
    linestyle="-", marker="x", markersize=ms);
    plot(sd2srs, lsdvs/LLCmax, label="Loss (onine)", color="red", linewidth=lw,
    linestyle="--", marker="o", markersize=ms);
    savefig(string(mypath,"lreg1Vssd2sr",filename,".pdf"));

end

pv, pgv, qv, qgv, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv,
qcv, objv, prv, qrv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv,
legv = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;

println(vv);println(kcv);println(kgv);println(kmv);println(semax);println(pdmax);println(RES);

begin
    clf();
    ax = gca(); grid("on"); ms = 10; lw = 2; fs = 20; lfs =20;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    # fig = figure("pyplot_subplot_mixed",figsize=(5,5))
    plot(sd2srs, objvs/LLCmax, linewidth=lw);
    mylegend = [L"\Delta \mathrm{v} = 0.00",L"\Delta \mathrm{v} = 0.03",L"\Delta \mathrm{v} = ",L"\Delta \mathrm{v} = "];
    myxlabel = "Disturbance to reserves ratio";
    xlabel(myxlabel,fontdict=font1); ylabel(L"L_{SD}",fontdict=font1)
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    # plot(sd2srs, oobjvs, label="Loss (onine)", color="red", linewidth=lw,
    # linestyle="--", marker="o", markersize=ms);
    legend(mylegend, loc = "upper left", fontsize=lfs); ylim(0,0.2);
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
    clf(); plot(sd2srs, lmgvs, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15)
    savefig(string(mypath,"lmgvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, lsdvs/LLCmax, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15)
    xlabel(myxlabel); ylabel(L"L_{SD}")
    savefig(string(mypath,"lsdvVssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, lcontvs + lsdvs, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15)
    savefig(string(mypath,"llcVssd2sr",filename,".pdf"));

    clf();
    # fig = figure("pyplot_subplot_mixed",figsize=(10,10))
    ax = gca();
    plot(sd2srs, (lacvs + lovrvs + lcontvs)/LLCmax, linewidth=lw);
    # xlabel(L"\frac{pd^a}{\overline{sr}}");
    xlabel(myxlabel,fontdict=font1);
    ylabel(L"L^{GC}",fontdict=font1);  ylim(0,0.02);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    legend(mylegend, loc = "upper left", fontsize=lfs);
    savefig(string(mypath,"lreg1Vssd2sr",filename,".pdf"));

    clf(); plot(sd2srs, (lacvs + lovrvs + lcontvs + lsdvs)/LLCmax, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15);
    plot(sd2srs, objvs/LLCmax, linewidth=lw);
    legend(mylegend, loc = "upper left", fontsize=15);
    savefig(string(mypath,"allVssd2sr",filename,".pdf"));


end

println(sumkevs)
println(sumkmvs)

begin
    clf(); plot(sd2srs, objvs, label=["a","b"], color="green", linewidth=lw,
    linestyle="-", marker="x", markersize=ms);
    plot(sd2srs, oobjvs, label=["Loss (onine)"], color="red", linewidth=lw,
    linestyle="--", marker="o", markersize=ms);
    legend(["1","1","1"])
    savefig(string(mypath,"random",filename,".pdf"));
end

begin
    a = [1,2,3]; a = a .* "why";
    println(string(a[1:3]))
end


# computational time for the network size
begin
    Ns = [6, 12, 37, 69]; Ts = [];
    for N in Ns
        setGlobalParameters();
        T = @time
    end
end

begin
    N = 12; setGlobalParameters(); LLCreq = 0.2LLCmax;
    delta, deltav0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
end

begin
    N = 6; setGlobalParameters();
    sm, p, pg, q, qg, pe, qe, beta, v, v0, t, P, Q, kc, kg, km, ke, ineq, ineqb,
    eq, eqb, deltaeq, deltab, deltaMul, subVolt, realCons,
    reacCons = getSlaveModelWithMicrogrid(zeros(Int64,N), 0);
    println(sm)
end

begin
    N = 12; setGlobalParameters(); M = 3; v0dis = 0;
    delta, pv, pgv, qv, qgv, betav, vv, tv, kcv, kgv, kmv, kev, pcv, qcv, objv,
    lacv, lovrv, lcontv, lsdv, lmgv,
    legv = getBendersMicrogridForFixedM(M);
end

begin
    a = [1;2;3]; b = a + 1; c = union(a,b); println(c)
end


begin
    N = 36; setGlobalParameters(); filename = string("N", N);

    fileprefix = "seqVsOnlineVsMic";
    resourceResponse = false;
    resourceResponse? fileprefix *= "Response" : nothing;
    UDG = resourceResponse? setdiff(DG, RES) : DG;
    # srmaxp = srmax; pc2sr = 1; pr2sr = 1; sd2sr = 0; se2sr = 0.6;
    # modifyLoadAndDGParameters(srmaxp, pc2sr, pr2sr, sd2sr, se2sr);

    afn = string(mypath,"bruteForceMicrogrid",filename,".jld"); # bendersCascadeFilename
    bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"]; bfDeltas = load(afn)["bfDeltas"];

    v0dis = 0.02; delta = zeros(Int64, N);
    deltav0 = -1; deltaf0 = 0;
    kmax = length(DG);
    # kmax = 0
    oobjvs = zeros(kmax+1,2);
    sobjvs = zeros(kmax+1,2);
    mobjvs = zeros(kmax+1,2);

    i = 0; j = 0;
    for deltav0 in [0,-1]
        j = deltav0 + 2;
        for M = 0:kmax
            i = M+1;
            delta = zeros(Int64, N); delta[DG[end-M+1:end]] = 1;
            delta = bfDeltas[:,i];
            # delta = copy(bfDeltas[:,i]);
            println("nodes attacked : ", nodes[delta .== 1]);
            pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, oobjv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltav0);

            pv, qv, pgv, qgv, betav, Pv, Qv, vv, tvoltv, kcv, kgv, pcv, qcv, sobjv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltav0);
            pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, mobjv, Pv, Qv, lovrv, lcontv, lsdv, lmgv = evaluateSlaveModelWithMicrogrid(delta,deltav0,deltaf0);

            oobjvs[i,j] = oobjv;
            sobjvs[i,j] = sobjv;
            mobjvs[i,j] = mobjv;

        end
    end

    clf();
    fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");
    ax = gca(); grid("on"); ms = 10; lw = 2; fs = 15;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    Ms = 100collect(0:kmax)/N;
    plot(Ms, 100-100oobjvs[:,2]/LLCmax, color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms);
    plot(Ms, 100-100sobjvs[:,2]/LLCmax, color="green",  linewidth=lw,
    linestyle="--", marker="s", markersize=ms);
    plot(Ms, 100-100mobjvs[:,2]/LLCmax, color="blue",  linewidth=lw,
    linestyle=":", marker="+", markersize=ms);
    mylegend = [
    "No response",
    "Response, without islanding",
    "Response, with islanding"];
    ylim(0,100); grid("on");

    legend(mylegend, loc = "lower left", fontsize=15);
    # println(oobjvs);
    # println(sobjvs);
    xlabel("% Number of nodes attacked", fontdict=font1);
    myylabel = "% System resilience" * L"\mathcal{R}";
    ylabel(myylabel, fontdict=font1);
    savefig(string(mypath, fileprefix, "_1_", filename,".pdf"));

    clf(); fig = figure(L"$L_{VR}$ vs $M$",figsize=(5,5)); ax = gca(); grid("on");

    plot(Ms, 100-100oobjvs[:,1]/LLCmax, color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms);
    plot(Ms, 100-100sobjvs[:,1]/LLCmax, color="green",  linewidth=lw,
    linestyle="--", marker="s", markersize=ms);
    plot(Ms, 100-100mobjvs[:,1]/LLCmax, color="blue",  linewidth=lw,
    linestyle=":", marker="+", markersize=ms);
    # mylegend = [
    # "Uncontrolled " * L"\Delta \mathrm{v} = " * string(v0dis),
    # "Cascade " * L"\Delta \mathrm{v} = " * string(v0dis),
    # "Islanding " * L"\Delta \mathrm{v} = " * string(v0dis)];
    # mylegend = [
    # "Worst attack, no response",
    # "Worst attack, no islanding",
    # "Worst attack, with islanding"];

    ylim(0,100); grid("on");

    legend(mylegend, loc = "lower left", fontsize=12);
    xlabel("% Number of nodes attacked", fontdict=font1);
    ylabel(myylabel,fontdict=font1);
    savefig(string(mypath, fileprefix, "_2_", filename,".pdf"));

    # printResults(pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    # objv, prv, qrv);
    # srmaxp = srmax; pc2sr = 2; pr2sr = 1; sd2sr = 0; se2sr = 1;
    # modifyLoadAndDGParameters(srmaxp, pc2sr, pr2sr, sd2sr, se2sr);
    # println(round(vv,3))
    # println("loads shed : ", nodes[kcv.==1])
    # println("DG shed : ", nodes[kgv.==1])
end

# pge-contingency values
begin
    pgvo = prmax; qgvo = qrmax;
    prvo = 0.9srmax * uv; qrvo = 0.3srmax * uv;
    pvo = pcmax - pgvo - prvo; qvo = qcmax - qrvo - qgvo;
    Pv = SuccMat * pv; Qv = SuccMat * qv;
    # vv = v0nom + deltav0 * v0dis - R * pv - X * qv;

    om = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(om, v[1:N]); vpar = zeros(AffExpr,N); vpar[2:N] = v[par[2:N]]; vpar[1] = v0nom + deltav0 * v0dis;
    @constraint(om, voltDrop, (vpar - v - resistance.*Pv - reactance.*Qv)[noLTC] .== 0);
    if length(LTC) > 0
        @constraint(om, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
    end
    solve(om); vvo = getvalue(v);
end

# diversification
begin
    N = 12; global diversification = true; setGlobalParameters();
    resourceResponse = true;
    filename = string("ResponseN", N); deltav0 = 0; vreg = 0.001;

    # First get original values
    prvo = 0.9srmax * uv; qrvo = 0.4srmax * uv; pgvo = prmax; qgvo = qrmax;
    pvo = pcmax - pgvo - prvo; qvo = qcmax - qrvo - qgvo;
    Pv = SuccMat * pvo; Qv = SuccMat * qvo;
    # vv = v0nom + deltav0 * v0dis - R * pv - X * qv;

    om = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(om, v[1:N]); vpar = zeros(AffExpr,N); vpar[2:N] = v[par[2:N]]; vpar[1] = v0nom + deltav0 * v0dis;
    @constraint(om, voltDrop, (vpar - v - resistance.*Pv - reactance.*Qv)[noLTC] .== 0);
    if length(LTC) > 0
        @constraint(om, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
    end
    solve(om); vvo = getvalue(v); Q1nom = Qv[1];
    println("vvo: ", round.(vvo,3))

    delta[9:12] = 1;

    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta,deltav0);

    clf(); fig = figure("diversification", figsize=(8,12)); ms = 10; fs = 15; lw = 4; lfs = 15;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    new_position = [0.10,0.10,0.8,0.8]; ax = gca(); ax[:set_position](new_position); xlim(0.5,8.5);

    ax = subplot(321); grid("on"); xx = collect(1:N); xx[9:12] -= 4; ylim(0,120);
    title("Uniform allocation", fontdict = font1);
    plot(xx, 100prvo/srmax, linestyle="", marker="x", markersize = ms, color="blue", linewidth=lw);
    plot(xx[1:8], 100prv[1:8]/srmax, linestyle="", marker="<", markersize = ms, color="gray");
    plot(xx[9:12], 100prv[9:12]/srmax, linestyle="", marker=">", markersize = ms, color="red");
    ylabel("Active power \n setpoint (in %)", fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    new_position = [0.1,0.7,0.35,0.25]; ax[:set_position](new_position); xlim(0.5,8.5);

    ax = subplot(323); grid("on"); ylim(0,120);
    plot(xx, 100qrvo/srmax, linestyle="", marker="x", markersize = ms, color="blue", linewidth=lw);
    plot(xx[1:8], 100qrv[1:8]/srmax, linestyle="", marker="<", markersize = ms, color="gray");
    plot(xx[9:12], 100qrv[9:12]/srmax, linestyle="", marker=">", markersize = ms, color="red");
    ylabel("Reactive power \n setpoint (in %)", fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    new_position = [0.1,0.4,0.35,0.25]; ax[:set_position](new_position);  xlim(0.5,8.5);

    ax = subplot(325); grid("on");  yy = collect(0:10);
    plot(yy, 0.96ones(yy), linewidth=lw, color="red");
    plot(xx, vvo, linestyle="", marker="x", markersize = ms, color="blue", linewidth=lw);
    plot(xx[1:8], vv[1:8], linestyle="", marker="<", markersize = ms, color="gray");
    plot(xx[9:12], vv[9:12], linestyle="", marker=">", markersize = ms, color="red");
    ylabel("Voltage (in pu)", fontdict=font1);
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    xlabel("Distance from substation", fontdict=font1);
    new_position = [0.1,0.1,0.35,0.25]; ax[:set_position](new_position); xlim(0.5,8.5);

    # First get original values
    mul = (1 - prvo[1]/srmax) / (mean(xx) - xx[1]);
    prvo = prvo - mul * (xx - mean(xx)) * srmax;
    mul = min(1 - qrvo[1]/srmax, qrvo[1]/srmax) / (mean(xx) - xx[1]);
    qrvo = qrvo + mul * (xx - mean(xx)) * srmax;
    pgvo = prmax; qgvo = qrmax;
    pvo = pcmax - pgvo - prvo; qvo = qcmax - qrvo - qgvo;
    Pv = SuccMat * pvo; Qv = SuccMat * qvo;
    # vv = v0nom + deltav0 * v0dis - R * pv - X * qv;

    om = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(om, v[1:N]); vpar = zeros(AffExpr,N); vpar[2:N] = v[par[2:N]]; vpar[1] = v0nom + deltav0 * v0dis;
    @constraint(om, voltDrop, (vpar - v - resistance.*Pv - reactance.*Qv)[noLTC] .== 0);
    if length(LTC) > 0
        @constraint(om, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
    end
    solve(om); vvo = getvalue(v); Q1nom = Qv[1];
    println("vvo: ", round.(vvo,3))

    pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta,deltav0);

    ax = subplot(322); grid("on"); xx = collect(1:N); xx[9:12] -= 4; ylim(0,120);
    title("Heterogeneous allocation", fontdict = font1);
    plot(xx, 100prvo/srmax, linestyle="", marker="x", markersize = ms, color="blue", linewidth=lw, label=L"$pg^o$");
    plot(xx[1:8], 100prv[1:8]/srmax, linestyle="", marker="<", markersize = ms, color="gray", label=L"$pg^c$, l");
    plot(xx[9:12], 100prv[9:12]/srmax, linestyle="", marker=">", markersize = ms, color="red", label=L"$pg^c$, r"); legend(fontsize=lfs);
    # ax[:legend](bbox_to_anchor=(1.01,1.0));
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    new_position = [0.52,0.7,0.35,0.25]; ax[:set_position](new_position); xlim(0.5,8.5);
    # ylabel("Active power \n setpoint (in %)", fontdict=font1);

    ax = subplot(324); grid("on"); ylim(0,120);
    plot(xx, 100qrvo/srmax, linestyle="", marker="x", markersize = ms, color="blue", linewidth=lw, label=L"$qg^o$");
    plot(xx[1:8], 100qrv[1:8]/srmax, linestyle="", marker="<", markersize = ms, color="gray", label=L"$qg^c$, l");
    plot(xx[9:12], 100qrv[9:12]/srmax, linestyle="", marker=">", markersize = ms, color="red", label=L"$qg^c$, r"); legend(fontsize=lfs);
    # ax[:legend](bbox_to_anchor=(1.01,1.0));
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    new_position = [0.52,0.4,0.35,0.25]; ax[:set_position](new_position); xlim(0.5,8.5);
    # ylabel("Reactive power \n setpoint (in %)", fontdict=font1);

    ax = subplot(326); grid("on");
    plot(yy, 0.96ones(yy), linewidth=lw, color="red");
    plot(xx, vvo, linestyle="", marker="x", markersize = ms, color="blue", linewidth=lw, label=L"$\mathrm{v}^o$");
    plot(xx[1:8], vv[1:8], linestyle="", marker="<", markersize = ms, color="gray", label=L"$\mathrm{v}^c$, l");
    plot(xx[9:12], vv[9:12], linestyle="", marker=">", markersize = ms, color="red", label=L"$\mathrm{v}^c$, r"); legend(fontsize=lfs);
    # ax[:legend](bbox_to_anchor=(1.01,1.0));
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    new_position = [0.52,0.1,0.35,0.25]; ax[:set_position](new_position); xlim(0.5,8.5);
    # ylabel("Voltage \\(in pu)", fontdict=font1);


    xlabel("Distance from substation", fontdict=font1);


    savefig(string(mypath, "diversification", filename,".pdf"));
end

println(Pv[1])

begin

    m = Model(solver=GurobiSolver());

    LOADS = [2,4]; TS = [1,2,3,4];
    @variable(m, xvar[LOADS,TS] >= 0);


    @constraint(m, [i in LOADS, j = 1:4], xvar[i,j] <= 2);
    @objective(m, :Min, sum(xvar[i,ts] for i in LOADS, ts in TS));

    println(m);

    solve(m);

    println(sum(LOADS[j] for j = 4))

end

begin
    # N = 12; setGlobalParameters(); filename = string("N",N);
    global T = 7, TS = collect(1:T), NG = 1;

    deltav0 = 1; delta = zeros(Int64,N); delta[DG] = 1;

    # getRecoveryModel(delta, deltav0, deltaf0)
  betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tvoltv, tfreqv, objv, lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv = evaluateRecoveryModel(delta, deltav0, deltaf0);

end

begin
    clf(); fig = figure(L"$objv$ vs $T$",figsize=(5,5)); ax = gca(); grid("on");
    fs = 15; lw = 2; ms = 10;
    font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

    xdata = zeros(Float64, T+1); ydata = zeros(Float64, T+1);
    for ts in TS
        xdata[ts+1] = ts; ydata[ts+1] = sysPerfv[ts];
    end
    xdata[1] = 0; ydata[1] = 100;

    plot(xdata, ydata, color="red", linewidth=lw,
    linestyle="-", marker="o", markersize=ms);
    mylegend = [
    "SP",
    # "Cascade " * L"\Delta \mathrm{v} = " * string(v0dis),
    # "Islanding " * L"\Delta \mathrm{v} = " * string(v0dis)
    ];
    ylim(0,110); grid("on");

    legend(mylegend, loc = "upper left", fontsize=15);
    xlabel("Time Period", fontdict=font1);
    ylabel("% System Performance",fontdict=font1);
    savefig(string(mypath, "recovery", filename,".pdf"));
end

begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    global T = 7;
    global TS = collect(1:T);
    global NG = 1;
    deltav0 = -1; deltaf0 = 0;
    # deltav0 = 0; deltaf0 = 0;
    delta = zeros(Int64,N); delta[DG] = 1;

    # sm, beta, kc, kg, pg, qg, p, q, km, P, Q, v, v0, tvolt, tfreq, Lovr, Lofr, Lsd, Lcont, Lmg, sysCost, sysPerf, voltDist, freqDist, resourceCons = getRecoveryModel(delta, deltav0, deltaf0);
    iter = 0;

    for v0disiter = 1:2
      if v0disiter == 1
        v0dis = 0.1;
      elseif v0disiter == 2
        v0dis = 0.2;
      end

      clf(); fig = figure(L"$objv$ vs $T$",figsize=(5,5)); ax = gca(); grid("on");
      fs = 18; lw = 2; ms = 10;
      font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

      ydata = zeros(Float64, 4+T+3);
      xdata = [collect(-2:0); collect(0:T+3)];
      ydata[1:3] = 100; #ydata[end-2:end] = 100;

      mycolors = ["red", "green", "blue"];
      mylinestyles = ["-", "--", ":"];
      mymarkers = ["o", "x", "+"];



      for NG = 1:3
        println("NG ",NG)

        # iter += 1;
        # for ts = 1:T
        #     iter == 1 ? println("ts ", ts) : nothing;
        #     if ts <= T-1
        #         JuMP.setRHS(voltDist[ts], v0dis * deltav0);
        #         JuMP.setRHS(freqDist[ts], f0dis * deltaf0);
        #     end
        #     if ts >= 2
        #         JuMP.setRHS(resourceCons[ts], NG);
        #     end
        # end
        #
        # solve(sm);
        # sysPerfv = getvalue(sysPerf);
          # betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tvoltv, tfreqv, objv, lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv = evaluateRecoveryModel(delta, deltav0, deltaf0);
          # lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv =
          sysCostv, sysPerfv = evaluateGreedyRecoveryModel(delta, deltav0, deltaf0);

          for ts in TS
              xdata[ts+4] = ts; ydata[ts+4] = sysPerfv[ts];
          end
          xdata[4] = 0; ydata[4] = sysPerfv[1]; ydata[end-2:end] = sysPerfv[T];

          plot(xdata, ydata, color = mycolors[NG], linewidth=lw, linestyle= mylinestyles[NG], marker = mymarkers[NG], markersize=ms);
      end

      mylegend = [
      L"\mathrm{G} = 1",
      L"\mathrm{G} = 2",
      L"\mathrm{G} = 3"
      ];
      ylim(20,110); grid("on");
      xlim(-2,T+3);

      legend(mylegend, loc = "lower right", fontsize=15);
      xlabel("Time Period", fontdict=font1);
      ylabel("% System Performance",fontdict=font1);
      setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);

      plotname = string(mypath, "recovery2", filename,"_", v0disiter,".pdf");
      # if v0disiter == 1
      #   plotname *= "1.pdf";
      # elseif deltav0 == -1
      #   plotname *= "2.pdf";
      # end
      savefig(plotname,bbox_inches="tight");
    end
end

begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    global T = 7, TS = collect(1:T), NG = 1;
    timeLimit = 7200;

    deltav0 = -1; deltaf0 = 0;
    # deltav0 = 0; deltaf0 = 0;
    delta = zeros(Int64,N); delta[DG] = 1;

    for v0disiter = 1:2
      if v0disiter == 1
        v0dis = 0.1;
      elseif v0disiter == 2
        v0dis = 0.2;
      end

      clf(); fig = figure(L"$objv$ vs $T$",figsize=(5,5)); ax = gca(); grid("on");
      fs = 18; lw = 2; ms = 10;
      font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

      ydata = zeros(Float64, 4+T+3);
      xdata = [collect(-2:0); collect(0:T+3)];
      ydata[1:3] = 100; ydata[end-2:end] = 100;

      mycolors = ["red", "green", "blue", "magenta"];
      mylinestyles = ["-", "--", ":","-."];
      mymarkers = ["o", "x", "+", "s"];

      for NG = 2:3
        println("NG ",NG)
          # betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tvoltv, tfreqv, objv, lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv = evaluateRecoveryModel(delta, deltav0, deltaf0);
          # lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv =
          @time betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tvoltv, tfreqv, objv, lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv = evaluateRecoveryModel(delta, deltav0, deltaf0, timeLimit);
          for ts in TS
              xdata[ts+4] = ts; ydata[ts+4] = sysPerfv[ts];
          end
          xdata[4] = 0; ydata[4] = sysPerfv[1]; ydata[end-2:end] = sysPerfv[T];

          plot(xdata, ydata, color = mycolors[(NG)%2+1], linewidth=lw, linestyle= mylinestyles[NG], marker = mymarkers[(NG)%2+1], markersize=ms);

          @time sysCostv, sysPerfv1 = evaluateGreedyRecoveryModel(delta, deltav0, deltaf0, timeLimit);
          for ts in TS
              xdata[ts+4] = ts; ydata[ts+4] = sysPerfv1[ts];
          end
          xdata[4] = 0; ydata[4] = sysPerfv1[1]; ydata[end-2:end] = sysPerfv1[T];

          plot(xdata, ydata, color = mycolors[NG], linewidth=lw, linestyle= mylinestyles[NG], marker = mymarkers[NG], markersize=ms);

          gc();

      end

      mylegend = [
      "MIP, " * L"\mathrm{G} = 2",
      # L"\mathrm{G} = 1",
      "Greedy, " * L"\mathrm{G} = 2",
      # L"\mathrm{G} = 2"
      "MIP, " * L"\mathrm{G} = 3",
      "Greedy, " * L"\mathrm{G} = 3"
      ];
      ylim(20,110); grid("on");
      xlim(-2,T+3);

      legend(mylegend, loc = "lower right", fontsize=15);
      xlabel("Time Period", fontdict=font1);
      ylabel("% System Performance",fontdict=font1);
      # fs = 15;
      setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
      plotname = string(mypath, "recovery2compare", filename,"_", v0disiter,".pdf");
      # if v0disiter == 1
      #   plotname *= "1.pdf";
      # elseif deltav0 == -1
      #   plotname *= "2.pdf";
      # end
      savefig(plotname,bbox_inches="tight");
    end
end

begin
    N = 36; setGlobalParameters(); filename = string("N",N);
    global T = 7, TS = collect(1:T), NG = 1;
    timeLimit = 600;

    deltav0 = -1; deltaf0 = 0;
    # deltav0 = 0; deltaf0 = 0;
    delta = zeros(Int64,N); delta[DG] = 1;

    for v0disiter = 2:2
      if v0disiter == 1
        v0dis = 0.1;
      elseif v0disiter == 2
        v0dis = 0.2;
      end

      clf(); fig = figure(L"$objv$ vs $T$",figsize=(5,5)); ax = gca(); grid("on");
      fs = 18; lw = 2; ms = 10;
      font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

      ydata = zeros(Float64, 4+T+3);
      xdata = [collect(-2:0); collect(0:T+3)];
      ydata[1:3] = 100; ydata[end-2:end] = 100;

      mycolors = ["red", "green", "blue", "magenta"];
      mylinestyles = ["-", "--", ":","-."];
      mymarkers = ["o", "x", "+", "s"];

      for NG = 3:3
        println("NG ",NG)
          # betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tvoltv, tfreqv, objv, lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv = evaluateRecoveryModel(delta, deltav0, deltaf0);
          # lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv =
          @time betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tvoltv, tfreqv, objv, lovrv, lofrv, lcontv, lsdv, sysCostv, sysPerfv = evaluateRecoveryModel(delta, deltav0, deltaf0, timeLimit);

      end
  end
end


[collect(0:2);collect(2:10)]
begin
    T = 4; TS = collect(1:T); setdiff(TS, [1]);

    LOADS = [2,4];
    delta = ones(Int64, LOADS);
    println(delta)

end

begin
    m = Model(solver=GurobiSolver(OutputFlag=0));
    x = 1; y = 1;
    @variable(m, z);
    @constraint(m, fc1, z <= x);
    @constraint(m, fc2, z <= y);
    @constraint(m, fc3, -z <= -y);
    @constraint(m, ch[i = 1:3], i == i);
    println(ch)

    @objective(m, :Min, -z)
    println(m)
    solve(m);

    dfc1 = getdual(fc1); dfc2 = getdual(fc2);
    println("obj, dfc1, dfc2, dsum: ", getobjectivevalue(m),", ", dfc1, ", ", dfc2, ", ", dfc1+dfc2);


    m = Model(solver=GurobiSolver(OutputFlag=0));
    x = 0; y = 1;
    @variable(m, z); @variable(m, y, Bin);

    @constraint(m, fc1, z == y);
    @constraint(m, fc2, y <= x);
    @objective(m, :Min, -z)
    println(m)
    solve(m);
    zv = getvalue(z); yv = getvalue(y);
    println("z, y: ", zv, ", ", yv);

    dfc1 = getdual(fc1); dfc2 = getdual(fc2);
    println("obj, dfc1, dfc2, dsum: ", getobjectivevalue(m),", ", dfc1, ", ", dfc2, ", ", dfc1+dfc2);

    m = Model(solver=GurobiSolver(OutputFlag=0));
    # x = 1; y = 1;
    @variable(m, z);
    @constraint(m, fc2, z <= x);
    @constraint(m, fc1, z <= yv);
    @constraint(m, fc1b, z >= yv);
    # @constraint(m, fc2, yv <= x);

    @objective(m, :Min, -z)
    println(m)
    solve(m);

    dfc1 = getdual(fc1); dfc2 = getdual(fc2);
    println("obj, dfc1, dfc2, dsum: ", getobjectivevalue(m),", ", dfc1, ", ", dfc2, ", ", dfc1+dfc2);

    m = Model(solver=GurobiSolver(OutputFlag=0));
    x = 0; y = 0;

    @variable(m, z);
    @constraint(m, fc2, z <= x);
    @constraint(m, fc1, z == y);
    # @constraint(m, fc1b, z >= yv);
    # @constraint(m, fc2, yv <= x);

    @objective(m, :Min, -z)
    println(m)
    solve(m);

    dfc1 = getdual(fc1); dfc2 = getdual(fc2);
    println("obj, dfc1, dfc2, dsum: ", getobjectivevalue(m),", ", dfc1, ", ", dfc2, ", ", dfc1+dfc2);
end

begin
  N = 24;
  setGlobalParameters(); #v0dis = 0.03;
  filename = string("N", N);
  afn = string(mypath,"bruteForceMicrogrid",filename,".jld");
  bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];
  bfDeltasM = load(afn)["bfDeltas"];

  afn = string(mypath,"bruteForceCascade",filename,".jld");
  bfMs = load(afn)["bfMs"]; bfobjvs = load(afn)["bfobjvs"];
  bfDeltas = load(afn)["bfDeltas"];

  for i = 1:13
    println("************************")
    println(nodes[bfDeltasM[:,i] .== 1]);
    println(nodes[bfDeltas[:,i] .== 1]);
  end

end

begin
    m = Model(solver=GurobiSolver());
    N=36;
    @variable(m, 1 >= x[1:N] >= 0);
    k = 1000000;
    @objective(m, :Min, 1000 * sum(x) + k * sum(x.*(1-x)));
    @constraint(m, cons[i = 1:N], x[i] >= 0.0012);
    # @constraint(m, cons2[i = 1:N], x[i] <= 0.99999);
    println(m)
    solve(m);
    println("x = ", getvalue(x))
    println("obj = ", getobjectivevalue(m));
    println("cons dual = ", getdual(cons));
    # println("cons2 dual = ", getdual(cons2));
end

begin
    k = 100;
    x = collect(0:0.001:1);
    y = 1000x + k * 1000x .* (1-x)

    myfig = figure("what is this");
    clf();
    plot(x,y);
    j = 11;
    println(x[j], " ",y[j])
    j = 1001;
    println(x[j], " ",y[j])

end

begin
    m = Model(solver=GurobiSolver());
    N=2;
    @variable(m, 1 >= x[1:N] >= 0);
    @objective(m, :Min, 1000 * sum(x));
    @constraint(m, cons[i = 1:N], x[i] <= 2);
    # @constraint(m, cons2[i = 1:N], x[i] <= 0.99999);
    println(m)
    solve(m);
    println("x = ", getvalue(x))
    println("obj = ", getobjectivevalue(m));
    println("cons dual = ", getdual(cons));
    # println("cons2 dual = ", getdual(cons2));
end

begin
    T = 4; yl = [0.5,0.2];
    m = Model(solver=GurobiSolver());
    @variable(m, x[1:3, 1:T], Bin);
    @variable(m, y[1:2, 1:T], Bin);
    @variable(m, z[1:T]);
    @constraint(m, [t = 1:T], z[t] == sum(yl .* y[:,t]));
    @constraint(m, [t = 1:T], y[1,t] <= x[2,t]);
    @constraint(m, [t = 1:T], y[1,t] <= x[1,t]);
    @constraint(m, [t = 1:T], y[2,t] <= x[3,t]);

    @constraint(m, [t = 1:T-1], sum(x[:,t+1]) - sum(x[:,t]) <= 1)
    @constraint(m, sum(x[:,1]) == 0);
    @constraint(m, [t = 1:T-1, i = 1:3], x[i,t+1] >= x[i,t]);
    @objective(m, :Max, sum(z));

    println(m);
    solve(m);

    println("x = ", getvalue(x)');
    println("z = ", getvalue(z));
end
