begin
    N = 69; setGlobalParameters(); filename = string("N", N);
    global printPerformance = true;
    global trial = false;
    ResWC = 80;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta));
    println("Resilience: ", 100(1-objv/LLCmax));
end


begin
  N = 69; setGlobalParameters(); filename = string("N", N);

  pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
  lacv, lovrv, lcontv, lsdv = getEasyBenchmarkValues();

  println("voltage : ", vv);
  println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);

  println("prv : ", prv)

  pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
  lacv, lovrv, lcontv, lsdv = getPrecontingencyValues();

  println("Pre-contingency values")
  println("voltage : ", vv)
  println("objv : ", objv)
  println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);
  println("prv : ", prv);
  println("qrv : ", qrv);

  delta = zeros(Int64,N); deltaV0 = 0;
  pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
  lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltaV0);

  println("Online cascade - no attack")
  println("voltage : ", vv)
  println("objv : ", objv)
  println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);
  println("prv : ", prv);
  println("qrv : ", qrv);

  delta[end-5+1:end] = 1;
  v0dis = 0.03; deltaV0 = 1;
  pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
  lacv, lovrv, lcontv, lsdv = getOnlineCascade(delta, deltaV0);

  println("Online cascade - all attack")
  println("voltage : ", vv)
  println("objv : ", objv)
  println("lacv, lovrv, lcontv, lsdv : ", lacv," ", lovrv," ", lcontv," ", lsdv);
  println("prv : ", prv);
  println("qrv : ", qrv);
end

println(LLCmax)

begin
  v0dis = 0.001;
  bnMs, bnminDeltav, bnobjvs, bnsumkcvs, bnsumkgvs,
  deltaV0s = getBendersCascadeMinNodesLLC();

  clf(); nd = length(bfobjvs); ms = 8; lw = 3;lw1 = 2;
  println(bnMs);

  oMsIndex = getIndexUniqueLargestIndex(bnMs);
  oMs = bnMs[oMsIndex]; odeltaV0s = deltaV0s[oMsIndex];
  ooDeltas = bnminDeltav[:,oMsIndex]; oobjvs = [];
  for i = 1:length(oMs)
    delta = ooDeltas[:,i]; deltaV0 = odeltaV0s[i];
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, oobjv, pgv, qgv,
    lacv, lovrv, lcontv, olsdv = getOnlineCascade(delta, deltaV0);
    oobjvs = [oobjvs; olsdv];
  end

  fig = figure(L"$objv$ vs $M$",figsize=(5,5));
  ax = gca(); grid("on"); ms = 10; lw = 2; fs = 12;
  font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);

  plot(bnMs, 100bnobjvs, label=L"benders", color="green", linewidth=lw,
  linestyle="--", marker="s", markersize=ms);
  plot(oMs, 100oobjvs/LLCmax, label=L"uncontrolled", color="red", linewidth=lw,
  linestyle="--", marker="o", markersize=ms);
  xlabel(L"\mathrm{k}",fontdict=font1); ylabel("% Load Shedding",fontdict=font1);

  legend(loc = "upper left", fontsize=15); ## legend position
  savefig(string(mypath,"bendersVsOnline",filename,".pdf"));

end


begin

end
