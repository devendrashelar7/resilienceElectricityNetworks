begin
    N = 24; setGlobalParameters(); filename = string("N", N);
    global printPerformance = true;
    global trial = false;
    ResWC = 50;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta));
    println("Resilience: ", 100(1-objv/LLCmax));
end

begin
    N = 24; setGlobalParameters(); filename = string("N", N);

    global trial = true;
    global printPerformance = true;
    ResWC = 65;
    LLCreq = (1 - ResWC/100) * LLCmax;
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta)); println("Resilience: ", 100(1-objv/LLCmax));

end

#microgrid
begin
    N = 24; setGlobalParameters(); filename = string("N", N);
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
