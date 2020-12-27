begin
    N = 118; setGlobalParameters(); filename = string("N", N);
end

begin
    N = 118; setGlobalParameters(); filename = string("N", N);
    global printPerformance = true;
    ResWC = 95; global trial = false;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    # delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getLinearAssistBendersCascade(LLCreq);
    println("LLCreq, LLCmax : ",LLCreq, " ", LLCmax)
    println("Cardinality: ", sum(delta));
    println("% DGs attacked : ", 100sum(delta)/N)
    println("Resilience: ", 100(1-objv/LLCmax));
end

begin
    N = 118; setGlobalParameters(); filename = string("N", N);
    global printPerformance = true;
    ResWC = 0;
    LLCreq = (1-ResWC/100) * LLCmax;
    delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    # delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, prv, qrv, lacv, lovrv, lcontv, lsdv = getLinearAssistBendersCascade(LLCreq);
    println("Cardinality: ", sum(delta));
end


begin
    m = Model(solver=GurobiSolver());
    @variable(m,x[1:3] >=0);
    @variable(m,y[1:3] >=0);
    @constraint(m, y .<= 3);
    # @constraint(m, con, sum(x.^2) <= 1);
    # @constraint(m, con, norm( x[i] for i=1:3 ) <= 3)
    # @constraint(m, con, norm( x) <= 3)
    @constraint(m, con, norm.([x y]) .<= 3)
    # @constraint(m, con, (x[1]^2+y[1]^2)^0.5 <= 3)
    @objective(m, :Min, sum(x));
    solve(m);
    println(m)
    println(getobjectivevalue(m))
    xd = getdual(con);

    println(xd);
    println(isnan.(xd))
    xd[isnan.(xd)] = 0
    println(xd);
    println(norm(getvalue(x)))
end

#microgrid
begin
    N = 118; setGlobalParameters(); filename = string("N", N);
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

begin
  N = 118; setGlobalParameters(); filename = string("N", N);
  deltav0 = -1; deltaf0 = 0; delta = ones(Int64,N);
  pv, qv, pgv, qgv, psynv, qsynv, pvsiv, qvsiv, ppqv, qpqv, betav, vv, v0v, fv, f0v, tvoltv, tfreqv, kcv, kgv, kmv, krsynv, krvsiv, pcv, qcv, objv, Pv, Qv, lovrv, lcontv, lsdv, lmgv = evaluateSlaveModelWithMicrogrid(delta, deltav0, deltaf0);

  println(100objv/LLCmax)

end

begin
    m = Model(solver=GurobiSolver());
    Pv = [2;1;1]; r = 0.01*[1;1;1];
    @variable(m,v[1:3] >=0);
    vpar = ones(AffExpr,3); vpar[2:3] = v[1];
    @constraint(m,vpar - v .== r .* Pv);
    # @objective(m, :Min, 0);
    solve(m);
    println(getvalue(v))

end
