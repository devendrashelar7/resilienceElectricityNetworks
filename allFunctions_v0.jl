# workspace();
# Pkg.add("JuMP")
# Pkg.add("Combinatorics");
# Pkg.add("Images");
# Pkg.add("JLD");
# Pkg.add("Gallium")
using PyPlot, JuMP, Gurobi, Mosek, Combinatorics, Images, JLD, CSV, DataFrames, Gallium

genv = Gurobi.Env()

global resourceResponse = false;
global deltaV0min = 0;
global printPerformance = true;
global trial = false;
global LTCSetting = 1.05;
global diversification = false;
global vreg = 0;
global P1nom = 0;
global Q1nom = 0;

global hasDroop = false;
homedir()
pwd()
cd("/Users/devendrashelar/Dropbox (MIT)/DNCS/papers/2019_TCNS_resubmission/code")
pwd()

mypath = "../results/";

function getLeaves(par)
    N = length(par); isLeaf = ones(Int64,N);
    for i = 1:length(par)
        par[i] != 0 ? isLeaf[par[i]] = 0 : continue;
    end
    leaves = find(isLeaf .== 1);
    leaves
end

function getPath(par)
    N = length(par);
    path = zeros(Int64,N);
    for i=1:N
        j = i;
        while j!=0
            path[i] += 1;
            j = par[j];
        end
    end
    path
end

function getChildrenMatrix(par)
    N = length(par);
    CHILD = zeros(Int64,N,N);
    for node=1:N
        parent = par[node];
        if parent != 0
            CHILD[parent, node] = 1;
        end
    end

    CHILD
end

function getSuccessorMatrix()
    SuccMat = zeros(Int64,N,N);
    for i = 1:N
        cur = i;
        while cur != 0
            SuccMat[cur, i] = 1;
            cur = par[cur];
        end
    end

    SuccMat
end

function getChildrenSets()
    CHILDSets = Array{Array{Int64}}(N);
    for i = 1:N
        CHILDSets[i] = nodes[CHILD[i,:] .== 1]
    end

    CHILDSets
end

# CHILDSets = getChildrenSets()
# println(CHILDSets)

# compute R and X
function getCommonPathVariables(N,par,r,x)
    # println("r, x: ", r, "\n",x)
    R = zeros(Float64,N,N); X = zeros(Float64,N,N); commonPath = zeros(Int64,N,N);
    for i=1:N
        for j=1:N
            i1 = i;
            j1 = j;
            lca = 0;
            #       println("going in ",i, " ", j);
            while lca==0 # least common ancestor not found
                #         println(i1, " ", j1);
                while j1!=0 # j1 has not reached the root
                    if j1 == i1 # if a common ancestor found, note and break
                        lca = j1;
                        break;
                    end
                    j1 = par[j1]; # go up path of j
                end

                if lca!=0
                    break;
                end
                i1 = par[i1]; j1 = j;
                if i1 == 0 && j1==0
                    println("Error in the tree. Check variable par.");
                    break;
                end
            end
            if lca!=0
                k = lca;
                #         println("lca ", lca);
                while k!=0 #populate values of R[i,j] and X[i,j]
                    commonPath[i,j] += 1;
                    R[i,j] += r[k]; X[i,j] += x[k];
                    k = par[k];
                end
            end
        end
    end
    commonPath, R, X
end

function getMasterModel()
    mm = Model(solver=GurobiSolver(OutputFlag = 0, genv));

    @variable(mm, dvar[1:length(UDG)], Bin);
    deltaVar = zeros(AffExpr,N);
    deltaVar[UDG] = dvar;

    @objective(mm, :Min, sum(deltaVar));
    # @constraint(mm, noDelta[i = 1:length(noSDI)], delta[noSDI[i]] == 0);
    # println(mm);

    mm, deltaVar
end

function addVariablesToBasicSlaveModel(sm)
    @variable(sm, p[1:N]); @variable(sm, q[1:N]);
    @variable(sm, prvar[1:length(DG)]); @variable(sm, qrvar[1:length(DG)]);

    @variable(sm, pgvar[1:length(RES)]); @variable(sm, qgvar[1:length(RES)]);
    @variable(sm, betavar[1:length(LOADS)]);

    #voltages are squared
    @variable(sm, v[1:N]); @variable(sm, v0);
    # @variable(sm, t[1:N]);
    @variable(sm, tvar);
    @variable(sm, P[1:N]); @variable(sm, Q[1:N]);


    pr = zeros(AffExpr,N); qr = zeros(AffExpr,N); beta = zeros(AffExpr,N); t = zeros(AffExpr,N); pg = zeros(AffExpr,N); qg = zeros(AffExpr,N);

    pr[DG] = prvar; qr[DG] = qrvar; beta[LOADS] = betavar; vpar = zeros(AffExpr,N); t[:] = tvar; pg[RES] = pgvar; qg[RES] = qgvar;
    for i = 1:N
        vpar[i] = par[i] == 0 ? v0 : v[par[i]];
    end

    xvar = [p; pr; q; qr; beta; P; Q; v; v0; t];

    p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar
end

function getBasicSlaveModel()
    sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));

    p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = addVariablesToBasicSlaveModel(sm);

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar
end

function getBasicSlaveSocpModel(withCascade=true)
    sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));
    if withCascade == false
        sm = Model(solver=MosekSolver(QUIET=true));
    end
    # sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));
    p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = addVariablesToBasicSlaveModel(sm);

    @variable(sm, ell[1:N]);
    xvar = vcat(xvar, ell);

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell
end

function getVariableValues(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v,  t, kc, kg, vpar)
    vv = zeros(Float64, length(v)); pv = zeros(Float64, length(p));
    prv = zeros(Float64, length(pr)); qv = zeros(Float64, length(q));
    qrv = zeros(Float64, length(qr)); betav = zeros(Float64, length(beta));
    tv = 0; kcv = zeros(Bool, length(kc)); kgv = zeros(Bool, length(kg));
    pcv = zeros(Float64, length(beta)); qcv = zeros(Float64, length(beta));
    pgv = zeros(Float64, length(pgmax)); qgv = zeros(Float64, length(qgmax));
    Pv = zeros(Float64, N); Qv = zeros(Float64, N);
    objv = -infty; lacv, lovrv, lcontv, lsdv = 0,0,0,0;

    if sstatus != :Optimal
        println("Slave model solve status : ", sstatus);
        return pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv
    end

    vv = getvalue(v); pv = getvalue(p); qv = getvalue(q);
    tv = getvalue(t);
    betav = getvalue(beta); pcv = pcmax .* betav; qcv = qcmax .* betav;
    prv = getvalue(pr); qrv = getvalue(qr);
    Pv = getvalue(P); Qv = getvalue(Q);
    pprv = 100getvalue(pr)/srmax; pqrv = 100getvalue(qr)/srmax;
    objv = getobjectivevalue(sm);
    # println("kc : ", getvalue(kc), ", kg :", getvalue(kg))
    kcv = round.(Int64, getvalue(kc)); kgv = round.(Int64, getvalue(kg));
    # pgv = (1-kgv).*pgmax; qgv = (1-kgv).*qgmax;
    pgv = getvalue(pg); qgv = getvalue(qg);

    vparv = getvalue(vpar);
    ellv = (Pv.^2+Qv.^2)./vparv;
    lacv = sum(WAC) * sum(resistance.*ellv);
    lovrv = WVR * sum(tv) / N;
    lcontv = sum(((1-betav).*WLC)[LOADS]); lsdv = sum((WSD-WLC).*kcv);

    # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv = round.([pv, prv[DG], qv, qrv[DG], betav, Pv, Qv, vv, tv],3);
    # pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = round.([pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv],3);


    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getVariableValuesSocp(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v, t, kc, kg, ell, vpar)
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getVariableValues(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v,  t, kc, kg, vpar);

    ellv = zeros(Float64, length(ell));
    if sstatus == :Optimal
        ellv = getvalue(ell)
    end
    lacv = sum(WAC) * sum(resistance.*ellv);

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function addObjectiveDisconnectModel(sm, t, beta, kcval, pr, P)
    @objective(sm, :Min,
    # sum(WAC .* pr) +
    # sum(WAC) * P[1] +
    WVR * sum(t) / N +
    sum(WLC[LOADS].*(1-beta)[LOADS]) +
    sum((WSD-WLC).*kcval));
end

function addObjectiveDisconnectSocpModel(sm, t, beta, kcval, pr, P, ell)
    @objective(sm, :Min,
    # sum(WAC .* pr) +
    # sum(WAC) * P[1] +
    WVR * sum(t) / N +
    sum(WLC[LOADS].*(1-beta)[LOADS]) +
    sum((WSD-WLC).*kcval) +
    sum(WAC)*sum(resistance.*ell));
end


function addBasicInequalities(sm, pr, qr, pg, qg, beta, v, t, kcval, kgval)
    @constraint(sm, betalb, beta[LOADS] .>= (1-kcval)[LOADS] .* betamin[LOADS]);
    @constraint(sm, betaub, beta[LOADS] .<= (1-kcval)[LOADS]);

    @constraint(sm, lovrMin, t .>= v0nom * uv - v);
    @constraint(sm, lovrMax, t .>= v - v0nom * uv);

    basicIneq = [betalb; betaub; lovrMin; lovrMax];
    basicIneqb = [((1-kcval) .* betamin)[LOADS]; (1-kcval)[LOADS]; v0nom * uv; -v0nom * uv];

    sm, basicIneq, basicIneqb
end

function addBasicEqualities(sm, v, vpar, p, q, pr, qr, pg, qg, beta, P, Q, kgval)
    @constraint(sm, voltDrop, (vpar - v - 2resistance.*P - 2reactance.*Q)[noLTC] .== 0)
    @constraint(sm, realFlow, P .== SuccMat * p);
    @constraint(sm, reacFlow, Q .== SuccMat * q);
    # @constraint(sm, realCons, p - beta.*pcmax + pr + pg .== 0);
    # @constraint(sm, reacCons, q - beta.*qcmax + qr + qg .== 0);
    @constraint(sm, realCons, p - beta.*pcmax + pr .== 0);
    @constraint(sm, reacCons, q - beta.*qcmax + qr .== 0);

    basicEq = [voltDrop; realFlow; reacFlow; realCons; reacCons];
    basicEqb = [ov[noLTC]; ov; ov; ov; ov];
    if length(LTC) > 0
        @constraint(sm, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
        basicEq = vcat(basicEq, voltDropLTC);
        basicEqb = vcat(basicEqb, ov[LTC]);
    end

    # UDG = resourceResponse ? setdiff(DG,RES) : DG;
    # if length(UDG) > 0
    #     @constraint(sm, UDG1, pr[UDG] .== (1-kgval[UDG]).* pgmax[UDG]);
    #     @constraint(sm, UDG2, qr[UDG] .== (1-kgval[UDG]).* qgmax[UDG]);
    #     basicEq = vcat(basicEq, [UDG1; UDG2]);
    #     basicEqb = vcat(basicEqb, [(1-kgval[UDG]).* pgmax[UDG]; (1-kgval[UDG]).* qgmax[UDG]]);
    # end

    sm, basicEq, basicEqb
end

function addBasicSocpEqualities(sm, v, vpar, p, q, pr, qr, pg, qg, beta, P, Q, kgval, ell)
    @constraint(sm, voltDrop, (vpar - v - 2resistance.*P - 2reactance.*Q + (resistance.^2+reactance.^2).*ell)[noLTC] .== 0)
    @constraint(sm, realFlow, P .== CHILD * P + p + resistance.*ell);
    @constraint(sm, reacFlow, Q .== CHILD * Q + q + reactance.*ell);
    # @constraint(sm, realCons, p - beta.*pcmax + pr + pg .== 0);
    # @constraint(sm, reacCons, q - beta.*qcmax + qr + qg .== 0);
    @constraint(sm, realCons, p - beta.*pcmax + pr .== 0);
    @constraint(sm, reacCons, q - beta.*qcmax + qr .== 0);

    basicSocpEq = [voltDrop; realFlow; reacFlow; realCons; reacCons];
    basicSocpEqb = [ov[noLTC]; ov; ov; ov; ov];
    if length(LTC) > 0
        @constraint(sm, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
        basicSocpEq = vcat(basicSocpEq, voltDropLTC);
        basicSocpEqb = vcat(basicSocpEqb, ov[LTC]);
    end

    sm, basicSocpEq, basicSocpEqb
end

function addBasicSocpInequalities(sm, pr, qr, pg, qg, beta, v, t, kcval, kgval, vpar, P, Q, ell)
    sm, basicIneq, basicIneqb = addBasicInequalities(sm, pr, qr, pg, qg, beta, v, t, kcval, kgval);

    @constraint(sm, socpIneq[i=1:N], norm([vpar[i], ell[i], sqrt(2)*P[i], sqrt(2)*Q[i]]) <= (vpar[i] + ell[i]));

    socpIneqb = ov;

    sm, socpIneq, socpIneqb, basicIneq, basicIneqb
end

function addAttackInequalities(sm, pr, qr, kgval, delta, deltaV0)

    # @constraint(sm, dgAttack, -kgval[UDG] .<= -delta[UDG]);
    # @variable(sm, kgconvar[1:length(UDG)]); kgcon = zeros(AffExpr,N); kgcon[UDG] = kgconvar;

    # @constraint(sm, dgAttack0, kgcon[UDG] .== delta[UDG]);
    @constraint(sm, dgAttack1, kgval[UDG] .>= delta[UDG]);
    global isCascadeMIP
    if isCascadeMIP
        global eta = 0;
    else
        global eta = 10epsilon;
    end

    @constraint(sm, dgAttack2a, pr[UDG] .<= (1 - delta[UDG]) .* pgmax[UDG]);
    @constraint(sm, dgAttack2b, qr[UDG] .<= (1 - delta[UDG]) .* qgmax[UDG]);
    @constraint(sm, dgAttack3a, pr[UDG] .<= (1 - kgval[UDG]) .* pgmax[UDG]);
    @constraint(sm, dgAttack3b, pr[UDG] .>= (1 - kgval[UDG]) .* pgmax[UDG]);
    @constraint(sm, dgAttack4a, qr[UDG] .<= (1 - kgval[UDG]) .* qgmax[UDG]);
    @constraint(sm, dgAttack4b, qr[UDG] .>= (1 - kgval[UDG]) .* qgmax[UDG]);

    dgAttack = [dgAttack1; dgAttack2a; dgAttack2b; dgAttack3a; dgAttack3b; dgAttack4a; dgAttack4b];
    # dgAttack = dgAttack1;
    deltaineq = dgAttack;
    # deltaineqb = ov[UDG];
    # deltaineqb = [ov[UDG]; pgmax[UDG]; qgmax[UDG]];
    deltaineqb = [ov[UDG]; pgmax[UDG]; qgmax[UDG]; (1 - kgval[UDG]) .* pgmax[UDG]; (1 - kgval[UDG]) .* pgmax[UDG]; (1 - kgval[UDG]) .* qgmax[UDG]; (1 - kgval[UDG]) .* qgmax[UDG]];

    sm, deltaineq, deltaineqb, dgAttack
end

function addAttackEqualities(sm, v0, vpar, p, q, pr, qr, beta, P, Q, kgval, delta, deltaV0)
    # println("addAttackEqualities: v0nom, deltaV0 : ", v0nom, ", ", deltaV0);
    # @constraint(sm, subvolt, v0 == v0nom + v0dis * deltaV0 - vreg * (Q[1] - Q1nom));
    @constraint(sm, subvolt, v0 == v0nom + v0dis * deltaV0);

    deltaeq = [subvolt];
    # deltaeqb = [v0nom + vreg * Q1nom];
    deltaeqb = [v0nom];

    sm, deltaeq, deltaeqb, subvolt
end

# function addAttackConstraints(sm, v0, vpar, p, q, pr, qr, beta, P, Q, kgval, delta, deltaV0)
#
# end


function addDisconnectivityConstraints(sm, kc, kg, v, lndgDisconnect = false)
    # println("lengths: ", length(kc)," ",length(vcmin)," ",length(v))
    @constraint(sm, loadLow, kc[LOADS] .>= (vcmin - v)[LOADS]);
    @constraint(sm, loadUp, kc[LOADS] .>= (v - vcmax)[LOADS]);

    disineq = [loadLow; loadUp]; disb = [vcmin[LOADS]; -vcmax[LOADS]];
    if length(UDG) > 0
        @constraint(sm, dgLow, kg[UDG] .>= (vgmin - v)[UDG]);
        @constraint(sm, dgUp, kg[UDG] .>= (v - vgmax)[UDG]);
        disineq = vcat(disineq, dgLow, dgUp); disb = vcat(disb, vgmin[UDG], -vgmax[UDG]);
    end

    LnDG = intersect(LOADS, DG);
    if lndgDisconnect && length(LnDG) > 0
        vcminLessVgmin = vcmin[LnDG] .<= vgmin[LnDG];
        length(LnDG) > 0 ? @constraint(sm, vcminLessVgmin .* (kc - kg)[LnDG] .<= 0): nothing;
    end

    sm, disineq, disb
end


function addCascadeConstraints(sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, kcv, kgv, delta, deltaV0)
    # Whenever you add constraints, add the update the list of constraints and right hand side ineq and ineqb

    sm, deltaineq, deltaineqb, dgAttack = addAttackInequalities(sm, pr, qr, kgv, delta, deltaV0);
    sm, deltaeq, deltaeqb, subvolt = addAttackEqualities(sm, v0, vpar, p, q, pr, qr, beta, P, Q, kgv, delta, deltaV0);
    sm, basicIneq, basicIneqb = addBasicInequalities(sm, pr, qr, pg, qg, beta, v, t, kcv, kgv);
    sm, basicEq, basicEqb = addBasicEqualities(sm, v, vpar, p, q, pr, qr, pg, qg, beta, P, Q, kgv);
    sm, disineq, disb = addDisconnectivityConstraints(sm, kcv, kgv, v, false);

    ineq = vcat(basicIneq, deltaineq, disineq); ineqb = vcat(basicIneqb, deltaineqb, disb);
    eq = [basicEq; deltaeq]; eqb = [basicEqb; deltaeqb];

    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack
end

function addSocpCascadeConstraints(sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell, kcv, kgv, delta, deltaV0)
    # Whenever you add constraints, add the update the list of constraints and right hand side ineq and ineqb

    sm, deltaineq, deltaineqb, dgAttack = addAttackInequalities(sm, pr, qr, kgv, delta, deltaV0);
    sm, deltaeq, deltaeqb, subvolt = addAttackEqualities(sm, v0, vpar, p, q, pr, qr, beta, P, Q, kgv, delta, deltaV0);
    sm, socpIneq, socpIneqb, basicIneq, basicIneqb = addBasicSocpInequalities(sm, pr, qr, pg, qg, beta, v, t, kcv, kgv, vpar, P, Q, ell);
    sm, basicSocpEq, basicSocpEqb = addBasicSocpEqualities(sm, v, vpar, p, q, pr, qr, pg, qg, beta, P, Q, kgv, ell);
    sm, disineq, disb = addDisconnectivityConstraints(sm, kcv, kgv, v, false);

    ineq = vcat(basicIneq, deltaineq, disineq); ineqb = vcat(basicIneqb, deltaineqb, disb);
    eq = [basicSocpEq; deltaeq]; eqb = [basicSocpEqb; deltaeqb];

    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack
end

# Given: attack vector, connectivity of LOADS and DERs that remain unchanged.
# Output: slave model which just an LP.
function getSlaveModelWithoutCascade(kcv, kgv, delta, deltaV0)
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = getBasicSlaveModel();

    addObjectiveDisconnectModel(sm, t, beta, kcv, pr, P);
    global isCascadeMIP = false;
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack = addCascadeConstraints(sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, kcv, kgv, delta, deltaV0);

    # println(sm);
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack
end

function getSlaveSocpModelWithoutCascade(kcv, kgv, delta, deltaV0)
    withCascade = false;
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell = getBasicSlaveSocpModel(withCascade);

    addObjectiveDisconnectSocpModel(sm, t, beta, kcv, pr, P, ell);

    global isCascadeMIP = false;
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack = addSocpCascadeConstraints(sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell, kcv, kgv, delta, deltaV0);

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kcv, kgv, xvar, vpar, ell, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack
end


function getSlaveModelWithCascade(delta, deltaV0)
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = getBasicSlaveModel();

    @variable(sm, kcvar[1:length(LOADS)], Bin); # kc
    @variable(sm, kgvar[1:length(DG)], Bin); # kg
    kc = zeros(AffExpr,N); kg = zeros(AffExpr,N);
    kc[LOADS] = kcvar; kg[DG] = kgvar;

    addObjectiveDisconnectModel(sm, t, beta, kc, pr, P);
    global isCascadeMIP = true;
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack  = addCascadeConstraints(sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, kc, kg, delta, deltaV0);
    # sm = addDisconnectivityConstraints(sm, kc, kg, v);

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack
end

function getSlaveSocpModelWithCascade(delta, deltaV0)
    withCascade=true;
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell = getBasicSlaveSocpModel(withCascade;

    @variable(sm, kcvar[1:length(LOADS)], Bin); # kc
    @variable(sm, kgvar[1:length(DG)], Bin); # kg
    kc = zeros(AffExpr,N); kg = zeros(AffExpr,N);
    kc[LOADS] = kcvar; kg[DG] = kgvar;

    addObjectiveDisconnectSocpModel(sm, t, beta, kc, pr, P, ell);

    global isCascadeMIP = true;
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack = addSocpCascadeConstraints(sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell, kc, kg, delta, deltaV0);
    # sm = addDisconnectivityConstraints(sm, kc, kg, v);

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack, ell
end

function evaluateSlaveModelWithCascade(delta, deltaV0)
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithCascade(delta, deltaV0);
    # println(sm);
    sstatus = solve(sm);

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v, t, kc, kg, vpar);
    # printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function evaluateSlaveSocpModelWithCascade(delta, deltaV0)
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack, ell = getSlaveSocpModelWithCascade(delta, deltaV0);

    # println(sm);
    sstatus = solve(sm);

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValuesSocp(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v, t, kc, kg, ell, vpar)
    # printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getEasyBenchmarkValues()
    betav = ones(Float64,N); kcv = zeros(Int64,N); kgv = copy(kcv);
    pcv = pcmax; qcv = qcmax;
    pgv = pgmax; qgv = qgmax;
    prv = srmax; qrv = srmax/3;

    pv = pcv - prv;
    qv = qcv - qrv;
    vv = v0nom - R * pv - X * qv;
    tv = sum(abs(vv-v0nom));
    Pv = SuccMat * pv; Qv = SuccMat * qv;

    lacv = 0; lovrv = WVR * tv; lcontv = sum(WLC .* (1-kcv) .* (1-betav));
    lsdv = sum((WSD - WLC) .* kcv);
    objv = lacv + lovrv + lcontv + lsdv;

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv
end

# The goal is to ensure that the voltages are within the bounds, but the loss
# of load control and load shedding is minimized as much as possible.
function getPrecontingencyValues()
    delta = zeros(Int64,N); deltaV0 = 0;
    # global betamin = uv;
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    # printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    # objv, pgv, qgv)

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv
end

function getPrecontingencyValuesSecondMethod()
    delta = zeros(Int64,N); deltaV0 = 0;

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv
end

function getOnlineCascade(delta, deltaV0)

    noAttackDelta = zeros(Int64,N); noDeltaV0 = 0
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(noAttackDelta, noDeltaV0);
    # println("vv: ",vv);
    # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    # lacv, lovrv, lcontv, lsdv = getPrecontingencyValues();
    # println("Loads disconnected: ", nodes[kcv.==1]);
    # println("DGs disconnected: ", nodes[kgv.==1]);
    # println("delta: ",delta)

    kgo = zeros(Int64,N); kgo[UDG] = delta[UDG];
    kco = kcv; myeps = 0.001;

    for i = 1:2N
        pcv = (1-kco) .* betav .* pcmax; qcv = (1-kco) .* betav .* qcmax;
        prv = (1-kgo) .* pgmax; qrv = (1-kgo) .* qgmax;
        pv = pcv - prv - pgv; qv = qcv - qrv - qgv;
        Pv = SuccMat * pv; Qv = SuccMat * qv;
        # vv = v0nom + deltaV0 * v0dis - R * pv - X * qv;

        om = Model(solver = GurobiSolver(OutputFlag = 0, genv));
        @variable(om, v[1:N]); vpar = zeros(AffExpr,N); vpar[2:N] = v[par[2:N]]; vpar[1] = v0nom + deltaV0 * v0dis;
        @constraint(om, voltDrop, (vpar - v - resistance.*Pv - reactance.*Qv)[noLTC] .== 0);
        if length(LTC) > 0
            @constraint(om, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
        end
        solve(om); vv = getvalue(v);

        # println(vv)
        genViolations = (1-kgo[UDG]).*(vgmin-vv - myeps)[UDG];
        genIndMax = indmax(genViolations);
        loadViolations = (1-kco[LOADS]).*(vcmin-vv - myeps)[LOADS];
        loadIndMax = indmax(loadViolations);
        # if genViolations[genIndMax] <= 0 && loadViolations[loadIndMax] <= 0
        #     break;
        # elseif genViolations[genIndMax] >= loadViolations[loadIndMax]
        #     kgo[UDG[genIndMax]] = 1;
        # else
        #     kco[LOADS[loadIndMax]] = 1;
        # end
        # if genViolations[genIndMax] > 0 && loadViolations[loadIndMax] <= 0
        if genViolations[genIndMax] > 0 && genViolations[genIndMax] >=0  loadViolations[loadIndMax]
            kgo[UDG[genIndMax]] = 1;
        elseif loadViolations[loadIndMax] > 0
            kco[LOADS[loadIndMax]] = 1;
        else
            break;
        end
        # println("Loads disconnected: ", nodes[kco.==1]);
        # println("DGs disconnected: ", nodes[kgo.==1]);
    end


    # println("vv : ", round(vv,2))
    # println("loads shed : ", nodes[kcnew.==1]);
    # println("DGs shed : ", nodes[kgnew.==1]);

    tv = maximum(abs.(vv - v0nom)); kcv, kgv = kco, kgo;
    P = SuccMat * pv; Q = SuccMat * qv;
    lacv = sum(WAC) * P[1];
    lovrv = WVR * sum(tv) / N;
    lcontv = sum((1 - kco) .* WLC .* (1 - betav));
    lsdv = sum(WSD .* kco);
    objv = lacv + lovrv + lcontv + lsdv;

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv
end

## To Be Done
function getOnlineCascadeSocp(delta, deltaV0)

    noAttackDelta = zeros(Int64,N); noDeltaV0 = 0
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(noAttackDelta, noDeltaV0);
    # println("vv: ",vv);
    # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    # lacv, lovrv, lcontv, lsdv = getPrecontingencyValues();
    # println("Loads disconnected: ", nodes[kcv.==1]);
    # println("DGs disconnected: ", nodes[kgv.==1]);
    # println("delta: ",delta)

    kgo = zeros(Int64,N); kgo[UDG] = delta[UDG];
    kco = kcv; myeps = 0.0001;

    for i = 1:2N
        pcv = (1-kco) .* betav .* pcmax; qcv = (1-kco) .* betav .* qcmax;
        prv = (1-kgo) .* pgmax; qrv = (1-kgo) .* qgmax;
        pv = pcv - prv - pgv; qv = qcv - qrv - qgv;
        Pv = SuccMat * pv; Qv = SuccMat * qv;
        # vv = v0nom + deltaV0 * v0dis - R * pv - X * qv;

        om = Model(solver = GurobiSolver(OutputFlag = 0, genv));
        @variable(om, v[1:N]); vpar = zeros(AffExpr,N); vpar[2:N] = v[par[2:N]]; vpar[1] = v0nom + deltaV0 * v0dis;
        @constraint(om, voltDrop, (vpar - v - resistance.*Pv - reactance.*Qv)[noLTC] .== 0);
        if length(LTC) > 0
            @constraint(om, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
        end
        solve(om); vv = getvalue(v);

        # println(vv)
        genViolations = (1-kgo[UDG]).*(vgmin-vv - myeps)[UDG];
        genIndMax = indmax(genViolations);
        loadViolations = (1-kco[LOADS]).*(vcmin-vv - myeps)[LOADS];
        loadIndMax = indmax(loadViolations);
        # if genViolations[genIndMax] <= 0 && loadViolations[loadIndMax] <= 0
        #     break;
        # elseif genViolations[genIndMax] >= loadViolations[loadIndMax]
        #     kgo[UDG[genIndMax]] = 1;
        # else
        #     kco[LOADS[loadIndMax]] = 1;
        # end
        # if genViolations[genIndMax] > 0 && loadViolations[loadIndMax] <= 0
        if genViolations[genIndMax] > 0 && genViolations[genIndMax] >=0  loadViolations[loadIndMax]
            kgo[UDG[genIndMax]] = 1;
        elseif loadViolations[loadIndMax] > 0
            kco[LOADS[loadIndMax]] = 1;
        else
            break;
        end
        # println("Loads disconnected: ", nodes[kco.==1]);
        # println("DGs disconnected: ", nodes[kgo.==1]);
    end


    # println("vv : ", round(vv,2))
    # println("loads shed : ", nodes[kcnew.==1]);
    # println("DGs shed : ", nodes[kgnew.==1]);

    tv = maximum(abs.(vv - v0nom)); kcv, kgv = kco, kgo;
    P = SuccMat * pv; Q = SuccMat * qv;
    lacv = sum(WAC) * P[1];
    lovrv = WVR * sum(tv) / N;
    lcontv = sum((1 - kco) .* WLC .* (1 - betav));
    lsdv = sum(WSD .* kco);
    objv = lacv + lovrv + lcontv + lsdv;

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv
end



function getOnlineCascadeAfterMicrogrid(delta, deltaV0)
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lovrv, lcontv, lsdv = evaluateSlaveModelWithMicrogridNoCascade(delta, deltaV0);

    kco = zeros(Int64,N); kgo = copy(delta);

    kcnew = vv .< vcmin; kgnew = vv .< vgmin;

    while (sum(kcnew) > sum(kco) || sum(kgnew) > sum(kgo))
        kco = kco | kcnew; kgo = kgo | kgnew;
        pcv = (1-kco).* betav .* pcmax; qcv = (1-kco).* betav .* qcmax;
        pgv = (1-kgo) .* pgmax; qgv = (1-kgo) .* qgmax;
        pv = pcv - prv;
        qv =  qcv - qrv;
        vv = v0nom - R * pv - X * qv;
        kcnew = kco | (vv .< vcmin); kgnew = kgo | (vv .< vgmin);
    end

    tv = max(vv-vmax,vmin - vv); kcv, kgv = kco, kgo;
    lovrv = WVR * sum(tv); lcontv = sum((1 - kco) .* WLC .* (1 - betav));
    lsdv = sum(WSD .* kco);
    objv = lovrv + lcontv + lsdv;

    pv, prv, qv, qrv, betav, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lovrv, lcontv, lsdv
end

function getDGUnaffected()
    delta = zeros(Int64,N); deltaV0 = 0; delta[DG] = 1;
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    # println("vv: ",vv);
    # pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    # lacv, lovrv, lcontv, lsdv = getPrecontingencyValues();
    # println("Loads disconnected: ", nodes[kcv.==1]);
    # println("DGs disconnected: ", nodes[kgv.==1]);
    # println("delta: ",delta)

    kco = kcv; kgo = copy(delta); myeps = 0.001;

    pcv = (1-kco) .* betamin .* pcmax; qcv = (1-kco) .* betav .* qcmax;
    pgv = (1-kgo) .* pgmax; qgv = (1-kgo) .* qgmax;
    pv = pcv - pgv; qv = qcv - qgv;
    Pv = SuccMat * pv; Qv = SuccMat * qv;
    # vv = v0nom + deltaV0 * v0dis - R * pv - X * qv;

    om = Model(solver = GurobiSolver(OutputFlag = 0, genv));
    @variable(om, v[1:N]); vpar = zeros(AffExpr,N); vpar[2:N] = v[par[2:N]]; vpar[1] = v0nom + deltaV0 * v0dis;
    @constraint(om, voltDrop, (vpar - v - resistance.*Pv - reactance.*Qv)[noLTC] .== 0);
    if length(LTC) > 0
        @constraint(om, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
    end
    solve(om); vv = getvalue(v);
    genViolations = (vgmin - vv)[DG];

    DG[genViolations .< 0]
end

function addDownstreamHeuristicCuts(mm, deltaVar)
    DGUF = getDGUnaffected();
    DGUF = DG;
    for k = 1:2 #length(DGUF)
        i = DGUF[k];
        for j = 1:N
            if SuccMat[i,j] == 1
                @constraint(mm, deltaVar[i] <= deltaVar[j]);
            end
        end
        # break;
    end
    mm
end


# Try a given number of random attacks of certain intensity, and evaluate
# the cascading impact
function getRandomAttackOnlineCascade(M, nAttacks)
    LLCbest = -infty;
    bestRandomDelta = zeros(Int64,N,1);
    bestkco = 0; bestkgo = 0; bestvv = 0; bestprv = 0; bestqrv = 0; bestbetav = 0;
    for i = 1:nAttacks
        delta = zeros(Int64,N);
        delta[randperm(N)[1:M]] = 1;
        kco, kgo, vv, prv, qrv, betav = getOnlineCascade(delta);
        LLCcur = sum(WSD .* kco + (1 - kco) .* WLC .* (1 - betav));
        if LLCcur > LLCbest
            LLCbest = LLCcur;
            bestkco, bestkgo, bestvv, bestprv, bestqrv, bestbetav = kco, kgo, vv, prv, qrv, betav;
            bestRandomDelta = delta;
        end
    end

    bestkco, bestkgo, bestvv, bestprv, bestqrv, bestbetav, bestRandomDelta
end

# Brute force algorithm -- needs to be reworked
function getBruteForceCascadeAttack(M)

    bestDelta = zeros(Int64,N); bestDeltaV0 = 0;
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithCascade(bestDelta, bestDeltaV0);

    bestobjv = -infty; bestDeltaV0 = 0; bestDelta = zeros(Int64,N);

    ludg = length(UDG);
    dgAttack1 = dgAttack[1:ludg]; dgAttack2 = dgAttack[ludg+1:2ludg]; dgAttack3 = dgAttack[2ludg+1:3ludg];
    bestObjv = -infty; nodes = collect(1:N);
    for curdeltaV0 = 0:0 #-1:1
        # JuMP.setRHS(subvolt, v0nom + v0dis * curdeltaV0);

        for attackSet in collect(combinations(UDG,M))
            delta = zeros(Int64,N); delta[attackSet] = 1;

            # for j = 1:length(UDG)
            #     i = UDG[j];
            #
            #     JuMP.setRHS(dgAttack1[j], -delta[i]);
            #     JuMP.setRHS(dgAttack2[j], pgmax[i] * (1-delta[i]));
            #     JuMP.setRHS(dgAttack3[j], qgmax[i] * (1-delta[i]));
            # end

            # sm, p, pr, q, qr, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithCascade(delta, curdeltaV0);
            # status = solve(sm);
            pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, curdeltaV0);

            if objv >= bestObjv
                bestObjv = objv; bestDelta = copy(delta); bestDeltaV0 = curdeltaV0;
            end
        end
    end

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(bestDelta, bestDeltaV0);
    println("Brute force results");
    println("Nodes attacked : ", nodes[bestDelta.>0])
    println("LOADS disconnected : ", nodes[kcv .> 0.5]);
    println("DG disconnected : ", nodes[kgv .> 0.5]);
    println("LOADS shed : ", nodes[betav .< 1]);
    println("LLC : ", sum((WSD-WLC) .* kcv) + sum(WLC .* (1-betav)));
    println("LLCmax : ", LLCmax);
    println("****************************************");

    bestDelta, bestDeltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv
end

function getBruteForceCascadeSocpAttack(M)

    bestDelta = zeros(Int64,N); bestDeltaV0 = 0;
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack, ell = getSlaveSocpModelWithCascade(bestDelta, bestDeltaV0);

    bestobjv = -infty; bestDeltaV0 = 0; bestDelta = zeros(Int64,N);

    ludg = length(UDG);
    dgAttack1 = dgAttack[1:ludg]; dgAttack2 = dgAttack[ludg+1:2ludg]; dgAttack3 = dgAttack[2ludg+1:3ludg];
    bestObjv = -infty; nodes = collect(1:N);
    for curdeltaV0 = 0:0 #-1:1
        # JuMP.setRHS(subvolt, v0nom + v0dis * curdeltaV0);

        for attackSet in collect(combinations(UDG,M))
            delta = zeros(Int64,N); delta[attackSet] = 1;

            # for j = 1:length(UDG)
            #     i = UDG[j];
            #
            #     JuMP.setRHS(dgAttack1[j], -delta[i]);
            #     JuMP.setRHS(dgAttack2[j], pgmax[i] * (1-delta[i]));
            #     JuMP.setRHS(dgAttack3[j], qgmax[i] * (1-delta[i]));
            # end

            # sm, p, pr, q, qr, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithCascade(delta, curdeltaV0);
            # status = solve(sm);
            pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveSocpModelWithCascade(delta, curdeltaV0);

            if objv >= bestObjv
                bestObjv = objv; bestDelta = copy(delta); bestDeltaV0 = curdeltaV0;
            end
        end
    end

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveSocpModelWithCascade(bestDelta, bestDeltaV0);
    println("Brute force results");
    println("Nodes attacked : ", nodes[bestDelta.>0])
    println("LOADS disconnected : ", nodes[kcv .> 0.5]);
    println("DG disconnected : ", nodes[kgv .> 0.5]);
    println("LOADS controlled : ", nodes[betav .< 1]);
    println("LLC : ", sum((WSD-WLC) .* kcv) + sum(WLC .* (1-betav)));
    println("LLCmax : ", LLCmax);
    println("****************************************");

    bestDelta, bestDeltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)
    nr = 2;
    println("voltage : ", round.(vv,nr));
    println("t : ", round.(tv,nr));
    println("objv : ", round(objv,nr));

    println("kcv : ", kcv);
    println("kgv : ", kgv);
    # println("dobjv : ", dobjv);

    println("betav : ", round.(betav,nr));
    println("Pv : ", round.(Pv,nr));
    println("Qv : ", round.(Qv,nr));

    println("pcv : ", round.(pcv,nr));
    println("pgv : ", round.(pgv,nr));
    println("prv : ", round.(prv,nr));
    println("qcv : ", round.(qcv,nr));
    println("qgv : ", round.(qgv,nr));
    println("qrv : ", round.(qrv,nr));

    println("pv : ", round.(pv,nr));
    println("qv : ", round.(qv,nr));

end

function printSocpResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, ellv)
    nr = 2;
    printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv);
    println("ellv : ", round.(ellv,nr));
end

global bestObjvUntilNow = 0;

function getBendersCascadeIteration(mm, deltaVar, deltaV0Var, verbose, LLCreq)
    hasLLCExceeded = false; isBendersCutAdded = false; isMasterModelNotSolvable = false;

    pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    # println(mm);
    if verbose
        println(mm);
    end

    mstatus = solve(mm);
    if mstatus != :Optimal
        println("Master model solve status : ", mstatus);
        isMasterModelNotSolvable = true;
        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv,
        qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv,
        kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0;
    end

    delta = round.(Int64,getvalue(deltaVar)); deltaV0 = getvalue(deltaV0Var);
    nNodesAttacked = sum(round.(Int64,delta));
    # println("delta : ", delta);
    if verbose
        println("no. of nodes attacked : ", nNodesAttacked);
        println("delta : ", nodes[delta .== 1]);
    end

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithCascade(delta, deltaV0);

    if verbose
        println("Updated slave model");
        println(sm);
    end

    sstatus = solve(sm);
    if sstatus != :Optimal
        println("Slave model solve status : ", sstatus);
        return;
    end

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objvMIP, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v,  t, kc, kg, vpar);

    verbose ? printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv) : nothing;

    # curLLC = sum((WSD-WLC).*kcv) + sum(WLC.*(1-betav))
    curLLC = objvMIP;
    if curLLC >= LLCreq
        verbose? println("LLCreq threshold exceeded : ", curLLC) : nothing;
        hasLLCExceeded = true;
    end

    # println("vv : ", vv);
    if verbose
        println("LOADS disconnected : ", nodes[kcv .> 0.5]);
        println("DG disconnected : ", nodes[kgv .> 0.5]);
        println("LOADS shed : ", nodes[betav .< 1], "\ntotal load shed : ", sum(WLC.*(1-betav)));
        println("LLC : ", curLLC);
        println("LLCreq : ", LLCreq);
        println("objv : ", objv);
    end

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithoutCascade(kcv, kgv, delta, deltaV0);

    if verbose || true
        println("****************************************");
        println("Printing slave model without cascade");
        println(sm);
    end

    sstatus = solve(sm);
    objv = getobjectivevalue(sm);

    ineqlbar = getdual(ineq);
    eqnubar = getdual(eq);
    deltaIneqlbar = getdual(deltaineq);
    deltaeqnubar = getdual(deltaeq);

    deltaIneqVar = [deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]];
    global eta
    deltaIneqMul = [ov[UDG]; -pgmax[UDG]; -qgmax[UDG]; ov[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];
    # println("sizes of deltaIneqvar", size(deltaIneqVar),", ", size(deltaIneqMul),", ", size(deltaIneqlbar));
    # deltaIneqMul = [ov[UDG]; -(1 - kgv[UDG]) .* pgmax[UDG]; -(1 - kgv[UDG]) .* qgmax[UDG]];
    # deltaIneqVar = [deltaVar[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];
    # deltaIneqMul = [uv[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];
    # println("deltaIneqlbar: ",deltaIneqlbar);

    deltaeqVar = [deltaV0Var];
    deltaeqVarMul = [v0dis];

    if verbose
        println("lengths :", length(deltaeqnubar), " ", length(deltaeqVarMul), " ", length(deltaeqVar));
        println("lengths :", length(deltaIneqlbar), " ", length(deltaIneqMul), " ", length(deltaIneqVar));
        println(length(eqnubar), " ", length(eqb), " ", length(ineqlbar), " ",
        length(ineqb));

        println("deltaeqnubar : ", round.(deltaeqnubar,3));
        println("deltaIneqlbar : ", round.(deltaIneqlbar,3));
        println("vv : ", vv);
        println("lacv, lovrv, lcontv, lsdv: ", lovrv, ", ", lcontv, ", ", lsdv)
        println("sumprod of nu and rhs : ", sum(eqnubar.*eqb));
        println("sumprod of lambda and rhs : ", sum(ineqlbar.*ineqb));
        println("deltaeqnubar: ",deltaeqnubar);
        println("deltaIneqlbar: ",deltaIneqlbar);
        println("sum(WLC) : ", sum(WLC));
        println("sum((WSD-WLC).*kcv) : ", sum((WSD-WLC).*kcv));
        println("sum of lhs : ", sum((WSD-WLC).*kcv) + sum(WLC[LOADS]) +sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb))
        println("objv : ", objv);
    end

    if maximum(abs.(deltaeqnubar)) + maximum(abs.(deltaIneqlbar)) >= myinf
        # bendersCutLHS = sum((WSD-WLC).*kcv) + sum(WLC[LOADS]) + sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb) + sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        bendersCutLHS = sum(WLC[LOADS]) + sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb) + sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        if !trial
            bendersCutLHS -= (objv + epsilon);
            # bendersCutLHS -= LLCreq;
            println("bendersCut\n",bendersCutLHS)
        else
            dconst = (sum((WSD-WLC).*kcv) + sum(WLC[LOADS]));
            diff = objv - dconst;
            eta = 0.98; rhs = epsilon;
            # rhs += dconst + max(diff, eta * diff + (1-eta)*300);
            rhs += objv - 5 * sum(deltaeqnubar .* deltaeqVarMul);
            # rhs += objv;
            bendersCutLHS -= rhs;
        end
        @constraint(mm, bendersCutLHS >= 0);
        # println("objv, lacv, lovrv, lcontv, lsdv : ", round.([objv, lacv, lovrv, lcontv, lsdv], 2));
        # println("nodes attacked: ", nodes[delta.== 1]);
        # println("bendersCut-objv: ", (bendersCutLHS));
        # @constraint(mm, bendersCutLHS >= objv + epsilon);
        # global bestObjvUntilNow
        global bestObjvUntilNow = max(bestObjvUntilNow, objv) + LLCmax / 1000;
        # @constraint(mm, bendersCutLHS >= LLCreq*(1 + epsilon));
        # println(100(1-bestObjvUntilNow/LLCmax));
        # epsilon = 0.1^5;
        # if N < 34
        # @constraint(mm, bendersCutLHS >= objv + epsilon);
        # elseif N < 69
        #     # @constraint(mm, bendersCutLHS >= objv + LLCmax / 500);
        #     @constraint(mm, bendersCutLHS >= objv + epsilon);
        # elseif N == 69
        #     # @constraint(mm, bendersCutLHS >= objv + LLCmax / 8000);
        #     @constraint(mm, bendersCutLHS >= objv + epsilon);
        # elseif N >= 118
        #     @constraint(mm, bendersCutLHS >= objv + LLCmax / 20000);
        # end
        # if sum(delta) <= length(UDG)-2
        #     @constraint(mm, bendersCutLHS >= objv + LLCmax / 500);
        #     # @constraint(mm, bendersCutLHS >= objv + epsilon);
        # else
        #     @constraint(mm, bendersCutLHS >= objv + epsilon);
        # end

        # @constraint(mm, bendersCutLHS >= bestObjvUntilNow );
        # @constraint(mm, bendersCutLHS >= LLCreq);

        isBendersCutAdded = true;
    else
        println("Not a good cut ", maximum(abs(eqnubar)));
        println("eqnubar: ", eqnubar);
        printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
        objv, pgv, qgv);
    end

    hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, delta, deltaV0, ellv
end

# global isKSpecified = true, M = 15;
function getBendersMethodCascade(LLCreq)
    kcval = ov; kgval = ov; delta = ov; deltaV0 = 0; cardinality = 0;
    withCascade = true; verbose = false; curLLC = 0;
    nNodesAttacked = zeros(Int64, N); printSummary = true;
    nIters = [0]; nSecs = [0]; nCurrentIter = 0; nCardinalities = [0]; objv = 0;
    resiliences = [100];
    global bestObjvUntilNow = 0;

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    # @constraint(mm, sum(deltaVar) >= kpercent/100 * length(UDG));
    N >= 36 ? mm = addDownstreamHeuristicCuts(mm, deltaVar) : nothing;
    # @constraint(mm, deltaV0)
    # println(mm)

    println("Entering loop");
    jmax = N >= 36? 5 : N >= 24? 800 : 100;
    lastChange = 0;
    j = 0; cumTime = 0; tic();
    for j = 1:jmax
        # println(mm);
        nCurrentIter += 1;

        if sum(delta) > cardinality || true

            currTime = toq(); cumTime += currTime;
            nSecs = [nSecs; cumTime];
            nIters = [nIters; j];
            cardinality = sum(delta);
            nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience]
            # println(resilience, ", (", j, "), ", round(cumTime,2), ", ", round(100cardinality/N,2));
            lastChange = nCurrentIter;


            tic();
        end

        if j - lastChange >= 150
            # @constraint(mm, sum(deltaVar) >= cardinality + 1);
            # println(mm);
        end

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, deltav, deltaV0v = getBendersCascadeIteration(mm, deltaVar, deltaV0Var, verbose, LLCreq);

        # println("Iteration : ", j);
        j % 50 == 0 ? print(j," ") : nothing;
        # isMasterModelNotSolvable? println(""): println("Nodes attacked: ", nodes[deltav .> 0.9]);
        # println("hasLLCExceeded, isMasterModelNotSolvable, isBendersCutAdded : ", hasLLCExceeded, " ", isMasterModelNotSolvable, " ", isBendersCutAdded);
        if isBendersCutAdded || hasLLCExceeded
            deltaV0 = deltaV0v; delta = deltav;
        end

        if isMasterModelNotSolvable || !isBendersCutAdded || hasLLCExceeded || j == jmax
            println("\nExit reasons : isMasterModelNotSolvable, isBendersCutAdded, hasLLCExceeded, IterLimitReached ", isMasterModelNotSolvable, " ", isBendersCutAdded, " ", hasLLCExceeded, " ", j == jmax);
            currTime = toq(); cumTime += currTime; nSecs = [nSecs; cumTime]; nIters = [nIters; j+1]; cardinality = sum(delta); nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience]
            println(resilience, ", (", j+1, "), ", round(cumTime,2), ", ", round(100cardinality/N,2));
            lastChange = nCurrentIter;
            break;
        end

    end
    println("Exited loop");

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);

    if printSummary
        println("Benders cascade method results");
        println("Nodes attacked : ", nodes[delta.>0])
        println("Loads disconnected : ", nodes[kcv .> 0.5]);
        println("DGs disconnected : ", nodes[kgv .> 0.5]);
        println("loads shed : ", nodes[betav .< 1]);
        println("LLC : ", objv);
        println("LLCreq : ", LLCreq);
        println("resilience: ",floor(100(1-objv/LLCmax),2));
        println("****************************************");
    end

    if printPerformance

        nSecs = round.(nSecs,2);
        println(nCardinalities)
        # nCardinalities = round.(100nCardinalities/N,1);
        nCardinalities = round.(Int64,nCardinalities);
        println(nCardinalities);
        ploss = round.(100 - resiliences,2);
        for ij = 1:length(nIters)
            println(ploss[ij], ", (",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], " \\\\");
        end
    end

    global nBendersCascadeIter = j;
    llc = objv;
    delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv
end

function getBendersCascadeIterationForSocp(mm, deltaVar, deltaV0Var, verbose, LLCreq)
    println("In getBendersCascadeIterationForSocp");
    hasLLCExceeded = false; isBendersCutAdded = false; isMasterModelNotSolvable = false;

    pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0, ellv = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    # println(mm);
    if verbose
        println(mm);
    end

    mstatus = solve(mm);

    if mstatus != :Optimal
        println("Master model solve status : ", mstatus);
        isMasterModelNotSolvable = true;
        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv,
        qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0;
    end

    delta = round.(Int64,getvalue(deltaVar)); deltaV0 = getvalue(deltaV0Var);
    nNodesAttacked = sum(round.(Int64,delta));
    # println("delta : ", delta);
    if verbose
        println("no. of nodes attacked : ", nNodesAttacked);
        println("delta : ", nodes[delta .== 1]);
    end

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack, ell = getSlaveSocpModelWithCascade(delta, deltaV0);

    if verbose
        println("Updated slave model");
        # println(sm);
    end

    sstatus = solve(sm);
    if sstatus != :Optimal
        println("Slave model solve status : ", sstatus);
        return delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv, ellv

    end

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objvMIP, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValuesSocp(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v, t, kc, kg, ell, vpar);

    verbose ? printSocpResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objvMIP, pgv, qgv, ellv) : nothing;

    # curLLC = sum((WSD-WLC).*kcv) + sum(WLC.*(1-betav))
    curLLC = objvMIP;
    if curLLC >= LLCreq
        verbose? println("LLCreq threshold exceeded : ", curLLC) : nothing;
        hasLLCExceeded = true;
    end

    # println("vv : ", vv);
    if verbose
        println("LOADS disconnected : ", nodes[kcv .> 0.5]);
        println("DG disconnected : ", nodes[kgv .> 0.5]);
        println("LOADS shed : ", nodes[betav .< 1], "\ntotal load shed : ", sum((WLC.*(1-betav))[LOADS]));
        println("LLC : ", curLLC);
        println("LLCreq : ", LLCreq);
        println("objv : ", objvMIP);
    end

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ell, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, socpIneq, socpIneqb, subvolt, deltaineq, deltaineqb, dgAttack = getSlaveSocpModelWithoutCascade(kcv, kgv, delta, deltaV0);

    if verbose
        println("****************************************");
        println("Printing slave model without cascade");
        println(sm);
    end

    sstatus = solve(sm);
    assert(sstatus == :Optimal)
    objv = getobjectivevalue(sm);

    verbose ? println("objvMIP = ", objvMIP, ", objv = ",objv) : nothing;
    verbose ? printSocpResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, ellv) : nothing;



    ineqlbar = getdual(ineq);
    eqnubar = getdual(eq);
    deltaIneqlbar = getdual(deltaineq);
    deltaeqnubar = getdual(deltaeq);
    socpIneqtbar = getdual(socpIneq);

    deltaIneqVar = [deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]];
    global eta
    deltaIneqMul = [ov[UDG]; -pgmax[UDG]; -qgmax[UDG]; ov[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];

    #preprocessing of deltaIneqlbar
    ludg = length(UDG);
    row = 2
    deltaIneqlbar[(row-1)*ludg+1:(row)*ludg] +=  deltaIneqlbar[(2row-1)*ludg+1:(2row)*ludg];
    row = 3
    deltaIneqlbar[(row-1)*ludg+1:(row)*ludg] +=  deltaIneqlbar[(2row-1)*ludg+1:(2row)*ludg];

    # deltaIneqMul = [ov[UDG]; -(1 - kgv[UDG]) .* pgmax[UDG]; -(1 - kgv[UDG]) .* qgmax[UDG]];
    # deltaIneqVar = [deltaVar[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];
    # deltaIneqMul = [uv[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];
    # println("deltaIneqlbar: ",deltaIneqlbar);

    deltaeqVar = [deltaV0Var];
    deltaeqVarMul = [v0dis];

    if verbose
        println("lengths :", length(deltaeqnubar), " ", length(deltaeqVarMul), " ", length(deltaeqVar));
        println("lengths :", length(deltaIneqlbar), " ", length(deltaIneqMul), " ", length(deltaIneqVar));
        println(length(eqnubar), " ", length(eqb), " ", length(ineqlbar), " ",
        length(ineqb));

        println("deltaeqnubar : ", round.(deltaeqnubar,3));
        println("deltaIneqlbar : ", round.(deltaIneqlbar,3));
        println("deltaIneqVar : ", deltaIneqVar);
        println("vv : ", vv);
        println("lacv, lovrv, lcontv, lsdv: ", lovrv, ", ", lcontv, ", ", lsdv)
        println("sumprod of nu and rhs : ", sum(eqnubar.*eqb));
        println("sumprod of lambda and rhs : ", sum(ineqlbar.*ineqb));
        println("deltaeqnubar: ",deltaeqnubar);
        println("deltaIneqlbar: ",deltaIneqlbar);
        println("sum(WLC)[LOADS] : ", sum(WLC[LOADS]));
        println("sum((WSD-WLC).*kcv) : ", sum((WSD-WLC).*kcv));
        println("sum of lhs : ", sum(((WSD-WLC).*kcv)[LOADS]) + sum(WLC[LOADS]) +sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb))
        println("objv : ", objv);
    end

    if maximum(abs.(deltaeqnubar)) + maximum(abs.(deltaIneqlbar)) >= myinf
        bendersCutLHS = sum((WSD-WLC).*kcv) + sum(WLC[LOADS]) + sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb) + sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);# + sum(socpIneqtbar.*socpIneqb);
        if !trial
            #epsilon = 20; # for N = 12 or 24
            epsilon = 20
            bendersCutLHS -= (objv + epsilon);
            # bendersCutLHS -= LLCreq;
        else
            dconst = (sum((WSD-WLC).*kcv) + sum(WLC[LOADS]));
            diff = objv - dconst;
            eta = 0.98; rhs = epsilon;
            # rhs += dconst + max(diff, eta * diff + (1-eta)*300);
            rhs += objv - 5 * sum(deltaeqnubar .* deltaeqVarMul);
            # rhs += objv;
            bendersCutLHS -= rhs;
        end
        # println("size of deltavar ",size(deltaVar),", ",size(delta))
        @constraint(mm, sum(deltaVar[delta.==0]) + sum(1-deltaVar[delta.==1]) >= 1)
        @constraint(mm, bendersCutLHS >= 0);
        # println("objv, lacv, lovrv, lcontv, lsdv : ", round.([objv, lacv, lovrv, lcontv, lsdv], 2));
        # println("nodes attacked: ", nodes[delta.== 1]);
        println("bendersCut-objv: ", (bendersCutLHS));
        # @constraint(mm, bendersCutLHS >= objv + epsilon);
        # global bestObjvUntilNow
        global bestObjvUntilNow = max(bestObjvUntilNow, objv) + LLCmax / 1000;
        # @constraint(mm, bendersCutLHS >= LLCreq*(1 + epsilon));
        # println(100(1-bestObjvUntilNow/LLCmax));
        # epsilon = 0.1^5;
        # if N < 34
        # @constraint(mm, bendersCutLHS >= objv + epsilon);
        # elseif N < 69
        #     # @constraint(mm, bendersCutLHS >= objv + LLCmax / 500);
        #     @constraint(mm, bendersCutLHS >= objv + epsilon);
        # elseif N == 69
        #     # @constraint(mm, bendersCutLHS >= objv + LLCmax / 8000);
        #     @constraint(mm, bendersCutLHS >= objv + epsilon);
        # elseif N >= 118
        #     @constraint(mm, bendersCutLHS >= objv + LLCmax / 20000);
        # end
        # if sum(delta) <= length(UDG)-2
        #     @constraint(mm, bendersCutLHS >= objv + LLCmax / 500);
        #     # @constraint(mm, bendersCutLHS >= objv + epsilon);
        # else
        #     @constraint(mm, bendersCutLHS >= objv + epsilon);
        # end

        # @constraint(mm, bendersCutLHS >= bestObjvUntilNow );
        # @constraint(mm, bendersCutLHS >= LLCreq);

        isBendersCutAdded = true;
    else
        println("Not a good cut ", maximum(abs(eqnubar)));
        println("eqnubar: ", eqnubar);
        printSocpResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, ellv);
    end

    hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, delta, deltaV0, ellv
end

function getBendersMethodCascadeForSocp(LLCreq)
    kcval = ov; kgval = ov; delta = ov; deltaV0 = 0; cardinality = 0;
    withCascade = true; verbose = false; curLLC = 0;
    nNodesAttacked = zeros(Int64, N); printSummary = true;
    nIters = [0]; nSecs = [0]; nCurrentIter = 0; nCardinalities = [0]; objv = 0;
    resiliences = [100];
    global bestObjvUntilNow = 0;

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    # @constraint(mm, sum(deltaVar) >= kpercent/100 * length(UDG));
    # N >= 36 ? mm = addDownstreamHeuristicCuts(mm, deltaVar) : nothing;
    # @constraint(mm, deltaV0)
    # println(mm)

    println("Entering loop");
    jmax = N >= 36? 50 : N >= 24? 800 : N == 12? 64 : 100;
    lastChange = 0;
    j = 0; cumTime = 0; tic();
    for j = 1:jmax
        # println(mm);
        nCurrentIter += 1;

        # bookkeeping
        if sum(delta) > cardinality || true

            currTime = toq(); cumTime += currTime;
            nSecs = [nSecs; cumTime];
            nIters = [nIters; j];
            cardinality = sum(delta);
            nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience]
            # println(resilience, ", (", j, "), ", round(cumTime,2), ", ", round(100cardinality/N,2));
            lastChange = nCurrentIter;

            tic();
        end

        if j - lastChange >= 150
            # @constraint(mm, sum(deltaVar) >= cardinality + 1);
            # println(mm);
        end

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, deltav, deltaV0v, ellv = getBendersCascadeIterationForSocp(mm, deltaVar, deltaV0Var, verbose, LLCreq);

        # println("Iteration : ", j);
        j % 50 == 0 ? print(j," ") : nothing;
        # isMasterModelNotSolvable? println(""): println("Nodes attacked: ", nodes[deltav .> 0.9]);
        # println("hasLLCExceeded, isMasterModelNotSolvable, isBendersCutAdded : ", hasLLCExceeded, " ", isMasterModelNotSolvable, " ", isBendersCutAdded);
        if isBendersCutAdded || hasLLCExceeded
            deltaV0 = deltaV0v; delta = deltav;
        end

        if isMasterModelNotSolvable || !isBendersCutAdded || hasLLCExceeded || j == jmax
            println("\nExit reasons : isMasterModelNotSolvable, isBendersCutAdded, hasLLCExceeded, IterLimitReached ", isMasterModelNotSolvable, " ", isBendersCutAdded, " ", hasLLCExceeded, " ", j == jmax);

            currTime = toq(); cumTime += currTime; nSecs = [nSecs; cumTime]; nIters = [nIters; j+1]; cardinality = sum(delta); nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience]

            println(resilience, ", (", j+1, "), ", round(cumTime,2), ", ", round(100cardinality/N,2));
            lastChange = nCurrentIter;
            break;
        end

    end
    println("Exited loop");

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveSocpModelWithCascade(delta, deltaV0);

    if printSummary
        println("Benders cascade method for Socp model results");
        println("Nodes attacked : ", nodes[delta.>0])
        println("Loads disconnected : ", nodes[kcv .> 0.5]);
        println("DGs disconnected : ", nodes[kgv .> 0.5]);
        println("loads shed : ", nodes[betav .< 1]);
        println("LLC : ", objv);
        println("LLCreq : ", LLCreq);
        println("resilience: ",floor(100(1-objv/LLCmax),2));
        println("****************************************");
    end

    if printPerformance

        nSecs = round.(nSecs,2);
        # println(nCardinalities)
        # nCardinalities = round.(100nCardinalities/N,1);
        nCardinalities = round.(Int64,nCardinalities);
        println("nCardinalities:\n", nCardinalities);
        ploss = round.(100 - resiliences,2);
        for ij = 1:length(nIters)
            println(ploss[ij], ", (",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], " \\\\");
        end
    end

    global nBendersCascadeIter = j;
    llc = objv;

    delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getLinearAssistBendersCascade(LLCreq)
    deltaV0 = 0; delta = zeros(Int64,N); delta[UDG] = 1; kL = 0; kR = 100; ludg = length(UDG); minkpercent = 100;
    # sm, p, pr, q, qr, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq,
    # eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, dgAttack =
    # while ludg * (kR - kL) >= 1
    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    if objv < LLCreq
        println("LLCreq not achievable");
        minkpercent = 100;
    else
        for i = 1:length(UDG)
            delta[UDG[i]] = 0;
            pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
            lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
            minkpercent = 100sum(delta)/length(UDG) - 0.0001;
            println("objv, llcreq: ",objv," ",LLCreq);
            if objv < LLCreq
                break;
            end
        end
    end
    println("minkpercent: ", minkpercent);
    delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq, minkpercent);

    delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv
end

function getBendersCascadeMinNodesLLC()
    kcval = ov; kgval = ov; delta = ov;
    withCascade = true; verbose = false; curLLC = 0;
    minNodesv = zeros(Int64,1); sumkcvs = zeros(Int64,1); sumkgvs = zeros(Int64,1);
    minDeltav = zeros(Int64,N);
    objvs = []; tempLLCreq = 0;
    printSummary = true; nIters = [0]; nSecs = [0]; nCurrentIter = 0; nCardinalities = [0]; lastChange = 0; cardinality = 0; objv = 0;

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    # mm = addDownstreamHeuristicCuts(mm, deltaVar);

    println("Entering loop");
    jmax = N > 12? 5000 : 1000; cumTime = 0; tic();
    for j = 1:jmax
        println(j);
        nCurrentIter += 1;
        if sum(delta) > cardinality

            currTime = toq(); cumTime += currTime;
            nSecs = [nSecs; cumTime];
            nIters = [nIters; j];
            cardinality = sum(delta);
            nCardinalities = [nCardinalities; cardinality];
            println("resilience: ",floor(100(1-objv/LLCmax),2));
            lastChange = nCurrentIter;
            tic();
        end

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, deltav, deltaV0v = getBendersCascadeIteration(mm, deltaVar,
        deltaV0Var, verbose, tempLLCreq);

        tempLLCreq = objv;

        if !isMasterModelNotSolvable
            delta = deltav; deltaV0 = deltaV0v;
        end
        if !isBendersCutAdded || isMasterModelNotSolvable
            println("Exit reasons : isMasterModelNotSolvable, isBendersCutAdded, hasLLCExceeded ", isMasterModelNotSolvable, " ", isBendersCutAdded, " ", hasLLCExceeded);
            currTime = toq(); cumTime += currTime; nSecs = [nSecs; cumTime]; nIters = [nIters; j]; cardinality = sum(delta); nCardinalities = [nCardinalities; cardinality];
            println("resilience: ",floor(100(1-objv/LLCmax),2));
            lastChange = nCurrentIter;
            break;
        end

        objvs = vcat(objvs, objv);

        if j > 1
            minDeltav = hcat(minDeltav, delta); sumkcvs = vcat(sumkcvs, sum(kcv));
            sumkgvs = vcat(sumkgvs, sum(kgv));
            nNodesAttacked = sum(delta);
            minNodesv = vcat(minNodesv,nNodesAttacked);
        end


    end

    if printPerformance

        nSecs = round.(nSecs,2);
        println(nCardinalities)
        nCardinalities = round.(100nCardinalities/N,1);
        println(nCardinalities)
        for ij = 1:length(nIters)
            println("(",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], " \\\\");
        end
    end

    minNodesv, minDeltav, objvs, sumkcvs, sumkgvs
end

function getBendersSocpCascadeMinNodesLLC()
    kcval = ov; kgval = ov; delta = ov;
    withCascade = true; verbose = false; curLLC = 0;
    minNodesv = zeros(Int64,1); sumkcvs = zeros(Int64,1); sumkgvs = zeros(Int64,1);
    minDeltav = zeros(Int64,N);
    objvs = []; tempLLCreq = 0;
    printSummary = true; nIters = [0]; nSecs = [0]; nCurrentIter = 0; nCardinalities = [0]; lastChange = 0; cardinality = 0; objv = 0;

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    # mm = addDownstreamHeuristicCuts(mm, deltaVar);

    println("Entering loop");
    jmax = N > 12? 5000 : 1000; cumTime = 0; tic();
    for j = 1:jmax
        println(j);
        nCurrentIter += 1;
        if sum(delta) > cardinality

            currTime = toq(); cumTime += currTime;
            nSecs = [nSecs; cumTime];
            nIters = [nIters; j];
            cardinality = sum(delta);
            nCardinalities = [nCardinalities; cardinality];
            println("resilience: ",floor(100(1-objv/LLCmax),2));
            lastChange = nCurrentIter;
            tic();
        end

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, deltav, deltaV0v, ellv = getBendersCascadeIterationForSocp(mm, deltaVar, deltaV0Var, verbose, tempLLCreq);

        tempLLCreq = objv;

        if !isMasterModelNotSolvable
            delta = deltav; deltaV0 = deltaV0v;
        end
        if !isBendersCutAdded || isMasterModelNotSolvable
            println("Exit reasons : isMasterModelNotSolvable, isBendersCutAdded, hasLLCExceeded ", isMasterModelNotSolvable, " ", isBendersCutAdded, " ", hasLLCExceeded);
            currTime = toq(); cumTime += currTime; nSecs = [nSecs; cumTime]; nIters = [nIters; j]; cardinality = sum(delta); nCardinalities = [nCardinalities; cardinality];
            println("resilience: ",floor(100(1-objv/LLCmax),2));
            lastChange = nCurrentIter;
            break;
        end

        objvs = vcat(objvs, objv);

        if j > 1
            minDeltav = hcat(minDeltav, delta); sumkcvs = vcat(sumkcvs, sum(kcv));
            sumkgvs = vcat(sumkgvs, sum(kgv));
            nNodesAttacked = sum(delta);
            minNodesv = vcat(minNodesv,nNodesAttacked);
        end


    end

    if printPerformance

        nSecs = round.(nSecs,2);
        println(nCardinalities)
        nCardinalities = round.(100nCardinalities/N,1);
        println(nCardinalities)
        for ij = 1:length(nIters)
            println("(",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], " \\\\");
        end
    end

    minNodesv, minDeltav, objvs, sumkcvs, sumkgvs
end

function setGlobalParameters()
    # println("Reached 1");
    global vmin, vmax, vmmin, vmmax, vgmin, vcmin, WAC, WLC, WVR, WSD, ov, uv
    global xrratio, rc, Cloadc, resistance, reactance, srmax
    global infty, tol, epsilon, I, O, fmin, fmax, WLC, WSD, nodes, WMG, WEG
    global path, commonPath, R, X, betamin, vgmax, vcmax, v0nom, v0dis
    global pcmax, qcmax, LLCmax, pgmax, qgmax, betamax, delta, pdmax, qdmax, Smax, pgnom, qgnom, vref, mp, mq

    global myinf, leaves, SDI, noSDI, MGL, noMGL, par, SuccMat, RES, EG, semax, DG, UDG, LOADS, CHILD
    global LTC, noLTC, N

    ov = zeros(Float64,N); uv = ones(Float64,N); v0dis = 0.03; epsilon = 0.1^5; infty = 10^5+1; tol = 0.1^3; myinf = 0.1^7; v0nom = 1; I = eye(N); O = zeros(Float64,N,N); vmin = 0.95uv; vmax = 1.05uv; WVR = 100; Smax = 2uv; vref = 1.00; vmmmax = 2uv;
    WMG = zeros(Float64,N); WEG = zeros(Float64,N); WLC = zeros(Float64,N); WSD = zeros(Float64,N); pcmax = zeros(Float64,N); qcmax = zeros(Float64,N); pgmax = zeros(Float64,N); qgmax = zeros(Float64,N);

    nodes = collect(1:N); EG = copy(nodes); LOADS = copy(nodes); LTC = zeros(Int64,0); DG = copy(nodes); SDI = nodes[nodes.%2.==1]; noSDI = setdiff(nodes,SDI); MGL = zeros(Int64,0);
    par = collect(0:N-1);
    if N%3 == 0
        MGL = [1; Int(N/3)+1; 2Int(N/3)+1];
    end
    if N%2 == 0
        EG = nodes[nodes.%2 .== 1];
    end
    # EG = union(EG, MGL);


    srand(716);
    # println(nodes[randperm(N)], " ", round(Int64,N/2));
    # DG = (nodes[randperm(N)])[1:round(Int64,N/2)];
    DG = sort((nodes[randperm(N)])[1:round(Int64,N/2)]);

    println("DG ", DG);

    if N == 3
        par[3] = 1;
        xrratio = 2; rc = 0.1; srmax = 0.2; Cloadc = 100;
        DG = copy(nodes);

        pcmax = 1srmax * uv; qcmax = pcmax / 3;
        pgmax = 0.5srmax * uv; qgmax = pgmax / 3;
        WLC[LOADS] = 100; WSD[LOADS] = 1000; WMG[MGL] = 400; WEG[EG] = 200;

        vmmin = vmin; vmmax = vmax; vgmin = 0.9uv; vgmax = 1.1uv; vcmin = 0.9uv; vcmax = 1.1uv;
        betamax = 0.2uv; betamin = 1 - betamax;
        pdmax = srmax * uv; qdmax = pdmax/3;
    elseif N == 6
        par[5] = 2;
        xrratio = 2; rc = 0.03; srmax = 0.1; Cloadc = 100;
        DG = copy(nodes);
        pcmax = 1.5srmax * uv; qcmax = pcmax / 3;
        pgmax = 0.5srmax * uv; qgmax = pgmax / 3; SDI = copy(nodes); noSDI =[];
        WLC[LOADS] = 100; WSD[LOADS] = 1000;
        MGL = copy(nodes);
        WMG[MGL] = 400; WEG[EG] = 200;
        vmmin = vmin; vmmax = vmax; vgmin = 0.95uv; vgmax = 1.1uv; vcmin = 0.9uv; vcmax = 1.1uv;
        LLCmax = sum(WSD); betamax = 0.2uv; betamin = 1 - betamax;
        pdmax = srmax * uv; qdmax = pdmax/3;
    elseif N == 12
        par[9] = 4;
        vmmin = vmin; vmmax = vmax; vgmin = 0.92uv; vgmax = 1.05uv; vcmin = 0.9uv; vcmax = 1.1uv; betamax = 0.2uv; betamin = 1 - betamax;
        xrratio = 2; rc = 0.01; srmax = 1/N; Cloadc = 100;
        WLC[LOADS] = 100; WSD[LOADS] = 1000;
        WMG[MGL] = 400; WEG[EG] = 200;

        if hasDroop
            srmax = 6/N; # for sequential vs online
            DG = nodes[nodes.%2 .== 1];
            # LOADS = setdiff(nodes, DG)
            pcmax[LOADS] = 1.25srmax; qcmax = pcmax / 3;
            pgmax[DG] = srmax; qgmax = pgmax / 3;
            SDI = copy(nodes); noSDI = [];
            RES = [2,7,11]; EG = RES;

            semax = zeros(Float64, N);
            # semax[EG] = 0.5(sum(pcmax))/length(EG);
            semax[RES] = 0.8srmax;

        elseif !diversification
            srmax = 6/N; # for sequential vs online
            LOADS = setdiff(nodes, DG)
            pcmax[LOADS] = 1.4srmax; qcmax = pcmax / 3;
            pgmax[DG] = srmax; qgmax = pgmax / 3;
            SDI = copy(nodes); noSDI = [];

            # ldg = length(DG); srand(167);
            # RES = sort((DG[randperm(ldg)])[1:round(Int64,ldg/2)]);
            RES = [2,7,11];
        else # to show diversification
            srmax = 3/N; # for sequential vs online
            # DG = copy(nodes);
            # LOADS = setdiff(nodes,DG)
            DG = copy(nodes); LOADS = copy(nodes);
            pcmax[LOADS] = 2.2srmax; qcmax = pcmax / 3;
            pgmax[DG] = srmax; qgmax = pgmax / 3;
            SDI = copy(nodes); noSDI = [];
            RES = copy(nodes)
        end
    elseif N == 24
        par[13] = 2; par[20] = 6; par[17] = 3;
        xrratio = 2; rc = 0.01; srmax = 6/N; Cloadc = 100;
        # DG = copy(nodes);
        LOADS = setdiff(nodes, DG);
        pcmax[LOADS] = 1.25srmax; qcmax = pcmax / 3;
        pgmax[DG] = srmax; qgmax = pgmax / 3;
        WLC[LOADS] = 100; WSD[LOADS] = 1000; MGL = [5,6,10];
        WMG[MGL] = 400; WEG[EG] = 200;
        betamax = 0.2uv; betamin = 1 - betamax;
        pdmax = 3srmax * uv; qdmax = pdmax/3; RES = copy(nodes);
        vmmin = vmin; vmmax = vmax; vgmin = 0.92uv; vgmax = 1.05uv; vcmin = 0.9uv; vcmax = 1.1uv;
    elseif N == 36
        par[5] = 3; par[6] = 2; par[10] = 7; par[13] = 10; par[15] = 13; par[16] = 2; par[20] = 18; par[21] = 16; par[24] = 22; par[25] = 22; par[27] = 25; par[31] = 29; par[32] = 28; par[36] = 34;
        xrratio = 2; rc = 0.01; srmax = 6/N; Cloadc = 100;
        # DG = copy(nodes);
        LOADS = setdiff(nodes, DG);
        pcmax[LOADS] = 1.25srmax; qcmax = pcmax / 3;
        pgmax[DG] = srmax; qgmax = pgmax / 3;
        WLC[LOADS] = 100; WSD[LOADS] = 1000; MGL = [5,6,10,13,15,25,32];
        WMG[MGL] = 400; WEG[EG] = 200;
        betamax = 0.2uv; betamin = 1 - betamax;
        pdmax = 3srmax * uv; qdmax = pdmax/3; RES = copy(nodes);
        vmmin = vmin; vmmax = vmax; vgmin = 0.92uv; vgmax = 1.05uv; vcmin = 0.9uv; vcmax = 1.1uv;
    elseif N == 69
        par[28] = 3; par[36] = 3; par[47] = 4; par[45] = 8; par[53] = 9; par[66] = 11; par[68] = 12;
        xrratio = 2; rc = 0.01; srmax = 6/N; Cloadc = 100; gammac = 1; lcfmin = 0.8ones(N,1);
        LOADS = setdiff(nodes, DG);

        pcmax[LOADS] = 1.25srmax; qcmax = pcmax / 3;
        pgmax[DG] = srmax; qgmax = pgmax / 3;
        WLC[LOADS] = 100; WSD[LOADS] = 1000;
        WMG[MGL] = 400; WEG[EG] = 200;
        vmmin = vmin; vmmax = vmax; vgmin = 0.92uv; vgmax = 1.05uv; vcmin = 0.9uv; vcmax = 1.1uv;

        betamax = 0.2uv; betamin = 1 - betamax;
        pdmax = 2srmax * uv; qdmax = pdmax/3;
    elseif N == 118
        # missing nodes are 28, 39, 57, 92, 104
        par[4] = 2; par[10] = 2; par[18] = 11; par[28] = 4; par[36] = 30; par[38] = 29; par[47] = 35; par[55] = 29; par[63] = 1; par[78] = 64; par[89] = 65; par[96] = 91; par[85] = 79; par[100] = 1; par[114] = 100;
        xrratio = 2; rc = 0.01; srmax = 6/N; Cloadc = 100; gammac = 1; lcfmin = 0.8ones(N,1);
        LOADS = setdiff(nodes, DG);
        pcmax[LOADS] = 1.25srmax; qcmax = pcmax / 3;
        pgmax[DG] = srmax; qgmax = pgmax / 3;
        WLC[LOADS] = 100; WSD[LOADS] = 1000;
        MGL = [28,36,47,45,53,66];
        WMG[MGL] = 400; WEG[EG] = 200;
        vmmin = vmin; vmmax = vmax; vgmin = 0.92uv; vgmax = 1.1uv; vcmin = 0.90uv; vcmax = 1.1uv;
        betamax = 0.2uv; betamin = 1 - betamax;
        pdmax = 2srmax * uv; qdmax = pdmax/3; RES = copy(nodes);
    elseif N == 124
        df = CSV.readtable("N123.csv");
        resistance = df[:R]; reactance = df[:X]; P = (df[:P])[1:85]; Q = (df[:Q])[1:85];
        loadNodes = (df[:Node])[1:85];
        loadNodes = map(loadNodes->parse(Int64,loadNodes),loadNodes)
        originalNodes = zeros(1:450);
        nodesInDoc = df[:nodeNoInDoc]; nodesInJulia = df[:nodeNoInJulia];
        to = df[:to]; from = df[:from];
        for i = 1:124
            j = nodesInDoc[i]; corr = nodesInJulia[i]; originalNodes[j] = corr;
        end
        nNodes = 124;
        N = nNodes; nodes = collect(1:nNodes); par = zeros(Int64,nNodes);

        for i = 1:122
            fromNode = from[i]; toNode = to[i];
            ofromNode = originalNodes[fromNode];
            otoNode = originalNodes[toNode];
            par[otoNode] = ofromNode;
        end

        pcmax = zeros(Float64,N); qcmax = zeros(Float64,N);
        pcmax[originalNodes[loadNodes]] = P;
        qcmax[originalNodes[loadNodes]] = Q;

        xrratio = 2; rc = 0.1; srmax = 1/N; Cloadc = 100; gammac = 1; lcfmin = 0.8ones(N,1);
        # pcmax = 2.5srmax * uv; qcmax = pcmax / 3;
        pgmax = 0.3pcmax; qgmax = 0.3qgmax;
        WLC[LOADS] = 100; WSD[LOADS] = 1000;
        MGL = [28,36,47,45,53,66];
        LTC = [13,33,73];
        WMG[MGL] = 400; WEG[EG] = 200;
        betamax = 0.2uv; betamin = 1 - betamax;
        pdmax = 2srmax * uv; qdmax = pdmax/3; RES = copy(MGL);
    end

    # RES = copy(MGL);
    srand(167); ldg = round(Int64,length(DG)/2);
    # RES = sort((DG[randperm(ldg)])[1:round(Int64,ldg/2)]);
    println("RES ", RES);

    # EG = copy(MGL);
    EG = RES;
    println("MGL ", MGL);
    println("*******************************************************")

    noMGL = setdiff(nodes, MGL); noLTC = setdiff(nodes, LTC);
    LLCmax = sum(WSD); # vmmin = vmin; vmmax = vmax; vgmin = 0.92uv; vgmax = 1.05uv; vcmin = 0.9uv; vcmax = 1.1uv;

    WAC = 10*ones(Float64,N);
    # WAC[1] = 0.1WLC[1];

    UDG = copy(DG);

    pgnom = 0.5pgmax; qgnom = 0.5qgmax; mq = 0.1;
    # semax[MGL] = 3semax[MGL];

    noMGL = setdiff(nodes, MGL);
    # println("Reached 2");
    # bilevel(N);
    leaves = getLeaves(par);
    CHILD = getChildrenMatrix(par);
    SuccMat = getSuccessorMatrix();

    r = rc * ones(Float64,N); x = xrratio * copy(r);

    if N != 124
        resistance = copy(r); reactance = copy(x);
    end
    resistance=resistance/2; reactance = reactance/2;

    path = getPath(par);
    commonPath, R, X = getCommonPathVariables(N,par,resistance,reactance);

    # println("R :", R);
    # println("X :", X);
    # println("size of R ", size(R));
    # N
end


function modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr)
    sruv = srmaxp * uv;
    pcmax = pc2sr * sruv; pgmax = pg2sr * sruv; pdmax = sd2sr * sruv;
    # semax[EG] = (sum(pcmax) + sum(pdmax) -
    # srmax * length(RES) - sum(pgmax))/length(MGL);
    # semax[MGL] = 3semax[MGL];
    semax[EG] = se2sr * sruv[EG];
    # semax[MGL] = 3semax[MGL];

    qcmax, qgmax, qdmax = pcmax/3, pgmax/3, pdmax/3;
end

function addObjectiveMicrogridModel(sm, t, beta, kcval, kmval, keval, pr, P)
    # println("Reached objective:")
    @objective(sm, :Min,
    # sum(WAC .* pr) +
    sum(WAC) * P[1] +
    WVR * sum(t) / N +
    sum(kcval[LOADS] .* (WSD-WLC)[LOADS]) +
    sum((1-beta)[LOADS] .* WLC[LOADS]) +
    # sum(keval[EG] .* WEG[EG]) +
    sum(kmval[MGL] .* WMG[MGL]));
end

function addMGBasicInequalities(sm, pr, qr, pg, qg, beta, v, t, kcval, kgval)
    @constraint(sm, betalb, beta[LOADS] .>= (1-kcval)[LOADS] .* betamin[LOADS]);
    @constraint(sm, betaub, beta[LOADS] .<= (1-kcval)[LOADS]);
    @constraint(sm, lovrMin, t .>= v0nom * uv - v);
    @constraint(sm, lovrMax, t .>= v - v0nom * uv);

    basicIneq = [betalb; betaub; lovrMin; lovrMax];
    basicIneqb = [((1-kcval) .* betamin)[LOADS]; (1-kcval)[LOADS]; v0nom * uv; -v0nom * uv];

    if length(RES) > 0
        @constraint(sm, prlb, pg[RES] .>= 0);
        @constraint(sm, prub, pg[RES] .<= semax[RES]);
        @constraint(sm, qrlb, qg[RES] .>= -semax[RES]);
        @constraint(sm, qrub, qg[RES] .<= semax[RES]);
        @constraint(sm, prqrdg1, 4pg[RES] + 3qg[RES] .<= 5semax[RES]);
        @constraint(sm, prqrdg2, 4pg[RES] - 3qg[RES] .<= 5semax[RES]);

        basicIneq = vcat(basicIneq, [prlb; prub; qrlb; qrub; prqrdg1; prqrdg2]);
        basicIneqb = vcat(basicIneqb, [ov[RES]; semax[RES]; -semax[RES]; semax[RES]; 5semax[RES]; 5semax[RES]]);
    end

    sm, basicIneq, basicIneqb
end

# Emergency generator is same as the storage devices.
function addEGInequalities(sm, pe, qe, keval)
    @constraint(sm, pelb, pe[EG] .>= 0);
    @constraint(sm, peub, pe[EG] .<= (keval.*semax)[EG]);
    @constraint(sm, qelb, qe[EG] .<= (keval.*semax)[EG]);
    @constraint(sm, qeub, qe[EG] .>= -(keval.*semax)[EG]);
    @constraint(sm, peqelb, (4pe + 3qe)[EG] .<= 5(keval.*semax)[EG]);
    @constraint(sm, peqeub, (4pe - 3qe)[EG] .<= 5(keval.*semax)[EG]);

    egIneq = [pelb; peub; qelb; qeub; peqelb; peqeub];
    egIneqb = [ov[EG]; (keval.*semax)[EG]; (keval.*semax)[EG];
    -(keval.*semax)[EG]; 5(keval.*semax)[EG]; 5(keval.*semax)[EG]];

    sm, egIneq, egIneqb
end

function addMGInequalities(sm, pg, qg, pe, qe, P, Q, v, vpar, kmval, keval)
    @constraint(sm, realCapUB, P[MGL] .<= 2(1-kmval)[MGL]);
    @constraint(sm, realCapLB, P[MGL] .>= -2(1-kmval)[MGL]);
    @constraint(sm, reacCapUB, Q[MGL] .<= 2(1-kmval)[MGL]);
    @constraint(sm, reacCapLB, Q[MGL] .>= -2(1-kmval)[MGL]);

    # @constraint(sm, voltMGUB, v[MGL] .<= v0nom + 2(1-kmval)[MGL]);
    # @constraint(sm, voltMGLB, v[MGL] .>= v0nom - 2(1-kmval)[MGL]);

    # Voltage droop equations
    @constraint(sm, voltMGUB, v[RES] .<= (vref - mq * (qg + qe - qgnom)[RES]) + 2(1-keval)[RES]);
    @constraint(sm, voltMGLB, v[RES] .>= (vref - mq * (qg + qe - qgnom)[RES]) - 2(1-keval)[RES]);

    @constraint(sm, dropMGLUB, v[MGL] - (vpar[MGL]
    - resistance[MGL].*P[MGL] - reactance[MGL].*Q[MGL]) .<= 2kmval[MGL]);
    @constraint(sm, dropMGLLB, v[MGL] - (vpar[MGL] -
    resistance[MGL].*P[MGL] - reactance[MGL].*Q[MGL]) .>= -2kmval[MGL]);

    mglineq = [realCapUB; realCapLB; reacCapUB; reacCapLB;
    voltMGUB;
    voltMGLB;
    dropMGLUB; dropMGLLB];
    mglineqb = [2(1-kmval)[MGL]; -2(1-kmval)[MGL]; 2(1-kmval)[MGL]; -2(1-kmval)[MGL];
    (vref - mq * (-qgnom)[RES]) +  2(1-keval)[RES];
    (vref - mq * (-qgnom)[RES]) - 2(1-keval)[RES];
    2kmval[MGL]; -2kmval[MGL]];

    sm, mglineq, mglineqb
end

function addMGAttackEqualities(sm, beta, v0, p, q, pr, qr, pe, qe, P, Q, kgval, kmval, delta, deltaV0)
    @constraint(sm, subVolt, v0 == v0nom + v0dis * deltaV0);

    deltaeq = [subVolt];
    deltab =  [v0nom];
    deltaMul = [v0dis];
    # deltaVecVar = [deltaV0Var; deltaVar; deltaVar];

    sm, deltaeq, deltab, deltaMul, subVolt
end

function addMGNoAttackEqualities(sm, beta, v, v0, vpar, p, q, pr, qr, pg, qg, pe, qe, P, Q, kgval, kmval, delta, deltaV0)

    @constraint(sm, realFlow, P .== CHILD * P + p);
    @constraint(sm, reacFlow, Q .== CHILD * Q + q);
    @constraint(sm, realCons, p - beta.*pcmax + pr + pe + pg .== 0);
    @constraint(sm, reacCons, q - beta.*qcmax + qr + qe + qg .== 0);

    noLM = setdiff(noMGL, LTC);
    @constraint(sm, dropNoLM, v[noLM] - (vpar - resistance.*P - reactance.*Q)[noLM] .== 0);
    @constraint(sm, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);

    mgleq = [realFlow; reacFlow; realCons; reacCons; dropNoLM; voltDropLTC];
    mgleqb = [ov; ov; ov; ov; ov[noLM]; ov[LTC]];

    sm, mgleq, mgleqb
end

function addMGOperatingConstraints(sm, kmval, keval)
    # @constraint(sm, (SuccMat' * kmval)[EG] .>= keval[EG]);
    for j in EG
        kmlines = 0;
        for i in MGL
            if SuccMat[i,j] == 1
                kmlines += kmval[i];
                @constraint(sm, kmval[i] <= keval[j]);
            end
        end
        # println(kmlines)
        kmlines != 0 ? @constraint(sm, kmlines >= keval[j]): nothing;
    end
end

function addMGConstraints(sm, p, pr, q, qr, pg, qg, pe, qe, beta, P, Q, v, t, v0, xvar, vpar, kcval, kgval, kmval, keval, delta, deltaV0)

    sm, basicIneq, basicIneqb = addMGBasicInequalities(sm, pr, qr, pg, qg, beta, v, t, kcval, kgval);
    sm, egIneq, egIneqb = addEGInequalities(sm, pe, qe, keval);
    sm, mglineq, mglineqb = addMGInequalities(sm, pg, qg, pe, qe, P, Q, v, vpar, kmval, keval);
    sm, deltaineq, deltaineqb, dgAttack = addAttackInequalities(sm, pr, qr, kgval, delta, deltaV0);

    sm, mgleq, mgleqb = addMGNoAttackEqualities(sm, beta, v, v0, vpar, p, q, pr, qr, pg, qg, pe, qe, P, Q, kgval, kmval, delta, deltaV0);
    sm, deltaeq, deltab, deltaMul, subVolt = addMGAttackEqualities(sm, beta, v0, p, q, pr, qr, pe, qe, P, Q, kgval, kmval, delta, deltaV0);
    sm, disineq, disb = addDisconnectivityConstraints(sm, kcval, kgval, v);
    # println("Added MGConstraints")

    ineq = [basicIneq; egIneq; mglineq; disineq; deltaineq]; ineqb = [basicIneqb; egIneqb; mglineqb; disb; deltaineqb];
    eq = [mgleq; deltaeq]; eqb = [mgleqb; deltab];

    ineq, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack
end

function getBasicMicrogridModel()
    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = getBasicSlaveModel();
    @variable(sm, pevar[1:length(EG)]); @variable(sm, qevar[1:length(EG)]);
    pe = zeros(AffExpr, N); pe[EG] = pevar;
    qe = zeros(AffExpr, N); qe[EG] = qevar;

    sm, p, pr, q, qr, pg, qg, pe, qe, beta, P, Q, v, t, v0, xvar, vpar
end

function getSlaveModelWithoutMicrogrid(kcval, kgval, kmval, keval, delta, deltaV0)
    global isCascadeMIP = false;
    sm, p, pr, q, qr, pg, qg, pe, qe, beta, P, Q, v, t, v0, xvar, vpar = getBasicMicrogridModel();

    addObjectiveMicrogridModel(sm, t, beta, kcval, kmval, keval, pr, P);

    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack = addMGConstraints(sm, p, pr, q, qr, pg, qg, pe, qe, beta, P, Q, v, t, v0, xvar, vpar, kcval, kgval, kmval, keval, delta, deltaV0);

    sm, p, pr, q, qr, beta, P, Q, v, t, v0, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack
end

function getSlaveModelWithMicrogrid(delta,deltaV0)
    global isCascadeMIP = true;

    sm, p, pr, q, qr, pg, qg, pe, qe, beta, P, Q, v, t, v0, xvar, vpar = getBasicMicrogridModel();

    @assert(length(UDG) > 0)
    @variable(sm, kcvar[1:length(LOADS)], Bin); @variable(sm, kgvar[1:length(UDG)], Bin); @variable(sm, kmvar[1:length(MGL)], Bin); @variable(sm, kevar[1:length(EG)], Bin);

    kc = zeros(AffExpr, N); kc[LOADS] = kcvar;
    kg = zeros(AffExpr, N); kg[UDG] = kgvar;
    km = zeros(AffExpr, N); km[MGL] = kmvar;
    ke = zeros(AffExpr, N); ke[EG] = kevar;

    addObjectiveMicrogridModel(sm, t, beta, kc, km, ke, pr, P);

    ineq, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack = addMGConstraints(sm, p, pr, q, qr, pg, qg, pe, qe, beta, P, Q, v, t, v0, xvar, vpar, kc, kg, km, ke, delta, deltaV0);

    addMGOperatingConstraints(sm, km, ke);
    println("addMGOperatingConstraints")

    sm, p, pr, q, qr, pg, qg, pe, qe, beta, v, v0, t, P, Q, kc, kg, km, ke, ineq, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack
end

function getMicrogridValues(sm, sstatus, p, pr, q, qr, pg, qg, pe, qe, beta, v, v0, t,
    P, Q, kc, kg, km, ke)

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, pr, q, qr, pg, qg, beta, P, Q, v, t, kc, kg, vpar);

    kmv = getvalue(km); v0v = getvalue(v0); pev = getvalue(pe); qev = getvalue(qe); kev = getvalue(ke);

    lmgv = sum(WMG[MGL] .* kmv[MGL]); lmgv = round(lmgv, 3);
    legv = sum(WEG[EG] .* kev[EG]); legv = round(legv, 3);
    kmv = round.(Int64, kmv); kev = round.(Int64, kev);

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv
end

function evaluateSlaveModelWithMicrogrid(delta,deltaV0)
    sm, p, pr, q, qr, pg, qg, pe, qe, beta, v, v0, t, P, Q, kc, kg, km, ke, ineq, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt  = getSlaveModelWithMicrogrid(delta, deltaV0);
    println(sm)
    sstatus = solve(sm);
    if sstatus != :Optimal
        println("evaluateSlaveModelWithMicrogrid: solve status ", sstatus);
    end

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv = getMicrogridValues(sm, sstatus, p, pr, q, qr, pg, qg, pe, qe, beta, v, v0, t, P, Q, kc, kg, km, ke);

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv
end

function getBruteForceMicrogridAttack(M)

    deltaV0 = 0; delta = zeros(Int64,N);
    sm, p, pr, q, qr, pg, qg, pe, qe, beta, v, v0, t, P, Q, kc, kg, km, ke, ineq, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithMicrogrid(delta, deltaV0);

    bestobjv = -infty; bestDeltaV0 = 0; bestDelta = zeros(Int64,N);
    delta = zeros(Int64,N); deltaV0 = 0;

    for curdeltaV0 = 0:0
        JuMP.setRHS(subVolt, v0nom + v0dis * curdeltaV0);

        for attackSet in collect(combinations(UDG,M))
            delta[:] = 0; delta[attackSet] = 1;
            # for j = 1:length(UDG)
            #     i = UDG[j];
            #     JuMP.setRHS(dgAttack[i], -curdelta[i]);
            # end
            #
            # status = solve(sm);
            # if status != :Optimal
            #     println("Not solved to optimality. Check for attacked nodes");
            #     println(attackSet);
            #     continue;
            # end
            pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv = evaluateSlaveModelWithMicrogrid(delta, deltaV0)
            if objv >= bestobjv
                bestobjv = objv; bestDelta = copy(delta); bestDeltaV0 = deltaV0;
            end
        end
    end

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv = evaluateSlaveModelWithMicrogrid(bestDelta, bestDeltaV0);
    println("Brute force with microgrid results");
    println("Nodes attacked : ", nodes[bestDelta.>0]);
    println("LOADS disconnected : ", nodes[kcv .> 0.5]);
    println("DG disconnected : ", nodes[kgv .> 0.5]);
    println("LOADS shed : ", nodes[betav .< 1]);
    println("lacv, lovrv, lcontv, lsdv, lmgv : ", lovrv, " ", lcontv, " ", lsdv, " ", lmgv);
    println("LLCmax : ", LLCmax);
    println("****************************************");

    bestDelta, bestDeltaV0, pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv
end

function printMGResults(pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv)
    nr = 2;
    println("voltage : ", round.([v0v, vv],nr));
    println("t : ", round.(tv,nr));
    println("objv : ", round(objv,nr));

    println("Loads : ", LOADS, "\nshed : ", nodes[kcv.==1], "\nbetav : ", betav[LOADS], "\npcv : ", round.(pcv[LOADS],nr), "\nqcv : ", round.(qcv[LOADS],nr));
    println("********************************");

    println("DG : ", DG, "\nkgv : ", nodes[kgv.==1], "\nprv : ", round.(prv[DG], nr), "\nqrv : ", round.(qrv[DG], nr));
    println("********************************");

    println("RES : ", RES, "\nkev : ", nodes[kev.==1], "\npgv : ", round.(pgv[RES], nr), "\nqgv : ", round.(qgv[RES], nr));
    println("********************************");

    println("EG : ", EG, "\npev : ", round.(pev[EG], nr), "\nqev : ", round.(qev[EG], nr));
    println("********************************");

    println("MGL : ", MGL, "\nkmv : ", nodes[kmv.==1]);
    # println("dobjv : ", dobjv);

    println("********************************");
    println("Pv : ", round.(Pv,nr));
    println("Qv : ", round.(Qv,nr));

    println("pv : ", round.(pv,nr));
    println("qv : ", round.(qv,nr));
end

function getBendersMicrogridIteration(mm, deltaVar, deltaV0Var, verbose, LLCreq)
    hasLLCExceeded = false; isBendersCutAdded = false;
    isMasterModelNotSolvable = false;

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, kev, lacv, lacv, lovrv, lcontv, lsdv, lmgv, legv, delta, deltaV0 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;


    # println("Iteration ", j);
    if verbose
        println(mm);
    end

    mstatus = solve(mm);
    if mstatus != :Optimal
        println("Master model solve status : ", mstatus);
        isMasterModelNotSolvable = true;
        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv,
        qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv, delta, deltaV0
    end

    delta = round.(Int64,getvalue(deltaVar)); deltaV0 = getvalue(deltaV0Var);
    nNodesAttacked = sum(round.(Int64,delta));
    # println("delta : ", delta);
    if verbose
        println("no. of nodes attacked : ", nNodesAttacked);
        println("delta : ", nodes[delta .== 1]);
    end

    sm, p, pr, q, qr, pg, qg, pe, qe, beta, v, v0, t, P, Q, kc, kg, km, ke, ineq, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithMicrogrid(delta, deltaV0);

    if verbose
        println("Updated slave model");
        println(sm);
    end

    sstatus = solve(sm);
    if sstatus != :Optimal
        println("Slave model solve status : ", sstatus);
        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv,
        qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv, delta, deltaV0;
    end

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objvMIP, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv = getMicrogridValues(sm, sstatus, p, pr, q, qr, pg, qg, pe, qe, beta, v, v0, t, P, Q, kc, kg, km, ke);

    verbose ? printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv) : nothing;

    # curLLC = sum((WSD-WLC).*kcv) + sum(WLC.*(1-betav))
    # curLLC = objvMIP;
    if objvMIP >= LLCreq
        println("LLCreq threshold exceeded : ", objvMIP);
        hasLLCExceeded = true;
    end

    # println("vv : ", vv);
    if verbose
        println("microgrids formed : ", nodes[kmv .> 0.5]);
        println("LOADS disconnected : ", nodes[kcv .> 0.5]);
        println("DG disconnected : ", nodes[kgv .> 0.5]);
        println("LOADS shed : ", nodes[betav .< 1], "\ntotal load shed : ", sum(WLC.*(1-betav)));
        println("LLC : ", curLLC);
        println("LLCreq : ", LLCreq);
        println("objv : ", objv);
    end

    sm, p, pr, q, qr, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltab, deltaMul, subVolt, deltaineq, deltaineqb, dgAttack = getSlaveModelWithoutMicrogrid(kcv, kgv, kmv, kev, delta, deltaV0);

    if verbose
        println("****************************************");
        println("Printing slave model without microgrid");
        println(sm);
    end

    sstatus = solve(sm);
    objv = getobjectivevalue(sm);

    ineqlbar = getdual(ineq);
    eqnubar = getdual(eq);
    deltaIneqlbar = getdual(deltaineq);
    deltaeqnubar = getdual(deltaeq);

    deltaIneqVar = [deltaVar[UDG]; deltaVar[UDG]; deltaVar[UDG]];
    global eta
    deltaIneqMul = [ov[UDG]; -eta * pgmax[UDG]; -eta .* qgmax[UDG]];
    # deltaIneqMul = [ov[UDG]; -(1 - kgv[UDG]) .* pgmax[UDG]; -(1 - kgv[UDG]) .* qgmax[UDG]];
    # deltaIneqVar = [deltaVar[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];
    # deltaIneqMul = [uv[UDG]; ov[UDG]; ov[UDG]; ov[UDG]];
    # println("deltaIneqlbar: ",deltaIneqlbar);

    deltaeqVar = [deltaV0Var];
    deltaeqVarMul = [v0dis];

    if verbose
        println("lengths :", length(deltaeqnubar), " ", length(deltaeqVarMul), " ", length(deltaeqVar));
        println("lengths :", length(deltaIneqlbar), " ", length(deltaIneqMul), " ", length(deltaIneqVar));
        println(length(eqnubar), " ", length(eqb), " ", length(ineqlbar), " ",
        length(ineqb));

        println("deltaeqnubar : ", round.(deltaeqnubar,3));
        println("deltaIneqlbar : ", round.(deltaIneqlbar,3));
        println("vv : ", vv);
        println("lacv, lovrv, lcontv, lsdv: ", lovrv, ", ", lcontv, ", ", lsdv)
        println("sumprod of nu and rhs : ", sum(eqnubar.*eqb));
        println("sumprod of lambda and rhs : ", sum(ineqlbar.*ineqb));
        println("deltaeqnubar: ",deltaeqnubar);
        println("deltaIneqlbar: ",deltaIneqlbar);
        println("sum(WLC) : ", sum(WLC));
        println("sum((WSD-WLC).*kcv) : ", sum((WSD-WLC).*kcv));
        println("sum of lhs : ", sum((WSD-WLC).*kcv) + sum(WLC) +sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb))
        println("objv : ", objv);
    end

    if maximum(abs.(deltaeqnubar)) + maximum(abs.(deltaIneqlbar)) >= myinf
        bendersCutLHS = sum(WMG[MGL].*kmv[MGL]) +
        # sum(WEG[EG].*kev[EG]) +
        sum((WSD-WLC).*kcv) + sum(WLC[LOADS]) + sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb) + sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        if !trial
            bendersCutLHS -= (objv + epsilon);
        else
            dconst = (sum((WSD-WLC).*kcv) + sum(WLC[LOADS]));
            diff = objv - dconst;
            eta = 0.98; rhs = epsilon;
            # rhs += dconst + max(diff, eta * diff + (1-eta)*300);
            rhs += objv - 5 * sum(deltaeqnubar .* deltaeqVarMul);
            # rhs += objv;
            bendersCutLHS -= rhs;
        end
        @constraint(mm, bendersCutLHS >= 0);
        # println("objv, lacv, lovrv, lcontv, lsdv : ", round.([objv, lacv, lovrv, lcontv, lsdv], 2));
        # println("nodes attacked: ", nodes[delta.== 1]);
        # println("bendersCut-objv: ", (bendersCutLHS));
        # @constraint(mm, bendersCutLHS >= objv + epsilon);

        isBendersCutAdded = true;
    else
        println("Not a good cut ", maximum(abs(eqnubar)));
        println("eqnubar: ", eqnubar);
        printResults(pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv);
    end

    hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv, delta, deltaV0
end

function getBendersMethodMicrogrid(LLCreq)
    kcv = ov; kgv = ov; delta = ov; deltaV0 = 0;
    withCascade = true; verbose = false; printSummary = true;
    curLLC = 0; nNodesAttacked = zeros(Int64, N);

    mm, deltaVar = getMasterModel();
    @variable(mm, 0 <= deltaV0Var <= 1, Int);
    # @constraint(mm, deltaV0)

    jmax = 100;
    for j = 1:jmax

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv, deltav, deltaV0v = getBendersMicrogridIteration(mm, deltaVar, deltaV0Var, verbose, LLCreq);

        println("Iteration : ", j);
        println("hasLLCExceeded, isMasterModelNotSolvable, isBendersCutAdded : ",
        hasLLCExceeded, " ", isMasterModelNotSolvable, " ", isBendersCutAdded);
        if isBendersCutAdded || hasLLCExceeded
            deltaV0 = deltaV0v; delta = deltav;
        end

        if isMasterModelNotSolvable || !isBendersCutAdded || hasLLCExceeded
            break;
        end
    end
    println("Exited loop");

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv = evaluateSlaveModelWithMicrogrid(delta, deltaV0);

    if printSummary
        println("Benders microgrid method results");
        println("Substation voltage : ", v0v);
        println("Nodes attacked : ", nodes[delta.>0])
        println("Microgrids formed : ", nodes[kmv.>0])
        println("LOADS disconnected : ", nodes[kcv .> 0.5]);
        println("DG disconnected : ", nodes[kgv .> 0.5]);
        println("LOADS shed : ", nodes[betav .< 1]);
        println("objv : ", objv);
        println("LLCreq : ", LLCreq);

        println("****************************************");
    end

    llc = lcontv + lsdv + lmgv + legv;
    delta, pv, prv, qv, qrv, betav, vv, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv, lmgv, legv
end

function getBendersMicrogridMinNodesLLC()
    kcval = ov; kgval = ov; kmval = zeros(Int64, length(MGL));
    delta = ov; deltaV0 = 0;
    withCascade = true; verbose = false; curLLC = 0;
    minNodesv = zeros(Int64,1); sumkcvs = zeros(Int64,1); sumkgvs = zeros(Int64,1);
    sumkmvs = zeros(Int64,1);
    minDeltav = zeros(Int64,N);
    objvs = []; tempLLCreq = 0;

    # Master model begins
    mm, deltaVar = getMasterModel();
    @variable(mm, -1 <= deltaV0Var <= 1, Int);

    jmax = N > 12? 15 : 100;
    for j = 1:jmax
        println("Iteration ", j);
        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv,
        qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv,
        qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv, deltav,
        deltaV0v = getBendersMicrogridIteration(mm, deltaVar, deltaV0Var, verbose,
        LLCmax);

        tempLLCreq = objv;

        if !isMasterModelNotSolvable
            delta = deltav; deltaV0 = deltaV0v;
        end
        if !isBendersCutAdded || isMasterModelNotSolvable
            break;
        end

        objvs = vcat(objvs, objv/LLCmax);

        if j > 1
            minDeltav = hcat(minDeltav, delta); sumkcvs = vcat(sumkcvs, sum(kcv));
            sumkgvs = vcat(sumkgvs, sum(kgv));
            nNodesAttacked = sum(delta);
            minNodesv = vcat(minNodesv,nNodesAttacked);
            sumkmvs = vcat(sumkmvs, sum(kmv));
        end

    end
    println("Exited loop");

    minNodesv, minDeltav, objvs, sumkcvs, sumkgvs, sumkmvs
end

function getBendersCascadeForFixedM(M)
    kcval = ov; kgval = ov; delta = ov; deltaV0 = 0;
    withCascade = true; verbose = false; curLLC = 0;
    nNodesAttacked = zeros(Int64, N); printSummary = false;

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    @constraint(mm, sum(deltaVar) == M);
    # @constraint(mm, deltaV0)

    println("Entering loop");
    jmax = N > 12? 15 : 100;
    for j = 1:jmax

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv,
        qrv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv,
        lcontv, lsdv, deltav, deltaV0v = getBendersCascadeIteration(mm, deltaVar,
        deltaV0Var, verbose, LLCreq);

        println("Iteration : ", j);
        # println("hasLLCExceeded, isMasterModelNotSolvable, isBendersCutAdded : ",
        # hasLLCExceeded, " ", isMasterModelNotSolvable, " ", isBendersCutAdded);
        if isBendersCutAdded || hasLLCExceeded
            deltaV0 = deltaV0v; delta = deltav;
        end

        if isMasterModelNotSolvable || !isBendersCutAdded || hasLLCExceeded
            break;
        end

    end
    println("Exited loop");

    pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);

    if printSummary
        println("Benders cascade method results");
        println("Nodes attacked : ", nodes[delta.>0])
        println("LOADS disconnected : ", nodes[kcv .> 0.5]);
        println("DG disconnected : ", nodes[kgv .> 0.5]);
        println("LOADS shed : ", nodes[betav .< 1]);
        println("LLC : ", objv);
        println("LLCreq : ", LLCreq);
        println("****************************************");
    end

    llc = objv;
    delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    objv, pgv, qgv, lacv, lovrv, lcontv, lsdv
end

function getBendersMicrogridForFixedM(M)
    kcv = ov; kgv = ov; delta = ov; deltaV0 = 0;
    withCascade = true; verbose = false; printSummary = true;
    curLLC = 0; nNodesAttacked = zeros(Int64, N);

    mm, deltaVar = getMasterModel();
    @variable(mm, -1 <= deltaV0Var <= 1, Int);
    @constraint(mm, sum(deltaVar) == M);

    jmax = 20;
    for j = 1:jmax

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, prv, qv,
        qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, pgv,
        qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv, legv, deltav,
        deltaV0v = getBendersMicrogridIteration(mm, deltaVar, deltaV0Var, verbose,
        LLCreq);

        println("Iteration : ", j);
        println("hasLLCExceeded, isMasterModelNotSolvable, isBendersCutAdded : ",
        hasLLCExceeded, " ", isMasterModelNotSolvable, " ", isBendersCutAdded);
        if isBendersCutAdded || hasLLCExceeded
            deltaV0 = deltaV0v; delta = deltav;
        end

        if isMasterModelNotSolvable || !isBendersCutAdded || hasLLCExceeded
            break;
        end
    end
    println("Exited loop");

    pv, prv, qv, qrv, pev, qev, betav, vv, v0v, tv, kcv, kgv, kmv, kev, pcv,
    qcv, objv, pgv, qgv, Pv, Qv, lacv, lovrv, lcontv, lsdv, lmgv,
    legv = evaluateSlaveModelWithMicrogrid(delta, deltaV0);

    if printSummary
        println("Benders microgrid method results");
        println("Substation voltage : ", v0v);
        println("Nodes attacked : ", nodes[delta.>0])
        println("Microgrids formed : ", nodes[kmv.>0])
        println("LOADS disconnected : ", nodes[kcv .> 0.5]);
        println("DG disconnected : ", nodes[kgv .> 0.5]);
        println("LOADS shed : ", nodes[betav .< 1]);
        println("objv : ", objv);
        println("LLCreq : ", LLCreq);

        println("****************************************");
    end

    llc = lcontv + lsdv + lmgv + legv;
    delta, pv, prv, qv, qrv, betav, vv, tv, kcv, kgv, kmv, kev, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv, lmgv, legv
end

function getUniqueFirstSortSecondArrays(a, b)

    na = zeros(Int64,length(unique(a)));
    nb = zeros(Float64, length(na));
    nind = zeros(Int64,length(na));

    i = 1; j = 1; cind = 0;
    while i <= length(a)
        # println("i,j:",i, " ", j);
        n = a[i];
        nbest = b[i]; cind = i;
        i = i + 1;
        while i <= length(a) && a[i] == n
            if b[i] >= nbest
                nbest = b[i]; cind = i;
            end
            i = i + 1;
        end

        na[j] = n; nb[j] = nbest; nind[j] = cind; j = j + 1;
    end
    # println(na, " ", nb);

    na, nb, nind
end

function getRecoveryModel(delta, deltaV0)
    sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));

    @variable(sm, pr[RES,TS]); @variable(sm, qr[RES,TS]);
    @variable(sm, pg[DG,TS]); @variable(sm, qg[DG,TS]);
    @variable(sm, kg[DG,TS]); @variable(sm, yg[DG,TS]);

    @variable(sm, beta[LOADS,TS]); @variable(sm, kc[LOADS,TS]); @variable(sm, pc[LOADS,TS]); @variable(sm, qc[LOADS,TS]);
    @variable(sm, km[MGL,TS]);


    @variable(sm, P[nodes,TS]); @variable(sm, Q[nodes,TS]);
    @variable(sm, v0[TS]); @variable(sm, v[nodes,TS]);
    # @variable(sm, t[1:N]);
    @variable(sm, t[nodes,TS]);

    vpar = zeros(AffExpr,N,T);

    p = zeros(AffExpr,N,T); q = zeros(AffExpr,N,T);

    for i in nodes, ts in TS
        vpar[i,ts] = par[i] == 0 ? v0[ts] : v[par[i],ts];
    end

    ### Load model
    # load control
    @constraint(sm, [i in LOADS, ts in TS], beta[i,ts] >= (1 - kc[i,ts]) * betamin[i]);
    @constraint(sm, [i in LOADS, ts in TS], beta[i,ts] <= 1 - kc[i,ts]);
    @constraint(sm, [i in LOADS, ts in TS], kc[i,ts] >= vcmin[i] - v[i,ts]);
    @constraint(sm, [i in LOADS, ts in TS], kc[i,ts] >= v[i,ts] - vcmax[i]);

    ### DG model
    # @constraint(sm, [i in DG, ts in TS], pg[i,ts] >= 0);
    @constraint(sm, [i in DG, ts in TS], pg[i,ts] == (1-kg[i,ts]) * pgmax[i]);
    @constraint(sm, [i in DG, ts in TS], qg[i,ts] == (1-kg[i,ts]) * qgmax[i]);
    # @constraint(sm, [i in DG, ts in TS], qg[i,ts] >= -(1-kg[i,ts]) * qgmax[i]);

    @constraint(sm, [i in DG, ts in TS], kg[i,ts] >= vgmin[i] - v[i,ts]);
    @constraint(sm, [i in DG, ts in TS], kg[i,ts] >= v[i,ts] - vgmax[i]);

    ### Resource model
    @constraint(sm, [i in RES, ts in TS], pr[i,ts] >= 0);
    @constraint(sm, [i in RES, ts in TS], pr[i,ts] <= semax[i]);
    @constraint(sm, [i in RES, ts in TS], qr[i,ts] >= -semax[i]);
    @constraint(sm, [i in RES, ts in TS], qr[i,ts] <= semax[i]);
    @constraint(sm, [i in RES, ts in TS], 4pr[i,ts] + 3qr[i,ts] <= 5semax[i]);
    @constraint(sm, [i in RES, ts in TS], 4pr[i,ts] - 3qr[i,ts] <= 5semax[i]);

    # Voltage droop equations
    @constraint(sm, [i in RES, ts in TS], v[i,ts] <= (vref - mq * (qr[i,ts] - qgnom[i])) + 2(1-km[1,ts]));
    @constraint(sm, [i in RES, ts in TS], v[i,ts] >= (vref - mq * (qr[i,ts] - qgnom[i])) - 2(1-km[1,ts]));

    # @constraint(sm, voltMGUB, v[MGL] .<= v0nom + 2(1-kmval)[MGL]);
    # @constraint(sm, voltMGLB, v[MGL] .>= v0nom - 2(1-kmval)[MGL]);

    for i in nodes
        for ts in TS
            if i in LOADS
                p[i,ts] += beta[i,ts] * pcmax[i]; q[i,ts] += beta[i,ts] * qcmax[i];
            end
            if i in RES
                p[i,ts] -= pr[i,ts]; q[i,ts] -= qr[i,ts];
            end
            if i in DG
                p[i,ts] -= pg[i,ts]; q[i,ts] -= qg[i,ts];
            end
        end
    end

    ### Power flow constraints
    ## power conservation
    @constraint(sm, [i in nodes, ts in TS], P[i,ts] == sum(CHILD[i,j] * P[j,ts] for j in nodes) + p[i,ts]);
    @constraint(sm, [i in nodes, ts in TS], Q[i,ts] == sum(CHILD[i,j] * Q[j,ts] for j in nodes) + q[i,ts]);

    ## Voltage constraints
    noLM = setdiff(noMGL, LTC);
    @constraint(sm, [i in noLM, ts in TS], v[i,ts] == vpar[i,ts] - resistance[i] * P[i,ts] - reactance[i] * Q[i,ts]);
    @constraint(sm, [i in LTC, ts in TS], v[i,ts] == LTCSetting * vpar[i,ts]);

    ### microgrid model
    # Power flow across microgrid lines
    @constraint(sm, [i in MGL, ts in TS], P[i,ts] <= 10(1-km[i,ts]));
    @constraint(sm, [i in MGL, ts in TS], P[i,ts] >= -10(1-km[i,ts]));
    @constraint(sm, [i in MGL, ts in TS], Q[i,ts] <= 10(1-km[i,ts]));
    @constraint(sm, [i in MGL, ts in TS], Q[i,ts] >= -10(1-km[i,ts]));
    # Voltage across microgrid lines
    @constraint(sm, [i in MGL, ts in TS], v[i,ts] - (vpar[i,ts]
    - resistance[i] * P[i,ts] - reactance[i] * Q[i,ts]) <= 2km[i,ts]);
    @constraint(sm, [i in MGL, ts in TS], v[i,ts] - (vpar[i,ts]
    - resistance[i] * P[i,ts] - reactance[i] * Q[i,ts]) >= -2km[i,ts]);

    ### Initial conditions
    @constraint(sm, [i in DG, ts = 1], kg[i,ts] == delta[i]);


    ### Recovery conditions
    @constraint(sm, [i in DG, ts = 1], yg[i,ts] == 0);
    @constraint(sm, [i in DG, ts = 1:T-1], kg[i,ts+1] == kg[i,ts] - yg[i,ts+1]);
    @constraint(sm, [i = 1, ts = 1:T-1], km[i,ts] == 1);

    @constraint(sm, [ts = 2:T], sum(yg[i,ts] for i in DG) <= NG);
    @constraint(sm, [ts = 1:T-1], v0[ts] == v0nom + v0dis * deltaV0);

    ### Final conditions
    @constraint(sm, v0[T] == v0nom);
    @constraint(sm, [i = 1, ts = T], km[i,ts] == 0);

    ### Objective
    @constraint(sm, [i in nodes, ts in TS], t[i,ts] >= v0nom - v[i,ts]);
    @constraint(sm, [i in nodes, ts in TS], t[i,ts] >= v[i,ts] - v0nom);

    @variable(sm, Lac[TS]); @variable(sm, Lovr[TS]); @variable(sm, Lsd[TS]); @variable(sm, Lcont[TS]); @variable(sm, Lmg[TS]); @variable(sm, sysPerf[TS]); @variable(sm, sysCost[TS]);

    @constraint(sm, [ts in TS], Lac[ts] == sum(WAC) * P[1,ts]);
    @constraint(sm, [ts in TS], Lovr[ts] == WVR/N * sum(t[i,ts] for i in nodes));
    @constraint(sm, [ts in TS], Lsd[ts] == sum(kc[i,ts] * (WSD-WLC)[i] for i in LOADS));
    @constraint(sm, [ts in TS], Lcont[ts] == sum((1-beta[i,ts]) * WLC[i] for i in LOADS));
    @constraint(sm, [ts in TS], Lmg[ts] == sum(km[i,ts] * WMG[i] for i in MGL));

    @constraint(sm, [ts in TS], sysCost[ts] == Lac[ts] + Lovr[ts] + Lsd[ts] + Lcont[ts] + Lmg[ts]);
    @constraint(sm, [ts in TS], sysPerf[ts] == 100 * (1 - sysCost[ts]/LLCmax));

    @objective(sm, :Min, sum(sysCost[ts] for ts in TS));

    println(sm);

    sm, beta, kc, kg, pg, qg, pr, qr, p, q, km, P, Q, v, v0, t, Lac, Lovr, Lsd, Lcont, Lmg, sysCost, sysPerf
end

function getRecoveryVariableValues(sm, sstatus, beta, kc, kg, pg, qg, pr, qr, p, q, km, P, Q, v, v0, t, Lac, Lovr, Lsd, Lcont, Lmg, sysCost, sysPerf)
    assert(sstatus == :Optimal);
    pcv = zeros(Float64, length(LOADS), T); qcv = zeros(Float64, length(LOADS), T);

    betav = getvalue(beta); kcv = getvalue(kc);
    # println(JuMP.size(betav), " ", size(pcmax));
    for i in LOADS, ts in TS
        pcv[i,ts] = betav[i,ts] * pcmax[i]; qcv[i,ts] = betav[i,ts] * qcmax[i];
    end

    kgv = getvalue(kg); pgv = getvalue(pg); qgv = getvalue(qg);
    prv = getvalue(pr); qrv = getvalue(qr);
    pv = getvalue(p); qv = getvalue(q);
    kmv = getvalue(km);

    Pv = getvalue(P); Qv = getvalue(Q);
    vv = getvalue(v); v0v = getvalue(v0);
    tv = getvalue(t);

    lacv = getvalue(Lac); lovrv = getvalue(Lovr); lsdv = getvalue(Lsd); lcontv = getvalue(Lcont); sysCostv = getvalue(sysCost); sysPerfv = getvalue(sysPerf);

    objv = getobjectivevalue(sm);

    # betav, kcv, pcv, qcv, kgv, pgv, qgv, prv, qrv, pv, qv, kmv, Pv, Qv, vv, v0v, tv, objv, lacv, lovrv, lcontv, lsdv = round.([betav, kcv, pcv, qcv, kgv, pgv, qgv, prv, qrv, pv, qv, kmv, Pv, Qv, vv, v0v, tv, objv, lacv, lovrv, lcontv, lsdv],3);

    # kcv = round.(Int64, kcv); kgv = round.(Int64,kgv); kmv = round.(Int64, kmv);

    betav, kcv, pcv, qcv, kgv, pgv, qgv, prv, qrv, pv, qv, kmv, Pv, Qv, vv, v0v, tv, objv, lacv, lovrv, lcontv, lsdv, sysCostv, sysPerfv
end

function evaluateRecoveryModel(delta, deltaV0)
    sm, beta, kc, kg, pg, qg, pr, qr, p, q, km, P, Q, v, v0, t, Lac, Lovr, Lsd, Lcont, Lmg, sysCost, sysPerf = getRecoveryModel(delta, deltaV0);
    sstatus = solve(sm);

    betav, kcv, pcv, qcv, kgv, pgv, qgv, prv, qrv, pv, qv, kmv, Pv, Qv, vv, v0v, tv, objv, lacv, lovrv, lcontv, lsdv, sysCostv, sysPerfv = getRecoveryVariableValues(sm, sstatus, beta, kc, kg, pg, qg, pr, qr, p, q, km, P, Q, v, v0, t, Lac, Lovr, Lsd, Lcont, Lmg, sysCost, sysPerf)
end
