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
a = zeros(Int64,0); a = [a;3]
function getAncestorSuccessorLoadPairs()
    ancestorLoads = zeros(Int64,0);
    succesorLoads = zeros(Int64,0);
    for i in LOADS
        for j in LOADS
            if i == j
                continue
            end
            if SuccMat[i,j] == 1
                ancestorLoads = [ancestorLoads;i]
                succesorLoads = [succesorLoads;j]
            end
        end
    end

    ancestorLoads, succesorLoads
end

function getAncestorSuccessorDGPairs()
    ancestorDGs = zeros(Int64,0);
    succesorDGs = zeros(Int64,0);
    for i in DG
        for j in DG
            if i == j
                continue
            end
            if SuccMat[i,j] == 1
                ancestorDGs = [ancestorDGs;i]
                succesorDGs = [succesorDGs;j]
            end
        end
    end

    ancestorDGs, succesorDGs
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
    @variable(sm, pgvar[1:length(DG)]); @variable(sm, qgvar[1:length(DG)]);

    @variable(sm, betavar[1:length(LOADS)]);

    #voltages are squared
    @variable(sm, v[1:N]); @variable(sm, v0);
    # @variable(sm, t[1:N]);
    @variable(sm, tvar);
    @variable(sm, P[1:N]); @variable(sm, Q[1:N]);


    beta = zeros(AffExpr,N); t = zeros(AffExpr,N); pg = zeros(AffExpr,N); qg = zeros(AffExpr,N);

    beta[LOADS] = betavar; vpar = zeros(AffExpr,N); t[:] = tvar; pg[DG] = pgvar; qg[DG] = qgvar;
    for i = 1:N
        vpar[i] = par[i] == 0 ? v0 : v[par[i]];
    end

    xvar = [p; q; pg; qg; beta; P; Q; v; v0; t];

    p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar
end

function getBasicSlaveModel()
    sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));

    p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = addVariablesToBasicSlaveModel(sm);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar
end

function getBasicSlaveSocpModel(withCascade=true)
    # sm = Model(solver=GurobiSolver(OutputFlag = 0, genv, Method=2,Cuts=2,Presolve=2,MIPFocus=1));
    sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));
    if withCascade == false
        sm = Model(solver=MosekSolver(QUIET=true));
    end
    # sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));
    p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = addVariablesToBasicSlaveModel(sm);

    @variable(sm, ell[1:N]);
    xvar = vcat(xvar, ell);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell
end



function getVariableValues(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, vpar)
    vv = zeros(Float64, length(v)); pv = zeros(Float64, length(p));
    qv = zeros(Float64, length(q));
    betav = zeros(Float64, length(beta));
    tv = 0; kcv = zeros(Bool, length(kc)); kgv = zeros(Bool, length(kg));
    pcv = zeros(Float64, length(beta)); qcv = zeros(Float64, length(beta));
    pgv = zeros(Float64, length(pgmax)); qgv = zeros(Float64, length(qgmax));
    Pv = zeros(Float64, N); Qv = zeros(Float64, N);
    objv = -infty; lacv, lovrv, lcontv, lsdv = 0,0,0,0;

    if sstatus != :Optimal
        println("Slave model solve status : ", sstatus);
        return pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
    end

    vv = getvalue(v); pv = getvalue(p); qv = getvalue(q);
    tv = getvalue(t);
    betav = getvalue(beta); pcv = pcmax .* betav; qcv = qcmax .* betav;
    Pv = getvalue(P); Qv = getvalue(Q);
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

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getVariableValuesSocp(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, ell, vpar)
    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, vpar);

    ellv = zeros(Float64, length(ell));
    if sstatus == :Optimal
        ellv = getvalue(ell)
    end
    lacv = sum(WAC) * sum(resistance.*ellv);

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function addObjectiveDisconnectModel(sm, t, beta, kcval)
    @objective(sm, :Min,
    WVR * sum(t) / N +
    sum(WLC[LOADS].*(1-beta[LOADS])) +
    sum((WSD[LOADS]-WLC[LOADS]).*kcval[LOADS]));
end

function addObjectiveDisconnectSocpModel(sm, t, beta, kcval, ell)
    @objective(sm, :Min,
    WVR * sum(t) / N +
    sum(WLC[LOADS].*(1-beta[LOADS])) +
    sum((WSD[LOADS]-WLC[LOADS]).*kcval[LOADS]) +
    sum(WAC)*sum(resistance.*ell));
end


function addBasicInequalities(sm, pg, qg, beta, v, t, kcval, kgval)
    @constraint(sm, betalb, beta[LOADS] .>= (1-kcval[LOADS]) .* betamin[LOADS]);
    @constraint(sm, betaub, beta[LOADS] .<= 1-kcval[LOADS]);

    @constraint(sm, lovrMin, t .>= v0nom * uv - v);
    @constraint(sm, lovrMax, t .>= v - v0nom * uv);

    basicIneq = [betalb; betaub; lovrMin; lovrMax];
    basicIneqb = [(1-kcval[LOADS]) .* betamin[LOADS]; 1-kcval[LOADS]; v0nom * uv; -v0nom * uv];

    sm, basicIneq, basicIneqb
end

function addBasicEqualities(sm, v, vpar, p, q, pg, qg, beta, P, Q, kgval)
    @constraint(sm, voltDrop, (vpar - v - 2resistance.*P - 2reactance.*Q)[noLTC] .== 0)
    @constraint(sm, realFlow, P .== SuccMat * p);
    @constraint(sm, reacFlow, Q .== SuccMat * q);
    @constraint(sm, realCons, p - beta.*pcmax + pg .== 0);
    @constraint(sm, reacCons, q - beta.*qcmax + qg .== 0);

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

function addBasicSocpEqualities(sm, v, vpar, p, q, pg, qg, beta, P, Q, kgval, ell)
    @constraint(sm, voltDrop, (vpar - v - 2resistance.*P - 2reactance.*Q + (resistance.^2+reactance.^2).*ell)[noLTC] .== 0)
    @constraint(sm, realFlow, P .== CHILD * P + p + resistance.*ell);
    @constraint(sm, reacFlow, Q .== CHILD * Q + q + reactance.*ell);

    @constraint(sm, realCons, p - beta.*pcmax + pg .== 0);
    @constraint(sm, reacCons, q - beta.*qcmax + qg .== 0);

    basicSocpEq = [voltDrop; realFlow; reacFlow; realCons; reacCons];
    basicSocpEqb = [ov[noLTC]; ov; ov; ov; ov];
    if length(LTC) > 0
        @constraint(sm, voltDropLTC, v[LTC] - LTCSetting * vpar[LTC] .== 0);
        basicSocpEq = vcat(basicSocpEq, voltDropLTC);
        basicSocpEqb = vcat(basicSocpEqb, ov[LTC]);
    end

    sm, basicSocpEq, basicSocpEqb
end

function addBasicSocpInequalities(sm, pg, qg, beta, v, t, kcval, kgval, vpar, P, Q, ell)
    sm, basicIneq, basicIneqb = addBasicInequalities(sm, pg, qg, beta, v, t, kcval, kgval);

    @constraint(sm, socpIneq[i=1:N], norm([vpar[i], ell[i], sqrt(2)*P[i], sqrt(2)*Q[i]]) <= (vpar[i] + ell[i]));

    socpIneqb = ov;

    sm,  basicIneq, basicIneqb
end

function addAttackInequalities(sm, pg, qg, kgval, delta, deltaV0)

    # @constraint(sm, dgAttack, -kgval[UDG] .<= -delta[UDG]);
    # @variable(sm, kgconvar[1:length(UDG)]); kgcon = zeros(AffExpr,N); kgcon[UDG] = kgconvar;

    # @constraint(sm, dgAttack0, kgcon[UDG] .== delta[UDG]);
    @constraint(sm, dgAttack1, kgval[UDG] .>= delta[UDG]);
    # global isCascadeMIP
    # if isCascadeMIP
    #     global eta = 0;
    # else
    #     global eta = 10epsilon;
    # end

    @constraint(sm, dgAttack2a, pg[UDG] .<= (1 - delta[UDG]) .* pgmax[UDG]);
    @constraint(sm, dgAttack2b, qg[UDG] .<= (1 - delta[UDG]) .* qgmax[UDG]);
    @constraint(sm, dgAttack3a, pg[UDG] .<= (1 - kgval[UDG]) .* pgmax[UDG]);
    @constraint(sm, dgAttack3b, pg[UDG] .>= (1 - kgval[UDG]) .* pgmax[UDG]);
    @constraint(sm, dgAttack4a, qg[UDG] .<= (1 - kgval[UDG]) .* qgmax[UDG]);
    @constraint(sm, dgAttack4b, qg[UDG] .>= (1 - kgval[UDG]) .* qgmax[UDG]);

    dgAttack = [dgAttack1; dgAttack2a; dgAttack2b; dgAttack3a; dgAttack3b; dgAttack4a; dgAttack4b];
    # dgAttack = dgAttack1;
    deltaineq = dgAttack;
    # deltaineqb = ov[UDG];
    # deltaineqb = [ov[UDG]; pgmax[UDG]; qgmax[UDG]];
    deltaineqb = [delta[UDG]; pgmax[UDG]; qgmax[UDG]; (1 - kgval[UDG]) .* pgmax[UDG]; (1 - kgval[UDG]) .* pgmax[UDG]; (1 - kgval[UDG]) .* qgmax[UDG]; (1 - kgval[UDG]) .* qgmax[UDG]];

    sm, deltaineq, deltaineqb, dgAttack, dgAttack1,  dgAttack2a, dgAttack2b, dgAttack3a, dgAttack3b, dgAttack4a, dgAttack4b
end

function addAttackEqualities(sm, v0, delta, deltaV0)
    # println("addAttackEqualities: v0nom, deltaV0 : ", v0nom, ", ", deltaV0);
    # @constraint(sm, subvolt, v0 == v0nom + v0dis * deltaV0 - vreg * (Q[1] - Q1nom));
    @constraint(sm, subvolt, v0 == v0nom + v0dis * deltaV0);

    deltaeq = [subvolt];
    deltaeqb = [v0nom];

    sm, deltaeq, deltaeqb, subvolt
end

function addDisconnectivityConstraints(sm, kc, kg, v)
    # println("lengths: ", length(kc)," ",length(vcmin)," ",length(v))
    @constraint(sm, loadLow, v[LOADS] .>= vcmin[LOADS] - kc[LOADS]);
    @constraint(sm, loadUp, v[LOADS] .<= vcmax[LOADS] + kc[LOADS]);

    disineq = [loadLow; loadUp]; disb = [vcmin[LOADS] - kc[LOADS]; vcmax[LOADS] + kc[LOADS]];
    # if length(UDG) > 0
    @constraint(sm, dgLow, v[UDG] .>= vgmin[UDG] - kg[UDG]);
    @constraint(sm, dgUp, v[UDG] .<= vgmax[UDG] + kg[UDG]);
    disineq = vcat(disineq, dgLow, dgUp); disb = vcat(disb, vgmin[UDG] - v[UDG], vgmax[UDG] + kg[UDG]);
    # end

    # LnDG = intersect(LOADS, DG);
    # if lndgDisconnect && length(LnDG) > 0
    #     vcminLessVgmin = vcmin[LnDG] .<= vgmin[LnDG];
    #     length(LnDG) > 0 ? @constraint(sm, vcminLessVgmin .* (kc - kg)[LnDG] .<= 0): nothing;
    # end

    sm, disineq, disb
end

function modifyBinaryAttackConstraints(delta, deltaV0, deltaineq, modifyinglp = true)
    nudg = length(UDG);

    rhsb = [delta[UDG]; (1 - delta[UDG]) .* pgmax[UDG]; (1 - delta[UDG]) .* qgmax[UDG]]; # only need to change rhs for operator MIP
    # println("nudg = ",nudg)
    istart = modifyinglp ? nudg + 1 : 1;
    # println("istart = ",istart)
    # println(nodes[delta.==1])
    for i = istart:3nudg
        JuMP.setRHS(deltaineq[i], rhsb[i])
    end
    deltaineq
end

function modifyAncestorSuccessorDGCuts(dgAncestorSuccessorCuts, delta)
    for i = 1:length(dgAncestorSuccessorCuts)
        j = ancestorDGs[i]; k = succesorDGs[i];
        JuMP.setRHS(dgAncestorSuccessorCuts[i], delta[j]+delta[k])
    end
    dgAncestorSuccessorCuts
end

function modifyAllBinaryOperatorConstraints(sm, kcval, kgval, delta, deltaV0, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta)
    # println("modifyAllBinaryOperatorConstraints\nkgval = ", kgval)
    # addAttackInequalities
    modifyinglp = true;
    deltaineq = modifyBinaryAttackConstraints(delta, deltaV0, deltaineq, modifyinglp);

    #modifying constraints for operator LP
    nudg = length(UDG);
    deltaineqb = [ov[UDG]; pgmax[UDG]; qgmax[UDG]; (1 - kgval[UDG]) .* pgmax[UDG]; (1 - kgval[UDG]) .* pgmax[UDG]; (1 - kgval[UDG]) .* qgmax[UDG]; (1 - kgval[UDG]) .* qgmax[UDG]];
    for i = 3nudg+1:7nudg # changing only kc/kg constraints
        JuMP.setRHS(deltaineq[i], deltaineqb[i])
    end

    # addAttackEqualities
    JuMP.setRHS(subvolt, v0nom + v0dis * deltaV0);
    deltaeqb = [v0nom + v0dis * deltaV0];

    # BasicInequalities
    nloads = length(LOADS);
    basicIneqb[1:2nloads] = [(1-kcval[LOADS]) .* betamin[LOADS]; 1-kcval[LOADS]];

    for i = 1:2nloads
        JuMP.setRHS(basicIneq[i], basicIneqb[i]);
    end

    # disconnect inequalities
    disb = [vcmin[LOADS] - kcval[LOADS]; vcmax[LOADS] + kcval[LOADS]; vgmin[UDG] - kgval[UDG]; vgmax[UDG] + kgval[UDG]];
    for i = 1:length(disb)
        JuMP.setRHS(disineq[i], disb[i]);
    end

    # modify all the rhs in (lambdaTimesb)
    ineqb = [basicIneqb; deltaineqb; disb]; eqb = [basicEqb; deltaeqb];

    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb
end

function modifyAllBinaryOperatorConstraintsForLinear(sm, kcval, kgval, delta, deltaV0, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta)
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb = modifyAllBinaryOperatorConstraints(sm, kcval, kgval, delta, deltaV0, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta)

    # Don't forget to modify the objective
    addObjectiveDisconnectModel(sm, t, beta, kcval);

    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb
end

function modifyAllBinaryOperatorConstraintsForSocp(sm, kcval, kgval, delta, deltaV0, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta, ell)
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb = modifyAllBinaryOperatorConstraints(sm, kcval, kgval, delta, deltaV0, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta)

    # Don't forget to modify the objective
    addObjectiveDisconnectSocpModel(sm, t, beta, kcval, ell);

    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb
end

function addCascadeConstraints(sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, kcv, kgv, delta, deltaV0)
    # Whenever you add constraints, add the update the list of constraints and right hand side ineq and ineqb

    sm, deltaineq, deltaineqb, dgAttack = addAttackInequalities(sm, pg, qg, kgv, delta, deltaV0);
    sm, deltaeq, deltaeqb, subvolt = addAttackEqualities(sm, v0, delta, deltaV0);
    sm, basicIneq, basicIneqb = addBasicInequalities(sm, pg, qg, beta, v, t, kcv, kgv);
    sm, basicEq, basicEqb = addBasicEqualities(sm, v, vpar, p, q, pg, qg, beta, P, Q, kgv);
    sm, disineq, disb = addDisconnectivityConstraints(sm, kcv, kgv, v);

    ineq = vcat(basicIneq, deltaineq, disineq); ineqb = vcat(basicIneqb, deltaineqb, disb);
    eq = [basicEq; deltaeq]; eqb = [basicEqb; deltaeqb];

    basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt
end

function addSocpCascadeConstraints(sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell, kcv, kgv, delta, deltaV0)
    # Whenever you add constraints, add the update the list of constraints and right hand side ineq and ineqb

    sm, deltaineq, deltaineqb, dgAttack = addAttackInequalities(sm, pg, qg, kgv, delta, deltaV0);
    sm, deltaeq, deltaeqb, subvolt = addAttackEqualities(sm, v0, delta, deltaV0);
    sm,  basicIneq, basicIneqb = addBasicSocpInequalities(sm, pg, qg, beta, v, t, kcv, kgv, vpar, P, Q, ell);
    sm, basicEq, basicEqb = addBasicSocpEqualities(sm, v, vpar, p, q, pg, qg, beta, P, Q, kgv, ell);
    sm, disineq, disb = addDisconnectivityConstraints(sm, kcv, kgv, v);

    ineq = vcat(basicIneq, deltaineq, disineq); ineqb = vcat(basicIneqb, deltaineqb, disb);
    eq = [basicEq; deltaeq]; eqb = [basicEqb; deltaeqb];

    basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt
end

function addSpecialOperatorBinaryConstraints(sm, kc, kg, delta)

    # println("kc = ",kc)
    # println("ancestorLoads = ", ancestorLoads);
    # println("successorLoads = ", succesorLoads);
    @constraint(sm, loadAncestorSuccessorCuts[i = 1:length(ancestorLoads)], kc[ancestorLoads[i]] <= kc[succesorLoads[i]]);

    @constraint(sm, dgAncestorSuccessorCuts[i = 1:length(ancestorDGs)], kg[ancestorDGs[i]] <= kg[succesorDGs[i]] + delta[ancestorDGs[i]] + delta[succesorDGs[i]]); # change the right hand side to delta[i] + delta[j] implying that the condition does not hold when either of the two DGs is attacked

    loadAncestorSuccessorCuts, dgAncestorSuccessorCuts
end
# Given: attack vector, connectivity of LOADS and DERs that remain unchanged.
# Output: slave model which just an LP.
function getSlaveModelWithoutCascade(kcv, kgv, delta, deltaV0)
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = getBasicSlaveModel();

    addObjectiveDisconnectModel(sm, t, beta, kcv);
    global isCascadeMIP = false;
    basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = addCascadeConstraints(sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, kcv, kgv, delta, deltaV0);

    # println(sm);
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt
end

function getSlaveSocpModelWithoutCascade(kcv, kgv, delta, deltaV0)
    withCascade = false;
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell = getBasicSlaveSocpModel(withCascade);

    addObjectiveDisconnectSocpModel(sm, t, beta, kcv, ell);

    global isCascadeMIP = false;
    basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = addSocpCascadeConstraints(sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell, kcv, kgv, delta, deltaV0);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell,  basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt
end


function getSlaveModelWithCascade(delta, deltaV0, withCuts = false)
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar = getBasicSlaveModel();

    @variable(sm, kcvar[1:length(LOADS)], Bin); # kc
    @variable(sm, kgvar[1:length(DG)], Bin); # kg
    kc = zeros(AffExpr,N); kg = zeros(AffExpr,N);
    kc[LOADS] = kcvar; kg[DG] = kgvar;

    addObjectiveDisconnectModel(sm, t, beta, kc);
    global isCascadeMIP = true;
    basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = addCascadeConstraints(sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, kc, kg, delta, deltaV0);

    loadAncestorSuccessorCuts, dgAncestorSuccessorCuts = nothing, nothing
    if withCuts
        loadAncestorSuccessorCuts, dgAncestorSuccessorCuts = addSpecialOperatorBinaryConstraints(sm, kc, kg, delta)
    end

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts
end

function getSlaveSocpModelWithCascade(delta, deltaV0, withCuts = false)
    withCascade=true;
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell = getBasicSlaveSocpModel(withCascade);

    @variable(sm, kcvar[1:length(LOADS)], Bin); # kc
    @variable(sm, kgvar[1:length(DG)], Bin); # kg
    kc = zeros(AffExpr,N); kg = zeros(AffExpr,N);
    kc[LOADS] = kcvar; kg[DG] = kgvar;

    addObjectiveDisconnectSocpModel(sm, t, beta, kc, ell);

    global isCascadeMIP = true;
    basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = addSocpCascadeConstraints(sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell, kc, kg, delta, deltaV0);

    loadAncestorSuccessorCuts, dgAncestorSuccessorCuts = nothing, nothing
    if withCuts
        loadAncestorSuccessorCuts, dgAncestorSuccessorCuts = addSpecialOperatorBinaryConstraints(sm, kc, kg, delta)
    end

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts, ell
end

function evaluateSlaveModelWithCascade(delta, deltaV0, withCuts = false)
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts = getSlaveModelWithCascade(delta, deltaV0, withCuts);
    # println(sm);
    sstatus = solve(sm);

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, vpar);
    # printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function evaluateSlaveSocpModelWithCascade(delta, deltaV0, withCuts = false)
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts, ell = getSlaveSocpModelWithCascade(delta, deltaV0, withCuts);

    # println(sm);
    sstatus = solve(sm);

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValuesSocp(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, ell, vpar)
    # printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getEasyBenchmarkValues()
    betav = ones(Float64,N); kcv = zeros(Int64,N); kgv = copy(kcv);
    pcv = pcmax; qcv = qcmax;
    pgv = pgmax; qgv = qgmax;

    pv = pcv - pgv;
    qv = qcv - qgv;
    vv = v0nom - R * pv - X * qv;
    tv = sum(abs(vv-v0nom));
    Pv = SuccMat * pv; Qv = SuccMat * qv;

    lacv = 0; lovrv = WVR * tv; lcontv = sum(WLC .* (1-kcv) .* (1-betav));
    lsdv = sum((WSD - WLC) .* kcv);
    objv = lacv + lovrv + lcontv + lsdv;

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

# The goal is to ensure that the voltages are within the bounds, but the loss
# of load control and load shedding is minimized as much as possible.
function getPrecontingencyValues()
    delta = zeros(Int64,N); deltaV0 = 0;
    # global betamin = uv;
    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveModelWithCascade(delta, deltaV0);
    # printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv,
    # objv, pgv, qgv)

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getPrecontingencyValuesSecondMethod()
    delta = zeros(Int64,N); deltaV0 = 0;

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv
end

function getOnlineCascade(delta, deltaV0)

    noAttackDelta = zeros(Int64,N); noDeltaV0 = 0
    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveModelWithCascade(noAttackDelta, noDeltaV0);
    # println("vv: ",vv);
    # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    # lacv, lovrv, lcontv, lsdv = getPrecontingencyValues();
    # println("Loads disconnected: ", nodes[kcv.==1]);
    # println("DGs disconnected: ", nodes[kgv.==1]);
    # println("delta: ",delta)

    kgo = zeros(Int64,N); kgo[UDG] = delta[UDG];
    kco = kcv; myeps = 0.001;

    for i = 1:2N
        pcv = (1-kco) .* betav .* pcmax; qcv = (1-kco) .* betav .* qcmax;
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

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end


function getOnlineCascadeSocp(delta, deltaV0)

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

    sstatus = solve(sm);
    if sstatus != :Optimal
        println("Model solve status : ", sstatus);
        return;
    end
    kginv = getvalue(kg);
    vinv = getvalue(v);

    @constraint(sm, kgintermediate[i=1:ludg], kg[UDG[i]] >= kginv[UDG[i]]);
    @constraint(sm, loadVoltInterMin[i=1:lloads], kc[LOADS[i]] >= (vcmin[LOADS[i]] - vinv[LOADS[i]]));
    @constraint(sm, loadVoltInterMax[i=1:lloads], kc[LOADS[i]] >= (vinv[LOADS[i]] - vcmax[LOADS[i]]));

    # println("After cascade\n",sm);
    # println("vv: ",vv);
    # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    # lacv, lovrv, lcontv, lsdv = getPrecontingencyValues();
    # println("Loads disconnected: ", nodes[kcv.==1]);
    # println("DGs disconnected: ", nodes[kgv.==1]);
    # println("delta: ",delta)
    sstatus = solve(sm);
    if sstatus != :Optimal
        println("Model solve status : ", sstatus);
        return;
    end

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, vpar);

    # println("vv : ", round(vv,2))
    # println("loads shed : ", nodes[kcnew.==1]);
    # println("DGs shed : ", nodes[kgnew.==1]);

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getDGUnaffected()
    delta = zeros(Int64,N); deltaV0 = 0; delta[DG] = 1;
    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    # println("vv: ",vv);
    # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
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
    bestRandomDelta = zeros(Int64,N,1); deltaV0 = 0;
    bestkco = 0; bestkgo = 0; bestvv = 0; bestpgv = 0; bestqgv = 0; bestbetav = 0;
    for i = 1:nAttacks
        delta = zeros(Int64,N);
        delta[randperm(N)[1:M]] = 1;
        kco, kgo, vv, pgv, qgv, betav = getOnlineCascadeSocp(delta);
        LLCcur = objv;
        if LLCcur > LLCbest
            LLCbest = LLCcur;
            bestkcv, bestkgv, bestvv, bestpgv, bestqgv, bestbetav = kcv, kgv, vv, pgv, qgv, betav;
            bestRandomDelta = delta;
        end
    end

    bestkcv, bestkgv, bestvv, bestpgv, bestqgv, bestbetav
end

function getRandomAttackOnlineCascadeForSocp(M, nAttacks)
    LLCbest = -infty;
    bestRandomDelta = zeros(Int64,N,1);
    bestkco = 0; bestkgo = 0; bestvv = 0; bestpgv = 0; bestqgv = 0; bestbetav = 0;
    for i = 1:nAttacks
        delta = zeros(Int64,N);
        delta[randperm(N)[1:M]] = 1;

        pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getOnlineCascadeSocp(delta,deltaV0);
        LLCcur = sum(WSD .* kco + (1 - kco) .* WLC .* (1 - betav));
        if LLCcur > LLCbest
            LLCbest = LLCcur;
            bestkco, bestkgo, bestvv, bestpgv, bestqgv, bestbetav = kco, kgo, vv, pgv, qgv, betav;
            bestRandomDelta = delta;
        end
    end

    bestkco, bestkgo, bestvv, bestpgv, bestqgv, bestbetav, bestRandomDelta
end

# Brute force algorithm -- needs to be reworked
function getBruteForceCascadeAttack(M, withCuts = true)

    bestDelta = zeros(Int64,N); bestDeltaV0 = 0;
        sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts = getSlaveModelWithCascade(bestDelta, bestDeltaV0);

    bestobjv = -infty; bestDeltaV0 = 0; bestDelta = zeros(Int64,N);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts = getSlaveModelWithCascade(delta, deltaV0, withCuts)

    ludg = length(UDG);
    dgAttack1 = dgAttack[1:ludg]; dgAttack2 = dgAttack[ludg+1:2ludg]; dgAttack3 = dgAttack[2ludg+1:3ludg];
    bestObjv = -infty; nodes = collect(1:N);
    for curdeltaV0 = 0:0 #-1:1
        # JuMP.setRHS(subvolt, v0nom + v0dis * curdeltaV0);

        for attackSet in collect(combinations(UDG,M))
            delta = zeros(Int64,N); delta[attackSet] = 1;
            modifyinglp = false;
            deltaineq = modifyBinaryAttackConstraints(delta, deltaV0, deltaineq, modifyinglp);

            if withCuts
                dgAncestorSuccessorCuts = modifyAncestorSuccessorDGCuts(dgAncestorSuccessorCuts,delta);
            end

            # println(sm)
            solve(sm);

            # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, vpar);
            objv = getobjectivevalue(sm)
            if objv >= bestObjv
                bestObjv = objv; bestDelta = copy(delta); bestDeltaV0 = curdeltaV0;
            end
        end
    end

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveModelWithCascade(bestDelta, bestDeltaV0);
    println("Brute force results");
    println("Nodes attacked : ", nodes[bestDelta.>0])
    println("LOADS disconnected : ", nodes[kcv .> 0.5]);
    println("DG disconnected : ", nodes[kgv .> 0.5]);
    println("LOADS shed : ", nodes[betav .< 1]);
    println("LLC : ", sum((WSD[LOADS]-WLC[LOADS]) .* kcv[LOADS]) + sum(WLC[LOADS] .* (1-betav[LOADS])) + sum(WAC) * sum(resistance.*ellv) + WVR * sum(tv) / N);
    println("LLCmax : ", LLCmax);
    println("****************************************");

    bestDelta, bestDeltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getBruteForceCascadeSocpAttack(M, withCuts = true)

    bestDelta = zeros(Int64,N); bestDeltaV0 = 0;
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, kc, kg, xvar, vpar, ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts, ell = getSlaveSocpModelWithCascade(bestDelta, bestDeltaV0, withCuts);

    bestObjv = -infty; bestDeltaV0 = 0; bestDelta = zeros(Int64,N);

    nodes = collect(1:N);
    for curdeltaV0 = 0:0 #-1:1
        # JuMP.setRHS(subvolt, v0nom + v0dis * curdeltaV0);

        for attackSet in collect(combinations(UDG,M))
            delta = zeros(Int64,N); delta[attackSet] = 1;

            modifyinglp = false;
            deltaineq = modifyBinaryAttackConstraints(delta, deltaV0, deltaineq, modifyinglp);

            if withCuts
                dgAncestorSuccessorCuts = modifyAncestorSuccessorDGCuts(dgAncestorSuccessorCuts,delta);
            end

            # println(sm)
            solve(sm);

            # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValuesSocp(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, ell, vpar);
            objv = getobjectivevalue(sm)

            # println("Nodes attacked : ", nodes[delta.>0], ", objv = ",objv)
            if objv >= bestObjv
                bestObjv = objv; bestDelta = copy(delta); bestDeltaV0 = curdeltaV0;
            end
        end
    end

    # deltaineq = modifyBinaryAttackConstraints(delta, deltaV0, deltaineq, modifyinglp);
    #
    # if withCuts
    #     dgAncestorSuccessorCuts = modifyAncestorSuccessorDGCuts(dgAncestorSuccessorCuts,bestDelta);
    # end
    # println(sm)

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveSocpModelWithCascade(bestDelta, bestDeltaV0);
    println("Brute force results");
    println("Nodes attacked : ", nodes[bestDelta.>0])
    println("LOADS disconnected : ", nodes[kcv .> 0.5]);
    println("DG disconnected : ", nodes[kgv .> 0.5]);
    println("LOADS controlled : ", nodes[betav .< 1]);
    println("LLC : ", sum((WSD[LOADS]-WLC[LOADS]) .* kcv[LOADS]) + sum(WLC[LOADS] .* (1-betav[LOADS])) + sum(WAC) * sum(resistance.*ellv) + WVR * sum(tv) / N);
    println("objv : ", objv)
    println("LLCmax : ", LLCmax);
    println("****************************************");

    bestDelta, bestDeltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)
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
    println("qcv : ", round.(qcv,nr));
    println("qgv : ", round.(qgv,nr));

    println("pv : ", round.(pv,nr));
    println("qv : ", round.(qv,nr));

end

function printSocpResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, ellv)
    nr = 2;
    printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv);
    println("ellv : ", round.(ellv,nr));
end

global bestObjvUntilNow = 0;

function getBendersCascadeIteration(mm, deltaVar, deltaV0Var, verbose, LLCreq, smmip, kcmip, kgmip, deltaineqmip, sm, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta, dgAncestorSuccessorCuts, withCuts, epsilon = 0, withVariableEpsilon = false, nhighest = 0)
    hasLLCExceeded = false; isBendersCutAdded = false; isMasterModelNotSolvable = false;

    pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    # println(mm);
    if verbose
        println(mm);
    end

    mstatus = solve(mm);
    if mstatus != :Optimal
        println("Master model solve status : ", mstatus);
        isMasterModelNotSolvable = true;
        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv,
        kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0;
    end

    delta = round.(Int64,getvalue(deltaVar)); deltaV0 = getvalue(deltaV0Var);
    nNodesAttacked = sum(round.(Int64,delta));
    # println("delta : ", delta);
    # Update the right hand side constraints of the operator MIP
    deltaineqmip = modifyBinaryAttackConstraints(delta, deltaV0, deltaineqmip);
    if verbose
        println("no. of nodes attacked : ", nNodesAttacked);
        println("delta : ", nodes[delta .== 1]);
    end

    if verbose
        println("Updated slave model");
        println(smmip);
    end

    sstatus = solve(smmip);
    kcv = getvalue(kcmip);
    kgv = getvalue(kgmip);
    objvMIP = getobjectivevalue(smmip);

    if sstatus != :Optimal
        # println(sm)
        println("kgmipv = ", kgv)
        println("kcmipv = ", kcv)
        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv = false, false, false, objvMIP
        println("Slave model solve status : ", sstatus);
        if sstatus == :Stall
            isBendersCutAdded = true;
            @constraint(mm, sum(deltaVar[delta.==0]) + sum(1-deltaVar[delta.==1]) >= 1)
            println("Stalled status, but delta eliminated.")
        end

        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, delta, deltaV0;
    end

    # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objvMIP, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, vpar);



    # verbose ? printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv) : nothing;

    # curLLC = sum((WSD-WLC).*kcv) + sum(WLC.*(1-betav))
    curLLC = objvMIP;
    if curLLC >= LLCreq
        verbose? println("LLCreq threshold exceeded : ", curLLC) : nothing;
        hasLLCExceeded = true;
    end

    # Update the right hand side constraints of the operator LP
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb = modifyAllBinaryOperatorConstraintsForLinear(sm, kcv, kgv, delta, deltaV0, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta);

    if verbose
        println("kgmipv = ", kgv)
        println("kcmipv = ", kcv)
        println("****************************************");
        println("Printing slave model without cascade");
        println(sm);
    end

    sstatus = solve(sm);
    objv = getobjectivevalue(sm);

    ludg = length(UDG);
    # ineqlbar = getdual(ineq);
    # eqnubar = getdual(eq);
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

    ludg = length(UDG);
    row = 2
    deltaIneqlbar[(row-1)*ludg+1:(row)*ludg] +=  deltaIneqlbar[(2row-1)*ludg+1:(2row)*ludg];
    row = 3
    deltaIneqlbar[(row-1)*ludg+1:(row)*ludg] +=  deltaIneqlbar[(2row-1)*ludg+1:(2row)*ludg];

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
        # bendersCutLHS = sum(WLC[LOADS]) + sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb) + sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        bendersCutLHS = sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        if !trial
            if withVariableEpsilon
                dualValues = deltaIneqlbar[ludg+1:3ludg] .* deltaIneqMul[ludg+1:3ludg];
                actualDualValues = zeros(Float64,ludg);

                for i = 1:ludg
                    actualDualValues[i] = dualValues[i] + dualValues[ludg+i];
                end
                
                sortDualValues = sort(actualDualValues,rev=true);
                ktemp = max(sum(delta),1);
                seqEnd = min(ludg,nhighest+ktemp);
                seqStart = seqEnd - ktemp + 1;
                epsilon = sum(sortDualValues[seqStart:seqEnd]);
                epsilon -= 0.0001 # added to avoid numerical issues
                print("ktemp = ",ktemp,", epsilon = ",round(epsilon,1), " ");
            end
            bendersCutLHS -= (epsilon);
            # println("bendersCut\n",bendersCutLHS)
            # for i = 1:ktemp
            #     j = (i+1)%ludg +1
            #     epsilon += sortDualValues[j]
            # end


            # println("dualValues = ",dualValues)
            # println("sortDualValues = ",sortDualValues)

            # epsilon = 20;
            # if N >=36
            #     if sum(delta) >= 2 && sum(delta) <= 15
            #         epsilon = 25
            #     else
            #         epsilon = 10;
            #     end
            # end
            # bendersCutLHS -= (objv + epsilon);

            # bendersCutLHS -= (objv + epsilon);
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

        println("objv = ",objv)
        @constraint(mm, bendersCutLHS >= 0);
        @constraint(mm, sum(deltaVar[delta.==0]) + sum(1-deltaVar[delta.==1]) >= 1)
        # println("objv, lacv, lovrv, lcontv, lsdv : ", round.([objv, lacv, lovrv, lcontv, lsdv], 2));
        # println("nodes attacked: ", nodes[delta.== 1]);
        # println("bendersCut-objv: ", (bendersCutLHS));
        # @constraint(mm, bendersCutLHS >= objv + epsilon);
        # global bestObjvUntilNow
        global bestObjvUntilNow = max(bestObjvUntilNow, objv) + LLCmax / 1000;

        isBendersCutAdded = true;
    else
        println("Not a good cut ", maximum(abs(eqnubar)));
        println("eqnubar: ", eqnubar);
        printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv);
    end

    hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, delta, deltaV0
end

# global isKSpecified = true, M = 15;
function getBendersMethodCascade(LLCreq)
    kcval = ov; kgval = ov; delta = ov; deltaV0 = 0;  cardinality = -1;
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

    ## operator MIP subproblem
    smmip, pmip, qmip, pgmip, qgmip, betamip, Pmip, Qmip, vmip, tmip, v0mip, kcmip, kgmip, xvarmip, vparmip, ineqmip, ineqbmip, eqmip, eqbmip, deltaeqmip, deltaeqbmip, subvoltmip, deltaineqmip, deltaineqbmip = getSlaveModelWithCascade(delta, deltaV0);

    ## operator LP subproblem
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = getSlaveModelWithoutCascade(kcval, kgval, delta, deltaV0)

    println("Entering loop");
    jmax = N >= 36? 5 : N >= 24? 800 : N == 12 ? 100 : 100;
    lastChange = 0;
    j = 0; cumTime = 0; tic();
    prevObjv = 0;
    for j = 1:jmax
        # println(mm);
        nCurrentIter += 1;

        if sum(delta) > cardinality

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

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, deltav, deltaV0v = getBendersCascadeIteration(mm, deltaVar, deltaV0Var, verbose, LLCreq, smmip, kcmip, kgmip, deltaineqmip, sm, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta);

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
            resilience = floor(100(1-prevObjv/LLCmax),2);
            resiliences = [resiliences; resilience]
            println(resilience, ", (", j+1, "), ", round(cumTime,2), ", ", round(100cardinality/N,2));
            lastChange = nCurrentIter;
            break;
        end

        prevObjv = objv;
    end
    println("Exited loop");

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveModelWithCascade(delta, deltaV0);

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
        for ij = 2:length(nIters)
            println(ploss[ij], ", (",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], " \\\\");
        end
    end

    global nBendersCascadeIter = j;
    llc = objv;
    delta, deltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getBendersCascadeIterationForSocp(mm, deltaVar, deltaV0Var, verbose, LLCreq, smmip, kcmip, kgmip, deltaineqmip, sm, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta, ell, dgAncestorSuccessorCuts, withCuts, epsilon = 0, withVariableEpsilon = false, nhighest = 0)
    hasLLCExceeded = false; isBendersCutAdded = false; isMasterModelNotSolvable = false;

    pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0, ellv = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    # println(mm);
    if verbose
        println(mm);
    end

    mstatus = solve(mm);

    if mstatus != :Optimal
        println("Master model solve status : ", mstatus);
        # println(mm);
        isMasterModelNotSolvable = true;
        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv,
        kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0;
    end

    delta = round.(Int64,getvalue(deltaVar)); deltaV0 = getvalue(deltaV0Var);
    nNodesAttacked = sum(round.(Int64,delta));
    # println("delta : ", delta);
    # Update the right hand side constraints of the operator MIP
    deltaineqmip = modifyBinaryAttackConstraints(delta, deltaV0, deltaineqmip);
    if withCuts
        modifyAncestorSuccessorDGCuts(dgAncestorSuccessorCuts, delta)
    end

    if verbose
        println("no. of nodes attacked : ", nNodesAttacked);
        println("delta : ", nodes[delta .== 1]);
    end

    if verbose
        println("Updated slave model");
        println(smmip);
    end

    sstatus = solve(smmip);
    if sstatus != :Optimal
        println("Slave model solve status : ", sstatus);
        return;
    end

    # pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objvMIP, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getVariableValues(sm, sstatus, p, q, pg, qg, beta, P, Q, v, t, kc, kg, vpar);

    objvMIP = getobjectivevalue(smmip);
    kcv = getvalue(kcmip); kcv = max.(kcv,0);
    kgv = getvalue(kgmip); kgv = max.(kgv,0);

    # verbose ? printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv) : nothing;

    # curLLC = sum((WSD-WLC).*kcv) + sum(WLC.*(1-betav))
    curLLC = objvMIP;
    if curLLC >= LLCreq
        verbose? println("LLCreq threshold exceeded : ", curLLC) : nothing;
        hasLLCExceeded = true;
    end

    # Update the right hand side constraints of the operator LP
    ineq, ineqb, eq, eqb, deltaeq, deltaeqb, subvolt, deltaineq, deltaineqb = modifyAllBinaryOperatorConstraintsForSocp(sm, kcv, kgv, delta, deltaV0, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta, ell);

    if verbose
        println("kgmipv = ", kgv)
        println("kcmipv = ", kcv)
        println("****************************************");
        println("Printing slave model without cascade");
        println(sm);
    end

    sstatus = solve(sm);
    objv = getobjectivevalue(sm);
    if sstatus != :Optimal
        # println(sm)
        println("kgmipv = ", kgv)
        println("kcmipv = ", kcv)
        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv = false, false, false, objvMIP
        println("Slave model solve status : ", sstatus);
        if sstatus == :Stall
            isBendersCutAdded = true;
            @constraint(mm, sum(deltaVar[delta.==0]) + sum(1-deltaVar[delta.==1]) >= 1)
            println("Stalled status, but delta eliminated.")
        end

        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, delta, deltaV0;
    end

    # ineqlbar = getdual(ineq);
    # eqnubar = getdual(eq);
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
    #preprocessing of deltaIneqlbar
    ludg = length(UDG);
    row = 2
    deltaIneqlbar[(row-1)*ludg+1:(row)*ludg] +=  deltaIneqlbar[(2row-1)*ludg+1:(2row)*ludg];
    row = 3
    deltaIneqlbar[(row-1)*ludg+1:(row)*ludg] +=  deltaIneqlbar[(2row-1)*ludg+1:(2row)*ludg];

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
        # bendersCutLHS = sum(WLC[LOADS]) + sum(ineqlbar.*ineqb) + sum(eqnubar.*eqb) + sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        # bendersCutLHS = sum((WSD-WLC).*kcv) + sum(WLC[LOADS]) + sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        bendersCutLHS = sum(deltaIneqlbar .* deltaIneqMul .* deltaIneqVar) + sum(deltaeqnubar .* deltaeqVarMul .* deltaeqVar);
        if !trial
            if withVariableEpsilon
                dualValues = deltaIneqlbar[ludg+1:3ludg] .* deltaIneqMul[ludg+1:3ludg];
                actualDualValues = zeros(Float64,ludg);

                for i = 1:ludg
                    actualDualValues[i] = dualValues[i] + dualValues[ludg+i];
                end
                sortDualValues = sort(actualDualValues,rev=true);
                ktemp = max(sum(delta),1);
                seqEnd = min(ludg,nhighest+ktemp);
                seqStart = seqEnd - ktemp + 1;
                epsilon = sum(sortDualValues[seqStart:seqEnd]);
                epsilon -= 0.0001 # added to avoid numerical issues
                print("ktemp = ",ktemp,", epsilon = ",round(epsilon,1), " ");
            end
            # for i = 1:ktemp
            #     j = (i+1)%ludg +1
            #     epsilon += sortDualValues[j]
            # end


            # println("dualValues = ",dualValues)
            # println("sortDualValues = ",sortDualValues)

            # epsilon = 20;
            # if N >=36
            #     if sum(delta) >= 2 && sum(delta) <= 15
            #         epsilon = 25
            #     else
            #         epsilon = 10;
            #     end
            # end
            # bendersCutLHS -= (objv + epsilon);
            bendersCutLHS -= (epsilon);
            # bendersCutLHS -= LLCreq;
            # println("bendersCut\n",bendersCutLHS)
        else
            dconst = (sum((WSD-WLC).*kcv) + sum(WLC[LOADS]));
            diff = objv - dconst;
            eta = 0.98; rhs = epsilon;
            # rhs += dconst + max(diff, eta * diff + (1-eta)*300);
            rhs += objv - 5 * sum(deltaeqnubar .* deltaeqVarMul);
            # rhs += objv;
            bendersCutLHS -= rhs;
        end

        # println("objv = ",objv)
        @constraint(mm, bendersCutLHS >= 0);
        @constraint(mm, sum(deltaVar[delta.==0]) + sum(1-deltaVar[delta.==1]) >= 1)
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
        printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv);
    end

    hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, delta, deltaV0
end

# global isKSpecified = true, M = 15;
function getBendersMethodCascadeForSocp(LLCreq, epsilon, withVariableEpsilon, nhighest, withCuts=true)
    kcval = ov; kgval = ov; delta = ov; deltaV0 = 0;  cardinality = -1;
    withCascade = true; verbose = false; curLLC = 0;
    nNodesAttacked = zeros(Int64, N); printSummary = true;
    nIters = [0]; nSecs = [0]; nCurrentIter = 0; nCardinalities = [0]; objv = 0;
    resiliences = [100];
    # withCuts = false;
    global bestObjvUntilNow = 0;

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    # @constraint(mm, sum(deltaVar) >= kpercent/100 * length(UDG));
    # N >= 36 ? mm = addDownstreamHeuristicCuts(mm, deltaVar) : nothing;
    # @constraint(mm, deltaV0)
    # println(mm)

    ## operator MIP subproblem
    smmip, pmip, qmip, pgmip, qgmip, betamip, Pmip, Qmip, vmip, tmip, v0mip, kcmip, kgmip, xvarmip, vparmip, ineqmip, ineqbmip, eqmip, eqbmip, deltaeqmip, deltaeqbmip, subvoltmip, deltaineqmip, deltaineqbmip, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts, ellmip  = getSlaveSocpModelWithCascade(delta, deltaV0, withCuts);

    ## operator LP subproblem
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell,  basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = getSlaveSocpModelWithoutCascade(kcval, kgval, delta, deltaV0)

    println("Entering loop");
    jmax = N >= 36? 5 : N >= 24? 800 : N == 12 ? 100 : 100;
    lastChange = 0;
    j = 0; cumTime = 0; tic();
    prevObjv = 0;
    for j = 1:jmax
        # println(mm);
        nCurrentIter += 1;

        if sum(delta) > cardinality

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

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, deltav, deltaV0v = getBendersCascadeIterationForSocp(mm, deltaVar, deltaV0Var, verbose, LLCreq, smmip, kcmip, kgmip, deltaineqmip, sm, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta, ell, dgAncestorSuccessorCuts, withCuts, epsilon, withVariableEpsilon, nhighest);

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
            resilience = floor(100(1-prevObjv/LLCmax),2);
            resiliences = [resiliences; resilience]
            println(resilience, ", (", j+1, "), ", round(cumTime,2), ", ", round(100cardinality/N,2));
            lastChange = nCurrentIter;
            break;
        end

        prevObjv = objv;
    end
    println("Exited loop");

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = evaluateSlaveSocpModelWithCascade(delta, deltaV0);

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
        for ij = 2:length(nIters)
            println(resiliences[ij], ", (",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], " \\\\");
        end
    end

    global nBendersCascadeIter = j;
    llc = objv;
    delta, deltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getLinearAssistBendersCascade(LLCreq)
    deltaV0 = 0; delta = zeros(Int64,N); delta[UDG] = 1; kL = 0; kR = 100; ludg = length(UDG); minkpercent = 100;

    pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
    lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
    if objv < LLCreq
        println("LLCreq not achievable");
        minkpercent = 100;
    else
        for i = 1:length(UDG)
            delta[UDG[i]] = 0;
            pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
            lacv, lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
            minkpercent = 100sum(delta)/length(UDG) - 0.0001;
            println("objv, llcreq: ",objv," ",LLCreq);
            if objv < LLCreq
                break;
            end
        end
    end
    println("minkpercent: ", minkpercent);
    delta, deltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv = getBendersMethodCascade(LLCreq, minkpercent);

    delta, deltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function getBendersCascadeMinNodesLLC(epsilon, withVariableEpsilon, nhighest)
    kcval = ov; kgval = ov; delta = ov; deltaV0 = 0;
    withCascade = true; verbose = false; curLLC = 0;
    minNodesv = zeros(Int64,1); sumkcvs = zeros(Int64,1); sumkgvs = zeros(Int64,1);
    minDeltav = zeros(Int64,N);
    objvs = []; tempLLCreq = 0;
    printSummary = true; nIters = [0]; nSecs = [0]; nCurrentIter = 0; nCardinalities = [0]; lastChange = 0; cardinality = 0; objv = 0; withCuts = true; resiliences = [100];

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    # mm = addDownstreamHeuristicCuts(mm, deltaVar);

    ## operator MIP subproblem
    smmip, pmip, qmip, pgmip, qgmip, betamip, Pmip, Qmip, vmip, tmip, v0mip, kcmip, kgmip, xvarmip, vparmip, ineqmip, ineqbmip, eqmip, eqbmip, deltaeqmip, deltaeqbmip, subvoltmip, deltaineqmip, deltaineqbmip, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts  =  getSlaveModelWithCascade(delta, deltaV0, withCuts);

    ## operator LP subproblem
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = getSlaveModelWithoutCascade(kcval, kgval, delta, deltaV0)

    println("Entering loop");
    print("j = ");
    #ktemp
    jmax = N > 12? 5000 : 1000; cumTime = 0; tic();
    for j = 1:jmax
        j % 20 == 0 ? print(j, " ") : nothing;
        nCurrentIter += 1;
        if sum(delta) > cardinality

            currTime = toq(); cumTime += currTime;
            nSecs = [nSecs; cumTime];
            nIters = [nIters; j];
            cardinality = sum(delta);
            nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience];
            println("resilience: ", resilience);
            println("j = ",j);
            lastChange = nCurrentIter;
            tic();
        end

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, deltav, deltaV0v = getBendersCascadeIteration(mm, deltaVar, deltaV0Var, verbose, tempLLCreq, smmip, kcmip, kgmip, deltaineqmip, sm, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta,  dgAncestorSuccessorCuts, withCuts, epsilon, withVariableEpsilon, nhighest);

        tempLLCreq = max(objv,tempLLCreq);

        if !isMasterModelNotSolvable
            delta = deltav; deltaV0 = deltaV0v;
        end
        if !isBendersCutAdded || isMasterModelNotSolvable
            println("Exit reasons : isMasterModelNotSolvable, isBendersCutAdded, hasLLCExceeded ", isMasterModelNotSolvable, " ", isBendersCutAdded, " ", hasLLCExceeded);
            currTime = toq(); cumTime += currTime; nSecs = [nSecs; cumTime]; nIters = [nIters; j]; cardinality = sum(delta); nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience];
            println("resilience: ",resilience);
            lastChange = nCurrentIter;
            break;
        end

        objvs = vcat(objvs, objv);

        if j > 1
            minDeltav = hcat(minDeltav, delta);
            # sumkcvs = vcat(sumkcvs, sum(kcv));
            # sumkgvs = vcat(sumkgvs, sum(kgv));
            nNodesAttacked = sum(delta);
            minNodesv = vcat(minNodesv,nNodesAttacked);
        end


    end

    if printPerformance

        nSecs = round.(nSecs,2);
        # println(nCardinalities)
        nNodesAttackedPercentage = round.(100nCardinalities/N,1);
        # println(nCardinalities)
        for ij = 1:length(nIters)
            println(resiliences[ij], ", (",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], ", ", nNodesAttackedPercentage[ij], " \\\\");
        end
    end

    minNodesv, minDeltav, objvs
end

function getBendersSocpCascadeMinNodesLLC(epsilon, withVariableEpsilon, nhighest)
    kcval = ov; kgval = ov; delta = ov; deltaV0 = 0;
    withCascade = true; verbose = false; curLLC = 0;
    minNodesv = zeros(Int64,1); sumkcvs = zeros(Int64,1); sumkgvs = zeros(Int64,1);
    minDeltav = zeros(Int64,N);
    objvs = []; tempLLCreq = 0;
    printSummary = true; nIters = [0]; nSecs = [0]; nCurrentIter = 0; nCardinalities = [0]; lastChange = 0; cardinality = 0; objv = 0; withCuts = true; resiliences = [100];

    mm, deltaVar = getMasterModel();
    @variable(mm, deltaV0min <= deltaV0Var <= 1, Int);
    # mm = addDownstreamHeuristicCuts(mm, deltaVar);

    ## operator MIP subproblem
    smmip, pmip, qmip, pgmip, qgmip, betamip, Pmip, Qmip, vmip, tmip, v0mip, kcmip, kgmip, xvarmip, vparmip, ineqmip, ineqbmip, eqmip, eqbmip, deltaeqmip, deltaeqbmip, subvoltmip, deltaineqmip, deltaineqbmip, loadAncestorSuccessorCuts, dgAncestorSuccessorCuts, ellmip  = getSlaveSocpModelWithCascade(delta, deltaV0, withCuts);

    ## operator LP subproblem
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell,  basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt = getSlaveSocpModelWithoutCascade(kcval, kgval, delta, deltaV0)

    println("Entering loop");
    print("j = ");
    #ktemp
    jmax = N > 12? 5000 : 1000; cumTime = 0; tic();
    for j = 1:jmax
        j % 20 == 0 ? print(j, " ") : nothing;
        nCurrentIter += 1;
        if sum(delta) > cardinality

            currTime = toq(); cumTime += currTime;
            nSecs = [nSecs; cumTime];
            nIters = [nIters; j];
            cardinality = sum(delta);
            nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience];
            println("resilience: ", resilience);
            println("j = ",j);
            lastChange = nCurrentIter;
            tic();
        end

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, objv, deltav, deltaV0v = getBendersCascadeIterationForSocp(mm, deltaVar, deltaV0Var, verbose, tempLLCreq, smmip, kcmip, kgmip, deltaineqmip, sm, basicIneq, basicIneqb, deltaineq, deltaineqb, disineq, disb, basicEq, basicEqb, deltaeq, deltaeqb, ineq, ineqb, eq, eqb, subvolt, t, beta, ell, dgAncestorSuccessorCuts, withCuts, epsilon, withVariableEpsilon, nhighest);

        tempLLCreq = max(objv,tempLLCreq);

        if !isMasterModelNotSolvable
            delta = deltav; deltaV0 = deltaV0v;
        end
        if !isBendersCutAdded || isMasterModelNotSolvable
            println("Exit reasons : isMasterModelNotSolvable, isBendersCutAdded, hasLLCExceeded ", isMasterModelNotSolvable, " ", isBendersCutAdded, " ", hasLLCExceeded);
            currTime = toq(); cumTime += currTime; nSecs = [nSecs; cumTime]; nIters = [nIters; j]; cardinality = sum(delta); nCardinalities = [nCardinalities; cardinality];
            resilience = floor(100(1-objv/LLCmax),2);
            resiliences = [resiliences; resilience];
            println("resilience: ",resilience);
            lastChange = nCurrentIter;
            break;
        end

        objvs = vcat(objvs, objv);

        if j > 1
            minDeltav = hcat(minDeltav, delta);
            # sumkcvs = vcat(sumkcvs, sum(kcv));
            # sumkgvs = vcat(sumkgvs, sum(kgv));
            nNodesAttacked = sum(delta);
            minNodesv = vcat(minNodesv,nNodesAttacked);
        end


    end

    if printPerformance

        nSecs = round.(nSecs,2);
        # println(nCardinalities)
        nNodesAttackedPercentage = round.(100nCardinalities/N,1);
        # println(nCardinalities)
        for ij = 1:length(nIters)
            println(resiliences[ij], ", (",nIters[ij],"), ", nSecs[ij],", ", nCardinalities[ij], ", ", nNodesAttackedPercentage[ij], " \\\\");
        end
    end

    minNodesv, minDeltav, objvs
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

    global ancestorLoads, succesorLoads, ancestorDGs, succesorDGs

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

    # println("DG ", DG);

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
    # println("RES ", RES);

    # EG = copy(MGL);
    EG = RES;
    # println("MGL ", MGL);
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

    ancestorLoads, succesorLoads = getAncestorSuccessorLoadPairs();
    ancestorDGs, succesorDGs = getAncestorSuccessorDGPairs();

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

function addVariablesToBasicSlaveRecoveryModel(sm)
    @variable(sm, p[1:N,1:T]); @variable(sm, q[1:N,1:T]);
    @variable(sm, pgvar[1:length(DG),1:T]); @variable(sm, qgvar[1:length(DG),1:T]);

    # @variable(sm, pgvar[1:length(RES)]); @variable(sm, qgvar[1:length(RES)]);
    @variable(sm, betavar[1:length(LOADS),1:T]);

    #voltages are squared
    @variable(sm, v[1:N,1:T]); @variable(sm, v0[1:T]);
    # @variable(sm, t[1:N]);
    @variable(sm, tvar[1:T]);
    @variable(sm, P[1:N,1:T]); @variable(sm, Q[1:N,1:T]);


    pg = zeros(AffExpr,N,T); qg = zeros(AffExpr,N,T); beta = zeros(AffExpr,N,T); t = zeros(AffExpr,N,T);

    pg[DG,:] = pgvar; qg[DG,:] = qgvar; beta[LOADS,:] = betavar;
    vpar = zeros(AffExpr,N,T);
    # pg[RES] = pgvar; qg[RES] = qgvar;

    for ts = 1:T
        t[:,ts] = tvar[ts];
    end
    for i = 1:N
        vpar[i,:] = par[i] == 0 ? v0' : v[par[i],:];
    end

    p, q, pg, qg, beta, P, Q, v, t, v0, vpar
end

function getBasicSlaveRecoverySocpModel(withBinary=true)
    sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));
    if withBinary == false
        sm = Model(solver=MosekSolver(QUIET=true));
    end
    # sm = Model(solver=GurobiSolver(OutputFlag = 0, genv));
    p, q, pg, qg, beta, P, Q, v, t, v0, vpar = addVariablesToBasicSlaveRecoveryModel(sm);

    @variable(sm, ell[1:N,1:T]);
    xvar = vcat(xvar, ell);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell
end

function getRecoveryCommonConstraintsObjectiveSocp(sm, delta, deltaV0, kcval, kgval)

    ### attack inequalities
    @constraint(sm, attackInitial[i in DG, ts = 1], kgval[i,ts] >=  delta[i]);
    @constraint(sm, dgactiveAttack[i in DG, ts = 1], pg[i,ts] <= (1-delta[i]) * pgmax[i]);
    @constraint(sm, dgreactiveAttack[i in DG, ts = 1], qg[i,ts] <= (1-delta[i]) * qgmax[i]);

    ### DG model
    @constraint(sm, dgactiveOperator[i in DG, ts in TS], pg[i,ts] <= (1-kgval[i,ts]) * pgmax[i]);
    @constraint(sm, dgreactiveOperator[i in DG, ts in TS], qg[i,ts] <= (1-kgval[i,ts]) * qgmax[i]);
    @constraint(sm, dgactiveLB[i in DG, ts in TS], pg[i,ts] >= 0);
    @constraint(sm, dgreactiveLB[i in DG, ts in TS], qg[i,ts] >= -(1-kgval[i,ts]) * qgmax[i]);
    @constraint(sm, dgDisconnectLB[i in DG, ts in TS], v[i,ts] >= vgmin[i] - kgval[i,ts]);
    @constraint(sm, dgDisconnectUB[i in DG, ts in TS], v[i,ts] <= vgmax[i] + kgval[i,ts]);

    ### Load model
    # load control
    @constraint(sm, loadLow[i in LOADS, ts in TS], beta[i,ts] >= (1 - kcval[i,ts]) * betamin[i]);
    @constraint(sm, loadUp[i in LOADS, ts in TS], beta[i,ts] <= 1 - kcval[i,ts]);
    @constraint(sm, loadDisconnectLB[i in LOADS, ts in TS], v[i,ts] >= vcmin[i] - kcval[i,ts]);
    @constraint(sm, loadDisconnectUB[i in LOADS, ts in TS], v[i,ts] <= vcmax[i] + kcval[i,ts]);

    ### Power flow constraints
    ## power conservation
    @constraint(sm, realCons[i in nodes, ts in TS], p[i,ts] - beta[i,ts] * pcmax[i] + pg[i,ts] == 0);
    @constraint(sm, reacCons[i in nodes, ts in TS], q[i,ts] - beta[i,ts] * qcmax[i] + qg[i,ts] == 0);

    @constraint(sm, realFlow[i in nodes, ts in TS], P[i,ts] == sum(CHILD[i,j] * P[j,ts] for j in nodes) + p[i,ts] + resistance[i] * ell[i]);
    @constraint(sm, reacFlow[i in nodes, ts in TS], Q[i,ts] == sum(CHILD[i,j] * Q[j,ts] for j in nodes) + q[i,ts] + reactance[i] * ell[i]);

    ## Voltage constraints
    if length(LTC) > 0
        @constraint(sm, voltDropLTC[i in LTC, ts in TS], v[i,ts] - LTCSetting * vpar[i,ts] .== 0);
    end
    @constraint(sm, voltDropNoLTC[i in noLTC, ts in TS], vpar[i,ts] - v[i,ts] - 2resistance[i] * P[i,ts] - 2reactance[i] * Q[i,ts] + (resistance[i]^2 + reactance[i]^2) * ell[i,ts] == 0)

    ## Current constraints
    @constraint(sm, socpIneq[i=1:N], norm([vpar[i], ell[i], sqrt(2)*P[i], sqrt(2)*Q[i]]) <= (vpar[i] + ell[i]));

    ### Initial conditions NOT required
    # @constraint(sm, initialAttack[i in DG, ts = 1], kg[i,ts] == delta[i]);

    ### Recovery conditions
    @constraint(sm, recoveryMonotonicity[i in DG, ts = 1:T-1], kgval[i,ts+1] <= kgval[i,ts]); # DGs once connected remain connected

    @constraint(sm, recoveryBudget[ts = 1:T-1], sum(kgval[i,ts] - kgval[i,ts+1] for i in DG) <= NG);
    @constraint(sm, subvoltAttack[ts = 1:T-1], v0[ts] == v0nom + v0dis * deltaV0);

    ### Final conditions
    @constraint(sm, subvoltFinal, v0[T] == v0nom);

    ### Objective
    @constraint(sm, voltageLB[i in nodes, ts in TS], t[i,ts] + v[i,ts] >= v0nom);
    @constraint(sm, voltageUB[i in nodes, ts in TS], t[i,ts] - v[i,ts] >= -v0nom);

    @variable(sm, sysPerf[TS]); @variable(sm, sysCost[TS]);

    @constraint(sm, performanceLoss[ts in TS], sysCost[ts] == WVR * sum(t[:,ts]) / N +
    sum(WLC[i]*(1-beta)[i,ts] for i in LOADS) +
    sum((WSD-WLC)[i] * kcval[i,ts] for i in LOADS) +
    sum(WAC)*sum(resistance[i] * ell[i,ts] for i in nodes));

    @constraint(sm, resilience[ts in TS], sysPerf[ts] == 100 * (1 - sysCost[ts]/LLCmax));

    @objective(sm, :Min, sum(sysCost[ts] for ts in TS));

    sysCost, sysPerf, attackInitial, dgactiveAttack, dgreactiveAttack, dgactiveOperator, dgreactiveOperator, dgactiveLB, dgreactiveLB, dgDisconnectLB, dgDisconnectUB, loadLow, loadUp, loadDisconnectLB, loadDisconnectUB, subvoltAttack, subvoltFinal, voltageLB, voltageUB, performanceLoss, resilience
end

function getRecoverySocpModelWithBinary(delta, deltaV0)
    withBinary=true;
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell = getBasicSlaveRecoverySocpModel(withBinary);

    @variable(sm, kcvar[1:length(LOADS),1:T], Bin);
    @variable(sm, kgvar[1:length(DG),1:T], Bin);
    kc = zeros(AffExpr,N,T); kg = zeros(AffExpr,N,T);
    kc[LOADS,:] = kcvar[:,:]; kg[DG,:] = kgvar[:,:];

    sysCost, sysPerf, attackInitial, dgactiveAttack, dgreactiveAttack, dgactiveOperator, dgreactiveOperator, dgactiveLB, dgreactiveLB, dgDisconnectLB, dgDisconnectUB, loadLow, loadUp, loadDisconnectLB, loadDisconnectUB, subvoltAttack, subvoltFinal, voltageLB, voltageUB, resilience = getRecoveryCommonConstraintsObjectiveSocp(sm, delta, deltaV0, kc, kg);

    println(sm);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell, sysCost, sysPerf
end

function getRecoverySocpModelWithoutBinary(delta, deltaV0, kcv, kgv)
    withBinary=false;
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell = getBasicSlaveRecoverySocpModel(withBinary);

    sysCost, sysPerf, attackInitial, dgactiveAttack, dgreactiveAttack, dgactiveOperator, dgreactiveOperator, dgactiveLB, dgreactiveLB, dgDisconnectLB, dgDisconnectUB, loadLow, loadUp, loadDisconnectLB, loadDisconnectUB, subvoltAttack, subvoltFinal, voltageLB, voltageUB, performanceLoss, resilience = getRecoveryCommonConstraintsObjectiveSocp(sm, delta, deltaV0, kcv, kgv);

    println(sm);

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell, sysCost, sysPerf, attackInitial, dgactiveAttack, dgreactiveAttack, dgactiveOperator, dgreactiveOperator, dgactiveLB, dgreactiveLB, dgDisconnectLB, dgDisconnectUB, loadLow, loadUp, loadDisconnectLB, loadDisconnectUB, subvoltAttack, subvoltFinal, voltageLB, voltageUB, resilience
end

function evaluateSlaveRecoverySocpModelWithBinary(delta, deltaV0)
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell, sysCost, sysPerf = getRecoverySocpModelWithBinary(delta, deltaV0);

    # println(sm);
    sstatus = solve(sm);

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, Pv, Qv, vv, v0v, tv, ellv, objv, sysCostv, sysPerfv = getRecoveryVariableValues(sm, sstatus, beta, kc, kg, pg, qg, p, q, P, Q, v, v0, t, sysCost, sysPerf)
    # printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, Pv, Qv, vv, v0v, tv, ellv, objv, sysCostv, sysPerfv
end

function evaluateSlaveRecoverySocpModelWithBinaryUsingGreedy(delta, deltaV0)
    withCascade=true;
    sm, p, q, pg, qg, beta, P, Q, v, t, v0, xvar, vpar, ell = getBasicSlaveSocpModel(withCascade);

    @variable(sm, kcvar[1:length(LOADS)], Bin); # kc
    @variable(sm, kgvar[1:length(DG)], Bin); # kg
    kc = zeros(AffExpr,N); kg = zeros(AffExpr,N); kc[LOADS] = kcvar; kg[DG] = kgvar;
    addObjectiveDisconnectSocpModel(sm, t, beta, kc, ell);

    # addAttackInequalities
    kgv = delta;
    @constraint(sm, dgAttack1[i in DG], kg[i] >= 0);
    @constraint(sm, dgAttack2, pg[DG] .== (1 - kg[DG]) .* pgmax[DG]);
    @constraint(sm, dgAttack3, qg[DG] .== (1 - kg[DG]) .* qgmax[DG]);

    # addAttackEqualities
    @constraint(sm, subvolt, v0 == v0nom + v0dis * deltaV0);

    # Basic inequality and equality constraints of SOCP model
    sm,  basicIneq, basicIneqb = addBasicSocpInequalities(sm, pg, qg, beta, v, t, kcv, kgv, vpar, P, Q, ell);
    sm, basicSocpEq, basicSocpEqb = addBasicSocpEqualities(sm, v, vpar, p, q, pg, qg, beta, P, Q, kgv, ell);
    sm, disineq, disb = addDisconnectivityConstraints(sm, kcv, kgv, v);

    @constraint(sm, recoveryMonotonicity[i in DG], kg[i] <=  kgv[i]); # DGs once connected remain connected

    @constraint(sm, recoveryBudget, sum(kgv[i] - kg[i] for i in DG) <= NG);

    sysCostv = zeros(Float64, T);
    for ts = 1:T
        if ts == 1
            for i in DG
                JuMP.setRHS(recoveryMonotonicity[i], 1)
            end
            for i in DG
                JuMP.setRHS(dgAttack1[i], delta[i])
            end
        elseif ts == T
            JuMP.setRHS(subvolt, v0nom)
        elseif sum(kgv) == 0
            for i in DG
                JuMP.setRHS(recoveryMonotonicity[i], kgv[i])
            end
        else # 1 < ts < T and all DGs not restored
            for i in DG
                JuMP.setRHS(recoveryMonotonicity[i], kgv[i])
            end
            nReparableDGs = min(NG, sum(kgv))
            @constraint(sm, recoveryBudget, sum(kgv[i] - kg[i] for i in DG) == nReparableDGs);
        end
    end

    # println(sm);
    sstatus = solve(sm);

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, Pv, Qv, vv, v0v, tv, ellv, objv, sysCostv, sysPerfv = getRecoveryVariableValues(sm, sstatus, beta, kc, kg, pg, qg, p, q, P, Q, v, v0, t, sysCost, sysPerf)
    # printResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv)

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, Pv, Qv, vv, v0v, tv, ellv, objv, sysCostv, sysPerfv
end

function getRecoveryVariableValues(sm, sstatus, beta, kc, kg, pg, qg, p, q, P, Q, v, v0, t, sysCost, sysPerf)
    assert(sstatus == :Optimal);
    pcv = zeros(Float64, length(LOADS), T); qcv = zeros(Float64, length(LOADS), T);

    betav = getvalue(beta); kcv = getvalue(kc);
    # println(JuMP.size(betav), " ", size(pcmax));
    for i in LOADS, ts in TS
        pcv[i,ts] = betav[i,ts] * pcmax[i]; qcv[i,ts] = betav[i,ts] * qcmax[i];
    end

    kgv = getvalue(kg); pgv = getvalue(pg); qgv = getvalue(qg);
    pv = getvalue(p); qv = getvalue(q);

    Pv = getvalue(P); Qv = getvalue(Q);
    vv = getvalue(v); v0v = getvalue(v0);
    tv = getvalue(t); ellv = getvalue(ell);

    sysCostv = getvalue(sysCost); sysPerfv = getvalue(sysPerf);

    objv = getobjectivevalue(sm);

    # betav, kcv, pcv, qcv, kgv, pgv, qgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tv, objv, lacv, lovrv, lcontv, lsdv = round.([betav, kcv, pcv, qcv, kgv, pgv, qgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tv, objv, lacv, lovrv, lcontv, lsdv],3);

    # kcv = round.(Int64, kcv); kgv = round.(Int64,kgv); kmv = round.(Int64, kmv);

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, Pv, Qv, vv, v0v, tv, ellv, objv, sysCostv, sysPerfv
end

function getBendersRecoveryIterationForSocp(mm, deltaVar, deltaV0Var, verbose, LLCreq)
    println("In getBendersCascadeIterationForSocp");
    hasLLCExceeded = false; isBendersCutAdded = false; isMasterModelNotSolvable = false;

    pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0, ellv = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    # println(mm);
    if verbose
        println(mm);
    end

    mstatus = solve(mm);

    if mstatus != :Optimal
        println("Master model solve status : ", mstatus);
        isMasterModelNotSolvable = true;
        return hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, kmv, lacv, lovrv, lcontv, lsdv, delta, deltaV0;
    end

    delta = round.(Int64,getvalue(deltaVar)); deltaV0 = getvalue(deltaV0Var);
    nNodesAttacked = sum(round.(Int64,delta));
    # println("delta : ", delta);
    if verbose
        println("no. of nodes attacked : ", nNodesAttacked);
        println("delta : ", nodes[delta .== 1]);
    end

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell, sysCost, sysPerf = getRecoverySocpModelWithBinary(delta, deltaV0);

    if verbose
        println("Updated slave model");
        # println(sm);
    end

    sstatus = solve(sm);
    if sstatus != :Optimal
        println("Slave model solve status : ", sstatus);
        return delta, deltaV0, pv, pgv, qv, qgv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, lacv, lovrv, lcontv, lsdv, ellv

    end

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, Pv, Qv, vv, v0v, tv, ellv, objvMIP, sysCostv, sysPerfv = getRecoveryVariableValues(sm, sstatus, beta, kc, kg, pg, qg, p, q, P, Q, v, v0, t, sysCost, sysPerf);

    verbose ? printSocpResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objvMIP, pgv, qgv, ellv) : nothing;

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

    sm, p, q, pg, qg, beta, P, Q, v, t, v0, vpar, ell, sysCost, sysPerf, attackInitial, dgactiveAttack, dgreactiveAttack, dgactiveOperator, dgreactiveOperator, dgactiveLB, dgreactiveLB, dgDisconnectLB, dgDisconnectUB, loadLow, loadUp, loadDisconnectLB, loadDisconnectUB, subvoltAttack, subvoltFinal, voltageLB, voltageUB, performanceLoss, resilience = getRecoverySocpModelWithoutBinary(delta, deltaV0, kcv, kgv);

    if verbose
        println("****************************************");
        println("Printing slave model without cascade");
        println(sm);
    end

    sstatus = solve(sm);
    assert(sstatus == :Optimal)
    objv = getobjectivevalue(sm);

    verbose ? println("objvMIP = ", objvMIP, ", objv = ",objv) : nothing;
    verbose ? printSocpResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, ellv) : nothing;

    lambdaTimesb =
    sum(getdual(dgactiveAttack)[i,ts] * pgmax[i] for i in DG, ts = 1) +
    sum(getdual(dgreactiveAttack)[i,ts] * qgmax[i] for i in DG, ts = 1) +
    sum(getdual(dgactiveOperator)[i,ts] * (1-kgv[i,ts]) * pgmax[i] for i in DG, ts in TS) +
    sum(getdual(dgreactiveOperator)[i,ts] * (1-kgv[i,ts]) * qgmax[i] for i in DG, ts in TS) +
    sum(getdual(dgreactiveLB)[i,ts] * -(1-kgv[i,ts]) * qgmax[i] for i in DG, ts in TS) +
    sum(getdual(dgDisconnectLB)[i,ts] * (vgmin[i] - kgv[i,ts]) for i in DG, ts in TS) +
    sum(getdual(dgDisconnectUB)[i,ts] * (vgmax[i] + kgv[i,ts])  for i in DG, ts in TS) +
    sum(getdual(loadLow)[i,ts] * (1 - kcv[i,ts]) * betamin[i] for i in LOADS, ts in TS) +
    sum(getdual(loadUp)[i,ts] * (1 - kcv[i,ts])  for i in LOADS, ts in TS) +
    sum(getdual(loadDisconnectLB)[i,ts] * (vcmin[i] - kcv[i,ts]) for i in LOADS, ts in TS) +1
    sum(getdual(loadDisconnectUB)[i,ts] * (vcmax[i] + kcv[i,ts]) for i in LOADS, ts in TS) +
    sum(getdual(subvoltAttack)[ts] * v0nom for ts = 1:T-1) +
    getdual(subvoltFinal) * v0nom +
    sum(getdual(voltageLB)[i,ts] * v0nom for i in nodes, ts in TS) +
    sum(getdual(voltageUB)[i,ts] * -v0nom for i in nodes, ts in TS) +
    sum(getdual(performanceLoss)[ts] * (sum(WLC[i] for i in LOADS) + sum((WSD[i]-WLC[i]) * kcv[i,ts] for i in LOADS)) for ts in TS) +
    sum(getdual(resilience)[ts] * 100 for ts in TS);

    lambdaTimesBTimesd =
    sum(getdual(dgactiveAttack)[i,ts] * pgmax[i] * deltaVar[i] for i in DG, ts = 1) +
    sum(getdual(dgreactiveAttack)[i,ts] * qgmax[i] * deltaVar[i] for i in DG, ts = 1) +
    sum(getdual(subvoltAttack)[ts] * v0dis * deltaV0Var for ts = 1:T-1);

    # if maximum(abs.(deltaeqnubar)) + maximum(abs.(deltaIneqlbar)) >= myinf
    bendersCutLHS = sum((WSD-WLC).*kcv) + sum(WLC[LOADS]) + lambdaTimesb + lambdaTimesBTimesd;# + sum(socpIneqtbar.*socpIneqb);
    if !trial
        #epsilon = 20; # for N = 12 or 24
        epsilon = 100
        bendersCutLHS -= (objv + epsilon);
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

    isBendersCutAdded = true;
    # else
    #     printSocpResults(pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, ellv);
    # end

    hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, delta, deltaV0, ellv
end

function getBendersMethodRecoverySocp(LLCreq)
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

        hasLLCExceeded, isBendersCutAdded, isMasterModelNotSolvable, pv, qv, betav, vv, v0v, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, Pv, Qv, lovrv, lcontv, lsdv, deltav, deltaV0v, ellv = getBendersCascadeIterationForSocp(mm, deltaVar, deltaV0Var, verbose, LLCreq);

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

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pv, qv, Pv, Qv, vv, v0v, tv, ellv, objv, sysCostv, sysPerfv = evaluateSlaveRecoverySocpModelWithBinary(delta, deltaV0);

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

    delta, deltaV0, pv, qv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv, ellv
end

function evaluateRecoveryModel(delta, deltaV0)
    sm, beta, kc, kg, pg, qg, p, q, km, P, Q, v, v0, t, Lac, Lovr, Lsd, Lcont, Lmg, sysCost, sysPerf = getRecoveryModel(delta, deltaV0);
    sstatus = solve(sm);

    betav, kcv, pcv, qcv, kgv, pgv, qgv, pgv, qgv, pv, qv, kmv, Pv, Qv, vv, v0v, tv, objv, lacv, lovrv, lcontv, lsdv, sysCostv, sysPerfv = getRecoveryVariableValues(sm, sstatus, beta, kc, kg, pg, qg, p, q, km, P, Q, v, v0, t, Lac, Lovr, Lsd, Lcont, Lmg, sysCost, sysPerf)
end

function computeOptimalityGap(bfMs, bfobjvs, bnMs, bnobjvs)
    fiter, biter = 1,1
    ludg = length(UDG);
    maxOptimalityGap = 0;
    avgOptimalityGap = 0;
    sumOptimalityGap = 0;
    for fiter = 1:length(bfMs)
        target = bfobjvs[fiter];
        while (biter <= length(bnMs) && bnobjvs[biter] < target)
            biter += 1
        end
        if bnobjvs[biter] >= target - 0.01 # avoiding numerical issues
            optimalityGap = 100*(bnMs[biter] - bfMs[fiter])/ludg;
        end
        maxOptimalityGap = max(maxOptimalityGap, optimalityGap)
        sumOptimalityGap += optimalityGap;
    end
    avgOptimalityGap = sumOptimalityGap/length(bfMs);
    return maxOptimalityGap, avgOptimalityGap
end
