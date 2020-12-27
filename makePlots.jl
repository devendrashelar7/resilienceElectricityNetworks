begin
  clf();
  fig = figure("pyplot_subplot_mixed",figsize=(10,10)) # Create a new blank figure
  x1 = [0,0.5,1,1,0.5,0,0]; y1 = [1,1,1/3,-1/3,-1,-1,1];
  plot(x1,y1, linewidth=4,linestyle="-",color="black");
  xlim(0,1.2);
  ylim(-1.2,1.2)
  savefig(string(mypath,"myderpolytope.pdf"));
end



begin
  N = 6; setGlobalParameters();
  delta = zeros(Int64,N); deltaV0 = 0;
  srmaxp = srmax; pc2sr = 0; pg2sr = 1; sd2sr = 0; se2sr = 1;
  modifyLoadAndDGParameters(srmaxp, pc2sr, pg2sr, sd2sr, se2sr);
  println(pcmax)

  clf();
  fig = figure("voltageInScenarios",figsize=(10,10)) # Create a new blank figure
  pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv,
  lovrv, lcontv, lsdv = evaluateSlaveModelWithCascade(delta, deltaV0);
  plot(linspace(0,8,9),[1;vv[1:8]], linewidth=4,linestyle="-",color="black");
  # xlim(0,1.2);
  # ylim(-1.2,1.2)
  println(pcv)
  savefig(string(mypath,"voltage12.pdf"));
end

# voltages for different scenarios - nominal - TN - DN disturbance - cascade
begin
    N = 15; fs = 25; font1 = Dict("family"=>"serif","color"=>"black","weight"=>"normal","size"=>fs);
    # p = 0.1ones(Float64,N);
    pcp = 0.09; dp = 0.03;
    r = 0.01ones(Float64,N);
    m = Model(solver=GurobiSolver()); @objective(m, :Min, 0);
    @variable(m, p[1:N]); @variable(m, v0)
    @variable(m, v[1:N]); @variable(m, P[1:N]); allV = [v0; v];
    vpar = ones(AffExpr,N); vpar[1] = v0; vpar[2:N] = v[1:N-1];
    LTC = [floor(Int64,N/3)]; noLTC = setdiff(collect(1:N), LTC);
    # for i = 1:N
    #     @constraint(m, P[i] == sum(p[i:N]));
    # end
    @constraint(m, cons[i=1:N], P[i] == sum(p[i:N]));
    @constraint(m, real[i=1:N], p[i] == pcp);
    @constraint(m, volt[i=1:length(noLTC)], (vpar - v)[noLTC[i]] == (r .* P)[noLTC[i]]);
    @constraint(m, (1.05*vpar-v)[LTC] .== 0);
    @constraint(m, subVolt, v0 == 1.04);

    solve(m); vv = getvalue(allV);
    # println(round(vv,2));

    clf();
    fig = figure("voltageInScenarios",figsize=(5,3));
    # Create a new blank figure
    new_position = [0.15,0.10,0.8,0.8] # Position Method 2
    # ylim(0,50)
    ax = gca();
    ax[:set_position](new_position)

    xx = collect(0:N); lw = 3;
    plot(xx, vv, linewidth = lw, label="nominal", linestyle="-",color="black");

    JuMP.setRHS(subVolt, 1.02);
    solve(m); vv = getvalue(allV);

    plot(xx, vv, linewidth = lw, label="only TN-side disruption", linestyle="--",color="green");

    for i = 1:N
        JuMP.setRHS(real[i], 0.12);
    end
    solve(m); vv = getvalue(allV);

    plot(xx, vv, linewidth = lw, label="TN/DN disruption", linestyle="--",color="blue");

    for i = 10:N
        JuMP.setRHS(real[i], 0.15);
    end
    solve(m); vv = getvalue(allV);
    plot(xx, vv, linewidth = lw, label="cascade", linestyle="-",color="red");

    JuMP.setRHS(volt[9], -0.15);
    solve(m); vv = getvalue(allV);
    # plot(xx, vv, linewidth = 2lw, label="microgrid", linestyle=":",color="blue");


    plot(xx, 1.05ones(xx), color="black")
    plot(xx, 0.95ones(xx), color="black")

    xlabel("Distance from substation", fontdict=font1); ylabel("Voltage (in pu)", fontdict=font1);
    # ticks(xmajor = false);
    legend(loc = "upper right", fontsize=fs); ## legend position
    # setp(ax[:get_xticklabels](),visible="off");
    setp(ax[:get_yticklabels](),fontsize=fs); setp(ax[:get_xticklabels](),fontsize=fs);
    xlim(0,N);
    ylim(0.80,1.2);
    savefig(string(mypath,"voltageN.pdf"));


end

Pkg.add("Convex")
using Convex

begin
    m = Model(solver=GurobiSolver()); x = 0;
    @variable(m, y[1:2] >= 1);
    @variable(m, t );

    @constraint(m, soc, norm(y[i] for i = 1:2) <= t);
    # @QuadConstraint(m, sum(y.^2) <= 2);
    # @constraint(m, ssrhs, t == 1);

    @objective(m, :Min, t);
    solve(m)
    println(m)
    yv = getvalue(y);
    println("sol: ", yv);
    println("dual: ", getdual(soc));
end
