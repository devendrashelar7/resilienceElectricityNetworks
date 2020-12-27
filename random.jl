attackSets = collect(combinations([1,2,3],2))
begin
  res = @parallel vcat for i = 1:10
    j = i * 2
    k = j * 2 - 2
    j, k
  end
  println(res)
end
# begin
  # addprocs(Sys.CPU_CORES - 1);
#   @everywhere using DistributedArrays
#   y = @DArray [@show i + j for i = 1:3, j = 4:6];
#   println(y)
# end

Pkg.update()
using JuMP
using Gurobi
# Pkg.add("Ipopt")
using Ipopt
ENV["MOSEKBINDIR"] = "/Users/devendrashelar/mosek/9.1/tools/platform/osx64x86/bin"

# Pkg.rm("Mosek")
Pkg.add("Mosek")
Pkg.build("Mosek")
using Mosek
begin
  # m = Model(solver=IpoptSolver(print_level=0));
  m = Model(solver=MosekSolver());
  # m = Model(solver=GurobiSolver());
  @variable(m, x[1:5] >= 0);
  y = zeros(AffExpr,2);
  y[1] = x[2];
  # @NLconstraint(m, xcon[i=1:5], x[i]^2 + x[6-i]^2 +(-1)^i*x[i]*x[6-i] <= 5*i);
  @constraint(m, xcon[i=1:2], norm([x[i]+x[6-i],x[i],x[6-i]]) <= 5*i);
  # @NLconstraint(m, xcon, x[1]^2 + x[3-1]^2 +(-1)^1*x[1]*x[3-1] <= 5*1);
  @objective(m, :Min, sum(x[i] for i = 1:5));
  # @NLobjective(m, :Min, sum(x));

  println(m)
  # println(xcon[1])
  # println(xcon[2])

  solve(m);
  println(getvalue(x));
  for i = 1:1
      println(getdual(xcon[i]))
  end
  # println(getdual(xcon[2]))
end

begin
  # m = Model(solver=IpoptSolver(print_level=0));
  m = Model(solver=MosekSolver());
  # m = Model(solver=GurobiSolver());
  @variable(m, x);
  @variable(m, y);
  @variable(m, z);

  @constraint(m, ylb, y >= 4);
  @constraint(m, zlb, z >= 4);
  @constraint(m, xcon, norm([y,z]) <= x);
  @objective(m, :Min, x);
  # @NLobjective(m, :Min, sum(x));

  println(m)
  # println(xcon[1])
  # println(xcon[2])

  solve(m);
  println("x = ",getvalue(x));
  println("dual = ",getdual(xcon))
  println("dual ylb = ",getdual(ylb))
  println("dual zlb = ",getdual(zlb))

  # println(getdual(xcon[2]))
end

begin
  # m = Model(solver=IpoptSolver(print_level=0));
  m = Model(solver=MosekSolver());
  # m = Model(solver=GurobiSolver());
  @variable(m, w);
  @variable(m, x);
  @variable(m, y);
  @variable(m, z);

  @constraint(m, ylb, y >= 4);
  @constraint(m, zlb, z >= 4);
  # @constraint(m, xcon, norm([y,z]) <= w*x);
  # @constraint(m, xcon, y^2+z^2 <= w*x);
  @constraint(m, xcon, norm([w x sqrt(2)*y sqrt(2)*z]) <= w + x);
  @constraint(m, quadcon, norm([w x]) <= 10);
  # @constraint(m, quadcon, w^2+x^2 <= 100);
  @objective(m, :Min, w+x);
  # @NLobjective(m, :Min, sum(x));

  println(m)
  # println(xcon[1])
  # println(xcon[2])

  solve(m);
  println("obj = ",getobjectivevalue(m)," = ",8*sqrt(2));
  println("x = ",getvalue(x));
  println("w,x,y,z = ",getvalue([w,x,y,z]));
  println("dual = ",getdual(xcon))
  println("dual ylb = ",getdual(ylb))
  println("dual zlb = ",getdual(zlb))
  println("dual quadcon = ",getdual(quadcon))

  # println(getdual(xcon[2]))
end


begin
  m = Model(solver=GurobiSolver());
  @variable(m, x[1:2]);
  y = zeros(AffExpr,2);
  y[1] = x[2];
  @constraint(m, xcon, x + y.>= 1);
  @constraint(m, xcon, x + 2y.>= 1);
  @constraint(m, xcon, x + 3y.>= 1);
  @objective(m, :Min, sum(x));
  println(m)
  println(xcon[1])
  println(xcon[2])
end

begin
    println("************************************")
    master = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(master, x[1:2], Bin);
    @objective(master, :Min, sum(x));

    mstatus = solve(master);
    xv = zeros(size(x));
    if mstatus == :Optimal
        xv = getvalue(x)
    end

    slave = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(slave, y[1:2]);
    @objective(slave, :Min, -3y[1] - 4y[2]);
    @constraint(slave, xcon[i=1:2], y[i] <= 1 - xv[i]);

    j = 0
    while j < 5
        println("j = ",j)
        println(slave)
        sstatus = solve(slave);
        assert(sstatus == :Optimal)
        objv = getobjectivevalue(slave)
        if objv >= -0.001
            break;
        end

        xcond = getdual(xcon);
        @constraint(master, sum(xcond .* (1-x)) >= 0)
        println(master)
        mstatus = solve(master);
        assert(mstatus == :Optimal)
        xv = getvalue(x)
        println("value of x = ", xv)

        for i = 1:2
            JuMP.setRHS(xcon[i], 1-xv[i])
        end
        j+=1
    end


end

begin
    println("************************************************************************")
    master = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(master, x[1:2], Bin);
    @objective(master, :Min, sum(x));
    println("\nMaster program is \n",master)

    mstatus = solve(master);
    xv = zeros(size(x));
    if mstatus == :Optimal
        xv = getvalue(x)
    end

    bias = -10
    slave = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(slave, y[1:2]);
    @variable(slave, z[1:2], Bin);
    @objective(slave, :Min, -3y[1] - 4y[2] - z[1] + bias);
    @constraint(slave, xcon[i=1:2], z[i] <= 1 - xv[i]);
    @constraint(slave, ycon[i=1:2], y[i] <= z[i]);

    slaveLP = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(slaveLP, yLP[1:2]);
    zv = zeros(size(z));
    @objective(slaveLP, :Min, -3yLP[1] - 4yLP[2] - zv[1] + bias);
    @constraint(slaveLP, zcon[i=1:2], yLP[i] <= zv[i]);
    target = -10
    j = 0
    while j < 3
        # Solving Slave MIP
        println("j = ",j)
        println("slave problem is:\n", slave)
        sstatus = solve(slave);
        assert(sstatus == :Optimal)
        objv = getobjectivevalue(slave)
        println("slave MIP obj = ", objv)
        if objv >= target
            break;
        end

        # Solving Slave LP
        zv = getvalue(z);
        for i = 1:2
            JuMP.setRHS(zcon[i], zv[i])
        end
        println("slaveLP problem is:\n", slaveLP)
        sstatus = solve(slaveLP);
        assert(sstatus == :Optimal)

        zcond = getdual(zcon);
        @constraint(master, sum(zcond .* (1-x)) >= target - bias)
        println("\nMaster program is \n",master)
        mstatus = solve(master);
        assert(mstatus == :Optimal)
        xv = getvalue(x)
        println("value of x = ", xv)

        for i = 1:2
            JuMP.setRHS(xcon[i], 1-xv[i])
        end
        j+=1
    end


end

begin
    println("************************************************************************")
    master = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(master, x[1:2], Bin);
    @objective(master, :Min, sum(x));
    println("\nMaster program is \n",master)

    mstatus = solve(master);
    xv = zeros(size(x));
    if mstatus == :Optimal
        xv = getvalue(x)
    end

    bias = -10
    slave = Model(solver=MosekSolver(QUIET=true));
    @variable(slave, w[1:2]);
    @variable(slave, y[1:2]);
    @variable(slave, z[1:2], Bin);
    @objective(slave, :Min, sum(w) - z[1] + bias);
    @constraint(slave, xcon[i=1:2], z[i] <= 1 - xv[i]);
    @constraint(slave, ycon[i=1:2], y[i] <= z[i]);
    @constraint(slave, wcon, norm(y) <= sum(w));


    slaveLP = Model(solver=MosekSolver(QUIET=true));
    # slaveLP = Model(solver=GurobiSolver(OutputFlag=0));
    @variable(slaveLP, yLP[1:2]);
    @variable(slaveLP, wLP[1:2]);
    zv = zeros(size(z));
    @objective(slaveLP, :Min, sum(wLP) - zv[1] + bias);
    @constraint(slaveLP, zconLP[i=1:2], yLP[i] <= zv[i]);
    @constraint(slaveLP, wconLP, norm(yLP) <= sum(wLP));

    target = -10
    j = 0

    function addBendersCut(cb)
        xv = getvalue(x)
        println("value of x = ", xv)

        for i = 1:2
            JuMP.setRHS(xcon[i], 1-xv[i])
        end

        println("slave problem is:\n", slave)
        sstatus = solve(slave);
        assert(sstatus == :Optimal)
        objv = getobjectivevalue(slave)
        println("slave MIP obj = ", objv)

        # Solving Slave LP
        zv = getvalue(z);
        for i = 1:2
            JuMP.setRHS(zcon[i], zv[i])
        end
        println("slaveLP problem is:\n", slaveLP)
        sstatus = solve(slaveLP);
        assert(sstatus == :Optimal)

        zcond = getdual(zconLP);
        # Wcond = getdual(wconLP);
        bendersCutLHS = sum(zcond .* (1-x))
        # + (-Wcond[1]*sum(wLP) + Wcond[2:end]' * yLP)
        println("Benders cut ", (bendersCutLHS - target + bias))
        @lazyconstraint(cb, bendersCutLHS >= target - bias)
        println("\nMaster program is \n",master)

        println("value of x = ", getvalue(x))
    end

    # addlazycallback(master, addBendersCut);
    # solve(master)

    while j < 3
        # Solving Slave MIP
        println("j = ",j)
        println("slave problem is:\n", slave)
        sstatus = solve(slave);
        assert(sstatus == :Optimal)
        objv = getobjectivevalue(slave)
        println("slave MIP obj = ", objv)
        if objv >= target
            break;
        end

        # Solving Slave LP
        zv = getvalue(z);
        for i = 1:2
            JuMP.setRHS(zcon[i], zv[i])
        end
        println("slaveLP problem is:\n", slaveLP)
        sstatus = solve(slaveLP);
        assert(sstatus == :Optimal)

        zcond = getdual(zconLP);
        # Wcond = getdual(wconLP);
        bendersCutLHS = sum(zcond .* (1-x))
        # + (-Wcond[1]*sum(wLP) + Wcond[2:end]' * yLP)

        @constraint(master, bendersCutLHS >= target - bias)
        println("\nMaster program is \n",master)
        mstatus = solve(master);
        assert(mstatus == :Optimal)
        xv = getvalue(x)
        println("value of x = ", xv)

        for i = 1:2
            JuMP.setRHS(xcon[i], 1-xv[i])
        end
        j+=1

    end

    println("value of x = ", getvalue(x))
    println("value of w = ", getvalue(wLP))
    println("value of y = ", getvalue(yLP))
    println("value of z = ", getvalue(z))
end
