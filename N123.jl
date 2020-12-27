using PyCall
ENV["PYTHON"]=""
    Pkg.build("PyCall")

Pkg.add("ExcelFiles")
Pkg.add("ExcelReaders")
using PyCall
using ExcelFiles, ExcelReaders
Pkg.add("CSV")
Pkg.add("DataFrames")
using CSV
using DataFrames

begin
    N = 123; setGlobalParameters(); filename = string("N", N);
    global printPerformance = true;
    ResWC = 80;
    LLCreq = (1-ResWC/100) * LLCmax;
    # delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getBendersMethodCascade(LLCreq);
    delta, deltaV0, pv, prv, qv, qrv, betav, Pv, Qv, vv, tv, kcv, kgv, pcv, qcv, objv, pgv, qgv, lacv, lovrv, lcontv, lsdv = getLinearAssistBendersCascade(LLCreq);
    println("Cardinality: ", sum(delta));
end


begin
    # df = CSV.read("N123.csv",null="\N");
    df = CSV.readtable("N123.csv")
    # println(df[1:4,1:end])
    # println(df[:R])
    resistance = df[:R];
    reactance = df[:R];
    P = (df[:P])[1:85]; Q = (df[:Q])[1:85];
    loadNodes = (df[:Node])[1:85];
    loadNodes = map(loadNodes->parse(Int64,loadNodes),loadNodes)
    originalNodes = zeros(1:450);
    nodesInDoc = df[:nodeNoInDoc]; nodesInJulia = df[:nodeNoInJulia];
    to = df[:to]; from = df[:from];
    for i = 1:124
        j = nodesInDoc[i];
        corr = nodesInJulia[i];
        originalNodes[j] = corr;
        # println("i ", i, " j ",j," corr ",corr);
    end

    # println(originalNodes[150]);
    # println(length(nodesInJulia))
    nNodes = 123; N = nNodes;
    nodes=collect(1:nNodes);
    par = zeros(Int64,nNodes);
    # println("nNodes ",nNodes);
    # println("lfrom ",length(from))

    for i = 1:122
        fromNode = from[i]; toNode = to[i];
        ofromNode = originalNodes[fromNode];
        otoNode = originalNodes[toNode];
        # println("i ",i," ",fromNode, " ", toNode," ", ofromNode," ",otoNode);
        par[otoNode] = ofromNode;
    end
    # println(length(nodesInJulia));
    # println(par)
    # println(length(par))
    # println(par[originalNodes[149]])

    # println(length(nodes[par.==0]))
    # println(nodes[par.==0])
    # println(par[nodes[par.==0]])
    # println(length(P))
    # println(typeof(P))
    pcmax = zeros(Float64,N); qcmax = zeros(Float64,N);
    # println(loadNodes)
    pcmax[originalNodes[loadNodes]] = P;
    qcmax[originalNodes[loadNodes]] = Q;
    # println(pcmax.*qcmax)
    # println(nodes[P.=="NA"])
    # println(P[124]," ",P[125])
    # println(originalNodes[par.==0])
    # println(par[nodes].==0)
    # resistance = convert(Array{Float64},resistance)
    # println(typeof(resistance))
    # a = float(resistance[1:4]).^2 + 4;
    # println(resistance[1:4])
    # println(reactance[1:4])
    # println(P[1:4].^4+2);
    # # acp = map(acp->parse(Float64,acp),acp[1:85])
    # # acp = map(acp->parse(Float64,acq),acp[1:85])
    # println(length(P));
    # println(acp[1:86])
    # println(acq[1:4]);
end

Pkg.dependents("ExcelReaders")
Pkg.installed()
using DataFrames

df = df = readxl(DataFrame, "123.xlsx", "Nodes!A1:C4")

Pkg.add("xlrd")

begin
  # Pkg.update()
  # Pkg.add("PyCall")
  #
  # Pkg.add("ExcelReaders")

  println(pwd())
  # Pkg.test("ExcelReaders")
  # f = openxl(string(pwd(),"/123.xlsx"));
  # openxl("/Users/devendrashelar/Dropbox (MIT)/DNCS/papers/secondJournal/code/123.xlsx")
  # @pyimport math
  # math.sin(math.pi / 4) - sin(pi / 4)
  f = openxl("N123.xlsx")

  rNodes = readxl(f, "nodes!A2:A124");
  # data = readxl("N123.xlsx", "Nodes!A1:A4")
  # data = readxlsheet("N123.xlsx", "Nodes")
  # println(length(rNodes))
  # println(VERSION)
end
