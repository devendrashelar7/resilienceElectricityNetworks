using JuMP

begin

    function hfunc(i, V, Y, S)

        rhs = conj(S[i])/conj(V[i])
        for k = 1:n
            if k != i
                rhs -= Y[i,k] * V[k]
            end
        end
        rhs /= Y[i,i]
    end

    n = 3
    V = (1+0im) * ones(n)
    Y = (0+0im) * ones(n,n) # admittance matrix
    p = 0.01+0.01im
    p /= 100
    S = [2p, -p, -p] # power at buses
    yp = 1+1im

    # generate the admittance matrix
    for i = 1:n
        for j = 1:n
            if i == j
                Y[i,j] = 2.01yp
            else
                Y[i,j] = -yp
            end
        end
    end

    for t = 1:3 # iterations
        for i = 2:n
            V[i] = hfunc(i, V, Y, S)  # V_i = h_i(V_i)
        end
        println(V)
    end

    println(abs.(V))
    println(angle.(V))

end
