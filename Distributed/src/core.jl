function workerPrinterCheck()
    f = @spawn myid()
    println("====")
    println(fetch(f))
end

function check_leaks()
    if length(DistributedArrays.refs) > 0
        sleep(1)  # allow time for any cleanup to complete and test again
        length(DistributedArrays.refs) > 0 && @warn("Probable leak of ", length(DistributedArrays.refs), " darrays")
    end
end

@inline function updateLocalGhost(loc_ghost_holder, dist_globaldata)
    localkeys = keys(loc_ghost_holder[1])
    # println(localkeys)
    for iter in localkeys
        #Dict To Array Equality
        loc_ghost_holder[1][iter] = dist_globaldata[loc_ghost_holder[1][iter].globalID]
        # if iter == 6305000
        #     print("Ghost>> ")
        #     println(loc_ghost_holder[1][iter])
        # end
    end
    return nothing
end

@inline function updateLocalGhostQ(loc_ghost_holder, dist_q)
    localkeys = keys(loc_ghost_holder[1])
    # println(localkeys)
    for iter in localkeys
        #Dict To Array Equality
        @. loc_ghost_holder[1][iter].q = dist_q[loc_ghost_holder[1][iter].globalID].q
        # if iter == 63050000
        #     print("Q>> ")
        #     println(loc_ghost_holder[1][iter])
        # end
    end
    return nothing
end

@inline function updateLocalGhostQPack(loc_ghost_holder, dist_qpack)
    localkeys = keys(loc_ghost_holder[1])
    # println(localkeys)
    for iter in localkeys
        #Dict To Array Equality
        merge_holder = dist_qpack[loc_ghost_holder[1][iter].globalID]
        for idx in 1:4
            loc_ghost_holder[1][iter].dq1[idx] = merge_holder.dq1[idx]
            loc_ghost_holder[1][iter].dq2[idx] = merge_holder.dq2[idx]
            loc_ghost_holder[1][iter].max_q[idx] = merge_holder.max_q[idx]
            loc_ghost_holder[1][iter].min_q[idx] = merge_holder.min_q[idx]
        end
    end
    return nothing
end

@inline function updateLocalGhostDQ(loc_ghost_holder, dist_qpack)
    localkeys = keys(loc_ghost_holder[1])
    for iter in localkeys
        merge_holder = dist_qpack[loc_ghost_holder[1][iter].globalID]
        for idx in 1:4
            loc_ghost_holder[1][iter].dq1[idx] = merge_holder.dq1[idx]
            loc_ghost_holder[1][iter].dq2[idx] = merge_holder.dq2[idx]
        end
    end
    return nothing
end

@inline function updateLocalGhostPrim(loc_ghost_holder, dist_prim)
    localkeys = keys(loc_ghost_holder[1])
    # println(localkeys)
    for iter in localkeys
        @. loc_ghost_holder[1][iter].prim = dist_prim[loc_ghost_holder[1][iter].globalID].prim
    end
    return nothing
end

function getInitialPrimitive(configData)
    rho_inf = configData["core"]["rho_inf"]
    mach = configData["core"]["mach"]
    machcos::Float64 = mach * cos(calculateTheta(configData))
    machsin::Float64 = mach * sin(calculateTheta(configData))
    pr_inf = configData["core"]["pr_inf"]
    primal = [rho_inf, machcos, machsin, pr_inf]
    return primal
end

function getInitialPrimitive2(configData)
    dataman = open("prim_soln_clean")
    data = read(dataman, String)
    data1 = split(data, "\n")
    finaldata = Array{Array{Float64,1},1}(undef, 0)
    for (idx,itm) in enumerate(data1)
        # try
        da = split(itm)
        da1 = parse.(Float64, da)
        push!(finaldata, da1)
    end
    close(dataman)
    return finaldata
end

function matchInitialQ(loc_q, loc_globaldata)
    local_size = length(loc_q)
    for idx in 1:local_size
        loc_q[idx].q = loc_globaldata[idx].q
    end
end

function placeNormals(loc_globaldata, globaldata, loc_ghost_holder, interior, wall, outer)
    local_size = length(loc_globaldata)
    # println(local_size, " is the local size")
    updateLocalGhost(loc_ghost_holder, globaldata)
    for idx in 1:local_size
        flag = loc_globaldata[idx].flag_1
        if flag == wall || flag == outer
            currpt = getxy(loc_globaldata[idx])
            leftpt = loc_globaldata[idx].left
            rightpt = loc_globaldata[idx].right

            if leftpt <= local_size
                leftpt = getxy(loc_globaldata[leftpt])
            else
                leftpt = getxy(loc_ghost_holder[1][leftpt])
            end

            if rightpt <= local_size
                rightpt = getxy(loc_globaldata[rightpt])
            else
                rightpt = getxy(loc_ghost_holder[1][rightpt])
            end

            normals = calculateNormals(leftpt, rightpt, currpt[1], currpt[2])
            setNormals(loc_globaldata[idx], normals)
        elseif flag == interior
            setNormals(loc_globaldata[idx], (0,1))
        else
            @warn "Illegal Point Type"
        end
    end
    return nothing
end

function calculateNormals(left, right, mx, my)
    lx = left[1]
    ly = left[2]

    rx = right[1]
    ry = right[2]

    nx1 = my - ly
    nx2 = ry - my

    ny1 = mx - lx
    ny2 = rx - mx

    nx = 0.5*(nx1 + nx2)
    ny = 0.5*(ny1 + ny2)

    det = hypot(nx, ny)

    nx = -nx/det
    ny = ny/det

    return (nx,ny)
end

function calculateConnectivity(loc_globaldata, globaldata, loc_ghost_holder)
    local_size = length(loc_globaldata)

    for idx in 1:local_size
        ptInterest = loc_globaldata[idx]
        currx = ptInterest.x
        curry = ptInterest.y
        nx = ptInterest.nx
        ny = ptInterest.ny

        flag = ptInterest.flag_1

        xpos_conn,xneg_conn,ypos_conn,yneg_conn = Array{Int32,1}(undef, 0),Array{Int32,1}(undef, 0),Array{Int32,1}(undef, 0),Array{Int32,1}(undef, 0)

        tx = ny
        ty = -nx

        for itm in ptInterest.conn
            if itm <= local_size
                itmx = loc_globaldata[itm].x
                itmy = loc_globaldata[itm].y
            else
                itmx = loc_ghost_holder[1][itm].x
                itmy = loc_ghost_holder[1][itm].y
            end

            delx = itmx - currx
            dely = itmy - curry

            dels = delx*tx + dely*ty
            deln = delx*nx + dely*ny
            if dels <= 0.0
                push!(xpos_conn, itm)
            end
            if dels >= 0.0
                push!(xneg_conn, itm)
            end
            if flag == 1
                if deln <= 0.0
                    push!(ypos_conn, itm)
                end
                if deln >= 0.0
                    push!(yneg_conn, itm)
                end
            elseif flag == 0
                push!(yneg_conn, itm)
            elseif flag == 2
                push!(ypos_conn, itm)
            end
        end
        setConnectivity(loc_globaldata[idx], (xpos_conn, xneg_conn, ypos_conn, yneg_conn))
    end
    return nothing
end

function fpi_solver(iter, ghost_holder, dist_globaldata, dist_q, dist_qpack, dist_prim, res_old, res_new, numPoints, main_store)
    # println(IOContext(stdout, :compact => false), globaldata[3].prim)
    # print(" 111\n")
    power = main_store[53]
    cfl = main_store[54]

    @sync for ip in procs(dist_globaldata)
        @spawnat ip begin
            updateLocalGhostPrim(ghost_holder[:L], dist_prim)
        end
    end

    if iter == 1
        println("Starting FuncDelta")
        @sync for ip in procs(dist_globaldata)
            @spawnat ip begin
                updateLocalGhost(ghost_holder[:L], dist_globaldata)
            end
        end
    end

    @sync for ip in procs(dist_globaldata)
        @spawnat ip begin
            func_delta(dist_globaldata[:L], ghost_holder[:L], cfl)
        end
    end

    phi_i = @view main_store[1:4]
	phi_k = @view main_store[5:8]
	G_i = @view main_store[9:12]
    G_k = @view main_store[13:16]
	result = @view main_store[17:20]
	qtilde_i = @view main_store[21:24]
	qtilde_k = @view main_store[25:28]
	Gxp = @view main_store[29:32]
	Gxn = @view main_store[33:36]
	Gyp = @view main_store[37:40]
	Gyn = @view main_store[41:44]
    ∑_Δx_Δf = @view main_store[45:48]
    ∑_Δy_Δf = @view main_store[49:52]

    for rk in 1:4
    #    # if iter == 1
    #        # println("Starting QVar")
    #    # end

        @sync for ip in procs(dist_globaldata)
            @spawnat ip begin
                q_variables(dist_globaldata[:L], dist_q[:L], result)
                updateLocalGhostQ(ghost_holder[:L], dist_q)
            end
        end

        @sync for ip in procs(dist_globaldata)
            @spawnat ip begin
                q_var_derivatives(dist_globaldata[:L], dist_qpack[:L], ghost_holder[:L], power, ∑_Δx_Δf, ∑_Δy_Δf, qtilde_i, qtilde_k)
                updateLocalGhostQPack(ghost_holder[:L], dist_qpack)
            end
        end

        for inner_iters in 1:3
            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    q_var_derivatives_innerloop(dist_globaldata[:L], dist_qpack[:L], ghost_holder[:L], power, ∑_Δx_Δf, ∑_Δy_Δf, qtilde_i, qtilde_k)
                    updateLocalGhostDQ(ghost_holder[:L], dist_qpack)
                end
            end
        end

        @sync for ip in procs(dist_globaldata)
            @spawnat ip begin
                cal_flux_residual(dist_globaldata[:L], dist_globaldata, ghost_holder[:L], Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k,
                result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, main_store)
            end
        end

        @sync for ip in procs(dist_globaldata)
            @spawnat ip begin
                state_update(dist_globaldata[:L], dist_prim[:L], iter, res_old[:L], res_new[:L], rk, ∑_Δx_Δf, ∑_Δy_Δf, main_store)
            end
        end

    end

    residue = 0
    if iter <= 2
        # res_old[1] = res_new[1]
        residue = 0
    else
        residue = log10(sqrt(sum(res_new))/sqrt(sum(res_old)))
    end

    println("Iteration Number ", iter, " ", residue)
    # println(IOContext(stdout, :compact => false), globaldata[3].prim)
    # residue = res_old
    return nothing
end

@inline function q_variables(loc_globaldata, loc_q, q_result)
    for (idx, itm) in enumerate(loc_globaldata)
        # if idx == 2000
        #     println("==========================================")
        # end
        rho = itm.prim[1]
        u1 = itm.prim[2]
        u2 = itm.prim[3]
        pr = itm.prim[4]
        beta = 0.5 * (rho / pr)
        two_times_beta = 2.0 * beta
        q_result[1] = log(rho) + log(beta) * 2.5 - (beta * ((u1 * u1) + (u2 * u2)))
        q_result[2] = two_times_beta * u1
        q_result[3] = two_times_beta * u2
        q_result[4] = -two_times_beta

        @. loc_globaldata[idx].q = q_result
        @. loc_q[idx].q = loc_globaldata[idx].q
    end
    return nothing
end

function q_var_derivatives(loc_globaldata, loc_qpack, loc_ghost_holder, power, ∑_Δx_Δq, ∑_Δy_Δq, max_q, min_q)
    dist_length = length(loc_globaldata)

    for (idx, itm) in enumerate(loc_globaldata)
        x_i = itm.x
        y_i = itm.y
        ∑_Δx_sqr = zero(Float64)
        ∑_Δy_sqr = zero(Float64)
        ∑_Δx_Δy = zero(Float64)
        fill!(∑_Δx_Δq, zero(Float64))
        fill!(∑_Δy_Δq, zero(Float64))
        
        @. max_q = loc_globaldata[idx].q
        @. min_q = loc_globaldata[idx].q

        for conn in itm.conn
            if conn <= dist_length
                globaldata_conn = loc_globaldata[conn]
            else
                # update_ghost_q_variables(loc_ghost_holder[1][conn])
                globaldata_conn = loc_ghost_holder[1][conn]
            end
            x_k = globaldata_conn.x
            y_k = globaldata_conn.y
            delx = x_k - x_i
            dely = y_k - y_i
            dist = hypot(delx, dely)
            weights = dist ^ power
            ∑_Δx_sqr += ((delx * delx) * weights)
            ∑_Δy_sqr += ((dely * dely) * weights)
            ∑_Δx_Δy += ((delx * dely) * weights)

            for i in 1:4
                ∑_Δx_Δq[i] += (weights * delx * (globaldata_conn.q[i] - loc_globaldata[idx].q[i]))
                ∑_Δy_Δq[i] += (weights * dely * (globaldata_conn.q[i] - loc_globaldata[idx].q[i]))
                if max_q[i] < globaldata_conn.q[i]
                    max_q[i] = globaldata_conn.q[i]
                end
                if min_q[i] > globaldata_conn.q[i]
                    min_q[i] = globaldata_conn.q[i]
                end
            end
        end
        @. loc_globaldata[idx].max_q = max_q
        @. loc_globaldata[idx].min_q = min_q
        @. loc_qpack[idx].max_q = loc_globaldata[idx].max_q
        @. loc_qpack[idx].min_q = loc_globaldata[idx].min_q
        det = (∑_Δx_sqr * ∑_Δy_sqr) - (∑_Δx_Δy * ∑_Δx_Δy)
        one_by_det = 1.0 / det
        @. loc_globaldata[idx].dq1 = one_by_det * (∑_Δx_Δq * ∑_Δy_sqr - ∑_Δy_Δq * ∑_Δx_Δy)
        @. loc_globaldata[idx].dq2 = one_by_det * (∑_Δy_Δq * ∑_Δx_sqr - ∑_Δx_Δq * ∑_Δx_Δy)
        @. loc_qpack[idx].dq1 = loc_globaldata[idx].dq1
        @. loc_qpack[idx].dq2 = loc_globaldata[idx].dq2
    end
    return nothing
end

function q_var_derivatives_innerloop(loc_globaldata, loc_qpack, loc_ghost_holder, power, ∑_Δx_Δq, ∑_Δy_Δq, qi_tilde, qk_tilde)
    dist_length = length(loc_globaldata)

    for (idx, itm) in enumerate(loc_globaldata)
        x_i = itm.x
        y_i = itm.y
        ∑_Δx_sqr = zero(Float64)
        ∑_Δy_sqr = zero(Float64)
        ∑_Δx_Δy = zero(Float64)
        fill!(∑_Δx_Δq, zero(Float64))
        fill!(∑_Δy_Δq, zero(Float64))

        for conn in itm.conn
            if conn <= dist_length
                globaldata_conn = loc_globaldata[conn]
            else
                globaldata_conn = loc_ghost_holder[1][conn]
            end
            x_k = globaldata_conn.x
            y_k = globaldata_conn.y
            delx = x_k - x_i
            dely = y_k - y_i
            dist = hypot(delx, dely)
            weights = dist ^ power
            ∑_Δx_sqr += ((delx * delx) * weights)
            ∑_Δy_sqr += ((dely * dely) * weights)
            ∑_Δx_Δy += ((delx * dely) * weights)

            for iter in 1:4
                qi_tilde[iter] = loc_globaldata[idx].q[iter] - 0.5 * (delx * loc_globaldata[idx].dq1[iter] + dely * loc_globaldata[idx].dq2[iter])
                qk_tilde[iter] = globaldata_conn.q[iter] - 0.5 * (delx * globaldata_conn.dq1[iter] + dely * globaldata_conn.dq2[iter])
                intermediate_var = weights * (qk_tilde[iter] - qi_tilde[iter])
                ∑_Δx_Δq[iter] += delx * intermediate_var
                ∑_Δy_Δq[iter] += dely * intermediate_var
            end
        end
        det = (∑_Δx_sqr * ∑_Δy_sqr) - (∑_Δx_Δy * ∑_Δx_Δy)
        one_by_det = 1.0 / det
        # for iter in 1:4
        @. loc_globaldata[idx].tempdq1 = one_by_det * (∑_Δx_Δq * ∑_Δy_sqr - ∑_Δy_Δq * ∑_Δx_Δy)
        @. loc_globaldata[idx].tempdq2 = one_by_det * (∑_Δy_Δq * ∑_Δx_sqr - ∑_Δx_Δq * ∑_Δx_Δy)
        # end 
    end
    for (idx, _) in enumerate(loc_globaldata)
        @. loc_globaldata[idx].dq1 = loc_globaldata[idx].tempdq1
        @. loc_globaldata[idx].dq2 = loc_globaldata[idx].tempdq2
        @. loc_qpack[idx].dq1 = loc_globaldata[idx].dq1
        @. loc_qpack[idx].dq2 = loc_globaldata[idx].dq2
    end
    return nothing
end