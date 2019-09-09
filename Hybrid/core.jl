function workerPrinterCheck()
    f = @spawn myid()
    println("====")
    println(fetch(f))
end

function returnArrayOne()
    return [myid(), myid() + 1]
end

function changeToOne(cutest)
    cutest[1] = 8
    cutest[2] = 9
    return nothing
end

function check_leaks()
    if length(DistributedArrays.refs) > 0
        sleep(1)  # allow time for any cleanup to complete and test again
        length(DistributedArrays.refs) > 0 && @warn("Probable leak of ", length(DistributedArrays.refs), " darrays")
    end
end

function writeToFile(loc_globaldata, numPoints)
    file = open("results/split/primvals"* "_" * string(getConfig()["core"]["max_iters"]) *
     "_" * string(numPoints) * "_" * string(myid()) * ".txt", "w")
    @showprogress 1 "This takes time" for (idx, _) in enumerate(loc_globaldata)
        primtowrite = loc_globaldata[idx].prim
        for element in primtowrite
            @printf(file,"%0.17f", element)
            @printf(file, " ")
        end
        print(file, "\n")
    end
    close(file)
end

@inline function updateLocalGhost(loc_ghost_holder, dist_globaldata)
    localkeys = keys(loc_ghost_holder[1])
    # println(localkeys)
    for iter in localkeys
        #Dict To Array Equality
        loc_ghost_holder[1][iter] = dist_globaldata[loc_ghost_holder[1][iter].globalID]
    end
    return nothing
end

@inline function updateLocalGhostQ(loc_ghost_holder, dist_q)
    localkeys = keys(loc_ghost_holder[1])
    # println(localkeys)
    for iter in localkeys
        #Dict To Array Equality
        @. loc_ghost_holder[1][iter].q = dist_q[loc_ghost_holder[1][iter].globalID].q
        # @. loc_ghost_holder_mutable[1][iter][9:12] = loc_ghost_holder[1][iter].q
    end
    return nothing
end

@inline function updateLocalGhostDQ(loc_ghost_holder, dist_dq)
    localkeys = keys(loc_ghost_holder[1])
    # println(localkeys)
    for iter in localkeys
        #Dict To Array Equality
        @. loc_ghost_holder[1][iter].dq = dist_dq[loc_ghost_holder[1][iter].globalID].dq
        # @. loc_ghost_holder_mutable[1][iter][13:16] = loc_ghost_holder[1][iter].dq[1]
        # @. loc_ghost_holder_mutable[1][iter][17:20] = loc_ghost_holder[1][iter].dq[2]
    end
    return nothing
end

function getInitialPrimitive(configData)
    rho_inf = configData["core"]["rho_inf"]::Float64
    mach = configData["core"]["mach"]::Float64
    machcos::Float64 = mach * cos(calculateTheta(configData))
    machsin::Float64 = mach * sin(calculateTheta(configData))
    pr_inf = configData["core"]["pr_inf"]::Float64
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

function fpi_solver(iter_store, ghost_holder, dist_globaldata, dist_q, dist_dq, dist_globaldata_mutable, ghost_holder_mutable, configData, res_old, numPoints)
    # println(IOContext(stdout, :compact => false), globaldata[3].prim)
    # print(" 111\n")
    str = CuStream()
    power = configData["core"]["power"]
    cfl = configData["core"]["cfl"]::Float64
    Mach::Float64 = configData["core"]["mach"]::Float64
    gamma::Float64 = configData["core"]["gamma"]::Float64
    pr_inf::Float64 = configData["core"]["pr_inf"]::Float64
    rho_inf::Float64 = configData["core"]["rho_inf"]::Float64
    theta = calculateTheta(configData)
    threadsperblock = parse(Int , ARGS[4])

    for iter in 1:iter_store
        if iter == 1
            println("Starting FuncDelta")
        end
    # @sync

        @sync for ip in procs(dist_globaldata)
            @spawnat ip begin
                updateLocalGhost(ghost_holder[:L], dist_globaldata)
            end
        end

        @sync for ip in procs(dist_globaldata)
            @spawnat ip begin
                # println(length(localpart(globaldata)))
                func_delta(dist_globaldata[:L], ghost_holder[:L], cfl, numPoints)
            end
        end

        for rk in 1:4
        #    # if iter == 1
        #        # println("Starting QVar")
        #    # end
            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    q_variables(dist_globaldata[:L], dist_q[:L], dist_globaldata_mutable[:L])
                end
            end

            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    updateLocalGhostQ(ghost_holder[:L], dist_q)
                end
            end

            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    q_var_derivatives(dist_globaldata[:L], dist_dq[:L], ghost_holder[:L], dist_globaldata_mutable[:L], power)
                end
            end
        # println(IOContext(stdout, :compact => false), dist_globaldata[3])
        #    # if iter == 1
        #        # println("Starting Calflux")
        #    # end

            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    updateLocalGhostDQ(ghost_holder[:L], dist_dq)
                end
            end

            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    cal_flux_residual(dist_globaldata[:L], ghost_holder[:L], dist_globaldata_mutable[:L], configData)
                end
            end

            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    gpuLocDataRest = CuArray(dist_globaldata_mutable[:L])
                    locBlocksPerGrid = Int(ceil(gpuLocNumPoints/threadsperblock))
                    @cuda blocks= locBlocksPerGrid threads= threadsperblock cal_flux_residual_kernel(gpuLocDataConn, gpuLocDataFixedPoint, gpuLocDataRest, gpuConfigData, gpuLocNumPoints)
                    dist_globaldata_mutable[:L] = Array(gpuLocDataRest)
                end
            end
        #    # println(IOContext(stdout, :compact => false), globaldata[3].prim)
        #    # println(IOContext(stdout, :compact => false), globaldata[3].prim)
        #    # residue = 0
        #    # if iter == 1
        #        # println("Starting StateUpdate")
        #    # end
            @sync for ip in procs(dist_globaldata)
                @spawnat ip begin
                    state_update(dist_globaldata[:L], dist_globaldata_mutable[:L], Mach, gamma, pr_inf, rho_inf, theta, iter, res_old, rk, numPoints)
                end
            end
        end
        println("Iteration Number ", iter)
    end
    synchronize()
    # println(IOContext(stdout, :compact => false), globaldata[3].prim)
    # residue = res_old
    return nothing
end

@inline function q_variables(loc_globaldata, loc_q, loc_globaldata_mutable)
    for (idx, itm) in enumerate(loc_globaldata)
        rho = itm.prim[1]
        u1 = itm.prim[2]
        u2 = itm.prim[3]
        pr = itm.prim[4]

        beta = 0.5 * (rho / pr)
        loc_globaldata_mutable[9, idx] = log(rho) + log(beta) * 2.5 - (beta * ((u1 * u1) + (u2 * u2)))
        two_times_beta = 2.0 * beta

        loc_globaldata_mutable[10, idx] = (two_times_beta * u1)
        loc_globaldata_mutable[11, idx] = (two_times_beta * u2)
        loc_globaldata_mutable[12, idx] = -two_times_beta
        loc_globaldata[idx].q = @views loc_globaldata_mutable[9:12, idx]
        loc_q[idx].q = loc_globaldata[idx].q
    end
    return nothing
end

function q_var_derivatives(loc_globaldata, loc_dq, loc_ghost_holder, loc_globaldata_mutable, power)
    sum_delx_delq = zeros(Float64, 4)
    sum_dely_delq = zeros(Float64, 4)
    dist_length = length(loc_globaldata)

    for (idx, itm) in enumerate(loc_globaldata)
        x_i = itm.x
        y_i = itm.y
        sum_delx_sqr = zero(Float64)
        sum_dely_sqr = zero(Float64)
        sum_delx_dely = zero(Float64)
        sum_delx_delq = fill!(sum_delx_delq, 0.0)
        sum_dely_delq = fill!(sum_dely_delq, 0.0)
        for i in 1:4
            loc_globaldata[idx].max_q[i] = loc_globaldata[idx].q[i]
            loc_globaldata_mutable[20 + i, idx] = loc_globaldata[idx].q[i]
            loc_globaldata[idx].min_q[i] = loc_globaldata[idx].q[i]
            loc_globaldata_mutable[24 + i, idx] = loc_globaldata[idx].q[i]
        end
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
            sum_delx_sqr += ((delx * delx) * weights)
            sum_dely_sqr += ((dely * dely) * weights)
            sum_delx_dely += ((delx * dely) * weights)

            for i in 1:4
                sum_delx_delq[i] += (weights * delx * (globaldata_conn.q[i] - loc_globaldata[idx].q[i]))
                sum_dely_delq[i] += (weights * dely * (globaldata_conn.q[i] - loc_globaldata[idx].q[i]))
                if loc_globaldata_mutable[20 + i, idx] < globaldata_conn.q[i]
                    loc_globaldata[idx].max_q[i] = globaldata_conn.q[i]
                    loc_globaldata_mutable[20 + i, idx] = globaldata_conn.q[i]
                end
                if loc_globaldata_mutable[24 + i, idx] > globaldata_conn.q[i]
                    loc_globaldata[idx].min_q[i] = globaldata_conn.q[i]
                    loc_globaldata_mutable[24 + i, idx] = globaldata_conn.q[i]
                end
            end
        end
        det = (sum_delx_sqr * sum_dely_sqr) - (sum_delx_dely * sum_delx_dely)
        one_by_det = 1.0 / det
        loc_globaldata_mutable[13:16, idx] = @. one_by_det * (sum_delx_delq * sum_dely_sqr - sum_dely_delq * sum_delx_dely)
        loc_globaldata[idx].dq[1] = @views loc_globaldata_mutable[13:16, idx]
        loc_globaldata_mutable[17:20, idx] = @. one_by_det * (sum_dely_delq * sum_delx_sqr - sum_delx_delq * sum_delx_dely)
        loc_globaldata[idx].dq[2] = @views loc_globaldata_mutable[17:20, idx]
        loc_dq[idx].dq[1] = loc_globaldata[idx].dq[1]
        loc_dq[idx].dq[2] = loc_globaldata[idx].dq[2]
        # globaldata[idx].dq = [tempsumx, tempsumy]
    end
    # println(IOContext(stdout, :compact => false), globaldata[3].dq)
    # println(IOContext(stdout, :compact => false), globaldata[3].max_q)
    # println(IOContext(stdout, :compact => false), globaldata[3].min_q)
    return nothing
end


# @inline function update_ghost_q_variables(ghost_point)
#     rho = ghost_point.prim[1]
#     u1 = ghost_point.prim[2]
#     u2 = ghost_point.prim[3]
#     pr = ghost_point.prim[4]

#     beta = 0.5 * (rho / pr)
#     ghost_point.q[1] = log(rho) + log(beta) * 2.5 - (beta * ((u1 * u1) + (u2 * u2)))
#     two_times_beta = 2.0 * beta

#     ghost_point.q[2] = (two_times_beta * u1)
#     ghost_point.q[3] = (two_times_beta * u2)
#     ghost_point.q[4] = -two_times_beta
#     return nothing
# end

# function updateLocalGhostNew(loc_ghost_holder, dist_globaldata)
#     localkeys = keys(loc_ghost_holder[1])
#     # println(localkeys)
#     for iter in localkeys
#         #Dict To Array Equality
#         loc_ghost_holder[1][iter].prim = dist_globaldata[loc_ghost_holder[1][iter].globalID].prim
#         loc_ghost_holder[1][iter].q = dist_globaldata[loc_ghost_holder[1][iter].globalID].q
#     end
#     return nothing
# end

# function updateGhost(ghost_holder, dist_globaldata)
#     for i in 1:nworkers()
#         localkeys = keys(ghost_holder[i])
#         for iter in localkeys
#             ghost_holder[i][iter] = dist_globaldata[ghost_holder[i][iter].globalID]
#         end
#     end
# end