function main()

    configData = getConfig()

    format = configData["format"]["type"]
    # file_name = string(ARGS[2])
    folder_name = string(ARGS[3])

    numPoints = parse(Int, ARGS[2])
    if format == "old"
        numPoints += 1
    end
    # numPoints = 39381464
    println(numPoints)
    # globaldata = Array{Point,1}(undef, numPoints)
    # table = Array{Int32,1}(undef, numPoints)
    # wallptsidx = Array{Int32,1}(undef, 0)
    # Interiorptsidx = Array{Int32,1}(undef, 0)
    # outerptsidx = Array{Int32,1}(undef, 0)
    # shapeptsidx = Array{Int32,1}(undef, 0)
    # print(splitdata[1:3])
    main_store = zeros(Float64, 62)
    defprimal = getInitialPrimitive(configData)

    # local_points_holder = Array{Array{Point,1},1}(undef, nworkers())
    ghost_holder = Array{Dict{Int64,Point},1}(undef, nworkers())
    # q_ghost_holder = Array{Dict{Int64,TempQ},1}(undef, nworkers())
    # dq_ghost_holder = Array{Dict{Int64,TempDQ},1}(undef, nworkers())
    # prim_ghost_holder = Array{Dict{Int64,TempPrim},1}(undef, nworkers())

    global_local_map_index = Dict{Tuple{Float64, Float64},Int64}()
    global_local_direct_index = Array{Int64, 1}(undef, numPoints)
    index_holder = dzeros(Int, nworkers())
    println("Indexing")

    createGlobalLocalMapIndex(global_local_map_index, global_local_direct_index, index_holder, folder_name::String)
    # println(global_local_map_index)

    println("Start Read")
    # count = 0
    # readFile(file_name::String, globaldata, table, defprimal, wallptsidx, outerptsidx, Interiorptsidx, shapeptsidx,
        # wallpts, Interiorpts, outerpts, shapepts, numPoints)

    # println("Read Local Files")
    # readDistribuedFile(folder_name::String, defprimal, 12, global_local_map_index)
    println("Reading multiple files")
    if format == "quadtree"
        globaldata_parts = [@spawnat p readDistribuedFileQuadtree(folder_name::String, defprimal, p, index_holder) for p in workers()]
    elseif format == "old"
        globaldata_parts = [@spawnat p readDistribuedFile(folder_name::String, defprimal, p, global_local_map_index) for p in workers()]
    end
    globaldata_parts = reshape(globaldata_parts, (nworkers()))
    dist_globaldata = DArray(globaldata_parts)

    println("Reading multiple files for Q")
    q_parts = [@spawnat p readDistribuedFileQ(folder_name::String, defprimal, p, global_local_map_index) for p in workers()]
    q_parts = reshape(q_parts, (nworkers()))
    dist_q = DArray(q_parts)

    println("Reading multiple files for QPack")
    dq_parts = [@spawnat p readDistribuedFileQPack(folder_name::String, defprimal, p, global_local_map_index) for p in workers()]
    dq_parts = reshape(dq_parts, (nworkers()))
    dist_qpack = DArray(dq_parts)

    println("Reading Ghost")
    readGhostFile(folder_name, ghost_holder, global_local_map_index, dist_globaldata)

    ghost_holder = distribute(ghost_holder, procs=workers(), dist=(length(workers()),))
    keys_holder = [@spawnat p returnKeys(ghost_holder[:L]) for p in workers()]
    keys_holder = reshape(keys_holder, (nworkers()))
    dist_keys = DArray(keys_holder)


    println("Reading Normals")
    interior = configData["point"]["interior"]
    wall = configData["point"]["wall"]
    outer = configData["point"]["outer"]
    @sync for pid in procs(dist_globaldata)
        @spawnat pid begin
            placeNormals(dist_globaldata[:L], dist_globaldata, ghost_holder[:L], dist_keys[:L], interior, wall, outer)
        end
    end

    println("Start table sorting")
    @sync for pid in procs(dist_globaldata)
        @spawnat pid begin
            calculateConnectivity(dist_globaldata[:L], ghost_holder[:L])
        end
    end

    main_store[53] = configData["core"]["power"]::Float64
    main_store[54] = configData["core"]["cfl"]::Float64 
    main_store[55] = configData["core"]["limiter_flag"]::Int64
    main_store[56] = configData["core"]["vl_const"]::Float64
    main_store[57] = configData["core"]["aoa"]::Float64
    main_store[58] = configData["core"]["mach"]::Float64
    main_store[59] = configData["core"]["gamma"]::Float64
    main_store[60] = configData["core"]["pr_inf"]::Float64
    main_store[61] = configData["core"]["rho_inf"]::Float64
    main_store[62] = calculateTheta(configData)::Float64

    println(Int(getConfig()["core"]["max_iters"]) + 1)
    function run_code(ghost_holder, dist_keys, dist_globaldata, dist_q, dist_qpack, res_old, res_new, main_store)
        max_iters = Int(getConfig()["core"]["max_iters"])
        for i in 2:max_iters
            fpi_solver(i, ghost_holder, dist_keys, dist_globaldata, dist_q, dist_qpack, res_old, res_new, main_store)
        end
    end

    res_old = dzeros(nworkers())
    res_new = dzeros(nworkers())
    function test_code(ghost_holder, dist_keys, dist_globaldata, dist_q, dist_qpack, res_old, res_new, main_store)
        println("! Starting warmup function")
        fpi_solver(1, ghost_holder, dist_keys, dist_globaldata, dist_q, dist_qpack, res_old, res_new, main_store)
        println("! Starting main function")
        @timeit to "nest 1" begin
            run_code(ghost_holder, dist_keys, dist_globaldata, dist_q, dist_qpack, res_old, res_new, main_store)
        end
    end


    test_code(ghost_holder, dist_keys, dist_globaldata, dist_q, dist_qpack, res_old, res_new, main_store)
    println("! Work Completed")
    # # println(to)
    open("../results/timer" * string(numPoints) * "_" * string(getConfig()["core"]["max_iters"]) *
        "_" * string(length(workers())) *".txt", "w") do io
        print_timer(io, to)
    end

    # t = rmprocs(2, 3, waitfor=0)
    # wait(t)
    # rmprocs(length(workers()))
    # @spawnat 2 begin
    # println(dist_globaldata[1])
    # println(dist_globaldata[end])
    # globaldata = makelocal(dist_globaldata)
    # close(dist_globaldata)
    # end
    # compute_cl_cd_cm(globaldata, configData, shapeptsidx)

    # println(IOContext(stdout, :compact => false), globaldata[1].q)
    # println(IOContext(stdout, :compact => false), globaldata[1].dq)
    # println(IOContext(stdout, :compact => false), globaldata[100].q)
    # println(IOContext(stdout, :compact => false), globaldata[100].dq)
    # println(IOContext(stdout, :compact => false), globaldata[1000].q)
    # println(IOContext(stdout, :compact => false), globaldata[1000].dq)
    # println()
    # println(IOContext(stdout, :compact => false), globaldata[1].flux_res)
    # println(IOContext(stdout, :compact => false), globaldata[100].flux_res)
    # println(IOContext(stdout, :compact => false), globaldata[1000].flux_res)
    # println()
    # println(IOContext(stdout, :compact => false), globaldata[1].delta)
    # println(IOContext(stdout, :compact => false), globaldata[100].delta)
    # println(IOContext(stdout, :compact => false), globaldata[1000].delta)
    # println()
    # println(IOContext(stdout, :compact => false), dist_globaldata[1].prim)
    # println(IOContext(stdout, :compact => false), dist_globaldata[100].prim)
    # println(IOContext(stdout, :compact => false), globaldata[1000].prim)
    # println(IOContext(stdout, :compact => false), globaldata[100].ypos_conn)
    # println(IOContext(stdout, :compact => false), globaldata[100].yneg_conn)
    # println(globaldata[1])

    # file = open("../results/primvals" * string(numPoints) * ".txt", "w")
    # @showprogress 1 "This takes time" for (idx, _) in enumerate(dist_globaldata)
    #     primtowrite = dist_globaldata[global_local_direct_index[idx]].prim
    #     for element in primtowrite
    #         @printf(file,"%0.17f", element)
    #         @printf(file, " ")
    #     end
    #     print(file, "\n")
    # end
    # close(file)
    close(ghost_holder)
    close(dist_qpack)
    close(dist_q)
    close(dist_keys)
    close(dist_globaldata)
    close(res_old)
    close(res_new)
    d_closeall()
end
