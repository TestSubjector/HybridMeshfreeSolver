function main()

    configData = getConfig()
    wallpts, Interiorpts, outerpts, shapepts = 0,0,0,0

    file_name = string(ARGS[2])
    folder_name = string(ARGS[3])

    numPoints = returnFileLength(file_name)
    println(numPoints)
    # globaldata = Array{Point,1}(undef, numPoints)
    # table = Array{Int32,1}(undef, numPoints)
    # wallptsidx = Array{Int32,1}(undef, 0)
    # Interiorptsidx = Array{Int32,1}(undef, 0)
    # outerptsidx = Array{Int32,1}(undef, 0)
    # shapeptsidx = Array{Int32,1}(undef, 0)
    # print(splitdata[1:3])
    defprimal = getInitialPrimitive(configData)

    # local_points_holder = Array{Array{Point,1},1}(undef, nworkers())
    # res_old_holder = Array{Float64,1}(undef, nworkers())
    ghost_holder = Array{Dict{Int32,Point},1}(undef, nworkers())
    ghost_holder_mutable = Array{Dict{Int32, Array{Float64,1}}, 1}(undef, nworkers())
    # q_ghost_holder = Array{Dict{Int64,TempQ},1}(undef, nworkers())
    # dq_ghost_holder = Array{Dict{Int64,TempDQ},1}(undef, nworkers())
    # prim_ghost_holder = Array{Dict{Int64,TempPrim},1}(undef, nworkers())

    global_local_map_index = Dict{Tuple{Float64, Float64},Int64}()
    global_local_direct_index = Array{Int32, 1}(undef, numPoints)
    files_length = Array{Int64, 1}(undef, numPoints)
    println("Indexing")

    createGlobalLocalMapIndex(global_local_map_index, global_local_direct_index, folder_name::String, files_length)
    # println(global_local_map_index)

    println("Start Read")
    # count = 0
    # readFile(file_name::String, globaldata, table, defprimal, wallptsidx, outerptsidx, Interiorptsidx, shapeptsidx,
        # wallpts, Interiorpts, outerpts, shapepts, numPoints)

    println("Read Local Files")
    # readDistribuedFile(folder_name::String, defprimal, 12, global_local_map_index)
    globaldata_parts = [@spawnat p readDistribuedFile(folder_name::String, defprimal, p, global_local_map_index) for p in workers()]
    globaldata_parts = reshape(globaldata_parts, (nworkers()))
    dist_globaldata = DArray(globaldata_parts)
    # globaldata_parts = nothing

    q_parts = [@spawnat p readDistribuedFileQ(folder_name::String, defprimal, p, global_local_map_index, files_length) for p in workers()]
    q_parts = reshape(q_parts, (nworkers()))
    dist_q = DArray(q_parts)
    # q_parts = nothing

    dq_parts = [@spawnat p readDistribuedFileDQ(folder_name::String, defprimal, p, global_local_map_index, files_length) for p in workers()]
    dq_parts = reshape(dq_parts, (nworkers()))
    dist_dq = DArray(dq_parts)
    # dq_parts = nothing

    globaldata_parts_mutable = [@spawnat p readDistribuedFileMutables(folder_name::String, p, files_length) for p in workers()]
    globaldata_parts_mutable = reshape(globaldata_parts_mutable, (1,(nworkers())))
    dist_globaldata_mutable = DArray(globaldata_parts_mutable)
    # globaldata_parts_mutable = nothing

    println("Reading Ghost")
    readGhostFile(folder_name, ghost_holder, global_local_map_index, dist_globaldata)
    ghost_holder = distribute(ghost_holder, procs=workers(), dist=(length(workers()),))

    readGhostFileMutables(folder_name, ghost_holder_mutable)
    ghost_holder_mutable = distribute(ghost_holder_mutable, procs=workers(), dist=(length(workers()),))

    format = configData["format"]["type"]
    if format == 1
        interior = configData["point"]["interior"]
        wall = configData["point"]["wall"]
        outer = configData["point"]["outer"]
        @sync for pid in procs(dist_globaldata)
            @spawnat pid begin
                placeNormals(dist_globaldata[:L], dist_globaldata, ghost_holder[:L], interior, wall, outer)
            end
        end
    end

    println("Start table sorting")
    @sync for pid in procs(dist_globaldata)
        @spawnat pid begin
            calculateConnectivity(dist_globaldata[:L], dist_globaldata, ghost_holder[:L])
            # setConnectivity(dist_globaldata[idx], connectivity)
        end
    end
    # @showprogress 3 "Computing Table" for idx in table
    #     connectivity = calculateConnectivity(globaldata, idx)
    #     setConnectivity(globaldata[idx], connectivity)
    #     # smallest_dist(globaldata, idx)
    #     # if idx % (length(table) * 0.25) == 0
    #     #     println("Bump In Table")
    #     # end
    # end

    # points_holder = DArray(reshape([lph[1]..., lph[2]...,lph[3]...,lph[4]...], :))
    @sync for pid in procs(dist_globaldata)
        @spawnat pid begin
            println(localindices(dist_globaldata))
        end
    end

    println("! Transferring FixedData to GPU")
    @sync for pid in workers()
        @spawnat pid begin
            numPoints = files_length[pid-1]
            global gpuLocNumPoints = numPoints
            locDataFixedPoint = Array{FixedPoint,1}(undef, numPoints)
            locDataConn = zeros(Int32, 55, numPoints)
            locGlobalData = dist_globaldata[:L]
            for idx in 1: numPoints
                convertToFixedArray(locDataFixedPoint, locGlobalData[idx], idx)
                convertToNeighbourArray(locDataConn, locGlobalData[idx], idx)
            end

            # locGhostGlobalData = ghost_holder[:L][1]
            # localkeys = keys(locGhostGlobalData)
            # locGhostDataFixedPoint = Array{Dict{Int32,FixedPoint},1}(undef, 1)
            # locGhostDataFixedPoint[1] = Dict{Int32,FixedPoint}()
            # for iter in localkeys
            #     convertToFixedArray(locGhostDataFixedPoint[1], locGhostGlobalData[iter], iter)
            # end

            global gpuLocDataFixedPoint = CuArray(locDataFixedPoint)
            global gpuLocDataConn = CuArray(locDataConn)
            # global gpuLocGhostDataFixedPoint = CuArray(locGhostDataFixedPoint)
            # global gpuLocGhostGlobalDataMutable = CuArray(locGhostGlobalDataMutable)
            # @cuda changeToOne(cutest)
            # part_test[:L] = Array(cutest)

            global gpuConfigData = CuArray([
                            getConfig()["core"]["points"],#1
                            getConfig()["core"]["cfl"],
                            getConfig()["core"]["max_iters"],
                            getConfig()["core"]["mach"],
                            getConfig()["core"]["aoa"],#5
                            getConfig()["core"]["power"],
                            getConfig()["core"]["limiter_flag"],
                            getConfig()["core"]["vl_const"],
                            getConfig()["core"]["initial_conditions_flag"],
                            getConfig()["core"]["interior_points_normal_flag"],#10
                            getConfig()["core"]["shapes"],
                            getConfig()["core"]["rho_inf"],
                            getConfig()["core"]["pr_inf"],
                            getConfig()["core"]["threadsperblock"],
                            getConfig()["core"]["gamma"],#15
                            getConfig()["core"]["clcd_flag"],
                            getConfig()["point"]["wall"],
                            getConfig()["point"]["interior"],
                            getConfig()["point"]["outer"]
                        ])
        end
    end

    # println(dist_globaldata[3])

    # println(wallptsidx)


    # return

    # configData = distribute(configData)
    # println(globaldata[1])
    # println(globaldata[end])
    # println("Size is before: ",length(globaldata))
    # dist_globaldata = distribute(globaldata, procs=workers(), dist=(length(workers()),))
    # println("Size is: ", length(globaldata))

    println(Int(getConfig()["core"]["max_iters"]) + 1)
    function run_code(ghost_holder, dist_globaldata, dist_q, dist_dq, dist_globaldata_mutable, ghost_holder_mutable, configData, res_old, numPoints)
        fpi_solver((Int(getConfig()["core"]["max_iters"])), ghost_holder, dist_globaldata, dist_q, dist_dq, dist_globaldata_mutable, ghost_holder_mutable, configData, res_old, numPoints)
    end

    res_old = zeros(Float64, 1)
    function test_code(ghost_holder, dist_globaldata, dist_q, dist_dq, dist_globaldata_mutable, ghost_holder_mutable, configData, res_old, numPoints)
        println("! Starting warmup function")
        fpi_solver(1, ghost_holder, dist_globaldata, dist_q, dist_dq, dist_globaldata_mutable, ghost_holder_mutable, configData, res_old, numPoints)
        res_old[1] = 0.0
        # Profile.clear_malloc_data()
        # @trace(fpi_solver(1, globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old), maxdepth = 3)
        # res_old[1] = 0.0
        # fpi_solver(1, globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old)
        # @profile fpi_solver(1, globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old)
        # Profile.print()
        # res_old[1] = 0.0
        println("! Starting main function")
        @timeit to "nest 1" begin
            run_code(ghost_holder, dist_globaldata, dist_q, dist_dq, dist_globaldata_mutable, ghost_holder_mutable, configData, res_old, numPoints)
        end
    end

    test_code(ghost_holder, dist_globaldata, dist_q, dist_dq, dist_globaldata_mutable, ghost_holder_mutable, configData, res_old, numPoints)
    println("! Work Completed")
    # # println(to)
    open("temp/timer" * string(numPoints) * "_" * string(getConfig()["core"]["max_iters"]) *
        "_" * string(length(workers())) *".txt", "w") do io
        print_timer(io, to)
    end

    println("! Free GPU")
    @sync for pid in workers()
        @spawnat pid begin
            Array(gpuLocDataFixedPoint)
            Array(gpuLocDataConn)
            Array(gpuLocGhostDataFixedPoint)
            # Array(gpuLocGhostGlobalDataMutable)
            # @cuda changeToOne(cutest)
            # part_test[:L] = Array(cutest)
        end
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
    # println(IOContext(stdout, :compact => false), globaldata[1].prim)
    # println(IOContext(stdout, :compact => false), globaldata[100].prim)
    # println(IOContext(stdout, :compact => false), globaldata[1000].prim)
    # println(IOContext(stdout, :compact => false), globaldata[100].ypos_conn)
    # println(IOContext(stdout, :compact => false), globaldata[100].yneg_conn)
    # println(globaldata[1])

    # file = open("results/primvals"* "_" * string(getConfig()["core"]["max_iters"]) *
    #  "_" * string(numPoints) * ".txt", "w")
    # @showprogress 1 "This takes time" for (idx, _) in enumerate(dist_globaldata)
    #     primtowrite = dist_globaldata[global_local_direct_index[idx]].prim
    #     for element in primtowrite
    #         @printf(file,"%0.17f", element)
    #         @printf(file, " ")
    #     end
    #     print(file, "\n")
    # end
    # close(file)

    println("! Write To Files")
    @sync for ip in procs(dist_globaldata)
        @spawnat ip begin
            writeToFile(dist_globaldata[:L], numPoints)
        end
    end

    println("! Close Distributed Arrays")
    # d_closeall()
    exit()
    close(ghost_holder)
    close(ghost_holder_mutable)
    close(dist_dq)
    close(dist_q)
    close(dist_globaldata)
    close(dist_globaldata_mutable)
end

function testraid()
    test = [@spawnat p returnArrayOne() for p in workers()]
    test_parts = reshape(test, (nworkers()))
    part_test = DArray(test_parts)
    println(part_test)
    @timeit to "nest 4" begin
        @sync for pid in workers()
            @spawnat pid begin
                global cutest = CuArray(part_test[:L])
                # @cuda changeToOne(cutest)
                # part_test[:L] = Array(cutest)
            end
        end

        @sync for pid in workers()
            @spawnat pid begin
                # println(part_test[:L])
                @cuda changeToOne(cutest)
                # part_test[:L] = Array(cutest)
            end
        end

        @sync for pid in workers()
            @spawnat pid begin
                # println(part_test[:L])
                # @cuda changeToOne(cutest)
                part_test[:L] = Array(cutest)
            end
        end
    end


    println(part_test)
    open("temp/timer_stop.txt", "w") do io
        print_timer(io, to)
    end
    close(part_test)
end
