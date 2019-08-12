function main()

    configData = getConfig()
    wallpts, Interiorpts, outerpts, shapepts = 0,0,0,0

    file_name = string(ARGS[1])
    ghost_folder_name = string(ARGS[2])

    numPoints = returnFileLength(file_name)
    println(numPoints)
    globaldata = Array{Point,1}(undef, numPoints)
    table = Array{Int32,1}(undef, numPoints)
    wallptsidx = Array{Int32,1}(undef, 0)
    Interiorptsidx = Array{Int32,1}(undef, 0)
    outerptsidx = Array{Int32,1}(undef, 0)
    shapeptsidx = Array{Int32,1}(undef, 0)
    # print(splitdata[1:3])
    defprimal = getInitialPrimitive(configData)

    local_points_holder = Array{Array{Point,1},1}(undef, nworkers())
    ghost_holder = Array{Dict{Int64,Point},1}(undef, nworkers())

    println("Start Read")


    # count = 0
    # readFile(file_name::String, globaldata, table, defprimal, wallptsidx, outerptsidx, Interiorptsidx, shapeptsidx,
        # wallpts, Interiorpts, outerpts, shapepts, numPoints)

    # format = configData["format"]["type"]
    # if format == 1
    #     interior = configData["point"]["interior"]
    #     wall = configData["point"]["wall"]
    #     outer = configData["point"]["outer"]
    #     @showprogress 2 "Computing Connectivity" for idx in 1:numPoints
    #         placeNormals(globaldata, idx, configData, interior, wall, outer)
    #     end
    # end

    # println("Start table sorting")
    # @showprogress 3 "Computing Table" for idx in table
    #     connectivity = calculateConnectivity(globaldata, idx)
    #     setConnectivity(globaldata[idx], connectivity)
    #     # smallest_dist(globaldata, idx)
    #     # if idx % (length(table) * 0.25) == 0
    #     #     println("Bump In Table")
    #     # end
    # end

    ras = [@spawnat p readDistribuedFile(ghost_folder_name::String, local_points_holder, defprimal, p) for p in workers()]
    rasa = reshape(ras, (nworkers()))
    points_holder = DArray(rasa)
    # readGhostFile(ghost_folder_name, ghost_holder, globaldata)
    # ghost_holder = distribute(ghost_holder, procs=workers(), dist=(length(workers()),))
    # lph = distribute(local_points_holder, procs=workers(), dist=(nworkers(),))
    # # @sync for pid in procs(local_points_holder)
    # #   @spawnat pid begin
    #     #   println(localindices(local_points_holder))
    # #   end
    # # end

    # points_holder = DArray(reshape([lph[1]..., lph[2]...,lph[3]...,lph[4]...], :))
    @sync for pid in procs(points_holder)
        @spawnat pid begin
            println(localindices(points_holder))
        end
    end

    exit()

    # println(wallptsidx)

    # println(globaldata[3])
    # return

    # configData = distribute(configData)
    # println(globaldata[1])
    # println(globaldata[end])
    # println("Size is before: ",length(globaldata))
    dist_globaldata = distribute(globaldata, procs=workers(), dist=(length(workers()),))
    # println("Size is: ", length(globaldata))

    println(Int(getConfig()["core"]["max_iters"]) + 1)
    function run_code(ghost_holder, dist_globaldata, configData, wallptsidx::Array{Int32,1}, outerptsidx::Array{Int32,1}, Interiorptsidx::Array{Int32,1}, res_old, numPoints)
        for i in 1:(Int(getConfig()["core"]["max_iters"]))
            fpi_solver(i, ghost_holder, dist_globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old, numPoints)
        end
    end

    res_old = zeros(Float64, 1)
    function test_code(ghost_holder, dist_globaldata, configData, wallptsidx::Array{Int32,1}, outerptsidx::Array{Int32,1}, Interiorptsidx::Array{Int32,1}, res_old, numPoints)
        println("! Starting warmup function")
        fpi_solver(1, ghost_holder, dist_globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old, numPoints)
        res_old[1] = 0.0
        # Profile.clear_malloc_data()
        # @trace(fpi_solver(1, globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old), maxdepth = 3)
        # res_old[1] = 0.0
        # fpi_solver(1, globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old)
        # @profile fpi_solver(1, globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old)
        # Profile.print()
        # res_old[1] = 0.0
        println("! Starting main function")
        @timeit to "nest 4" begin
            run_code(ghost_holder, dist_globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old, numPoints)
        end
    end


    test_code(ghost_holder, dist_globaldata, configData, wallptsidx, outerptsidx, Interiorptsidx, res_old, numPoints)
    println("! Work Completed")
    # # println(to)
    open("results/timer" * string(numPoints) * "_" * string(getConfig()["core"]["max_iters"]) *
        "_" * string(length(workers())) *".txt", "w") do io
        print_timer(io, to)
    end

    # t = rmprocs(2, 3, waitfor=0)
    # wait(t)
    # rmprocs(length(workers()))
    # @spawnat 2 begin
    println(dist_globaldata[1])
    println(dist_globaldata[end])
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

    file = open("results/primvals" * string(numPoints) * ".txt", "w")
    for (idx, _) in enumerate(dist_globaldata)
        primtowrite = dist_globaldata[idx].prim
        for element in primtowrite
            @printf(file,"%0.17f", element)
            @printf(file, " ")
        end
        print(file, "\n")
    end
    close(file)
    close(ghost_holder)
    close(dist_globaldata)
end
