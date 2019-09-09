function returnFileLength(file_name::String)
    data1 = read(file_name, String)
    splitdata = split(data1, "\n")
    return length(splitdata) - 1
end

function createGlobalLocalMapIndex(global_local_map_index, global_local_direct_index, folder_name::String, files_length, actual_files_length)
    index_flag = 1
    for iter in 1:length(workers())
        if iter - 1 < 10
            filename = folder_name * "/" * "partGrid000" * string(iter-1)
        elseif iter - 1 < 100
            filename = folder_name * "/" * "partGrid00" * string(iter-1)
        elseif iter - 1 < 1000
            filename = folder_name * "/" * "partGrid0" * string(iter-1)
        else
            filename = folder_name * "/" * "partGrid" * string(iter-1)
        end
        data = read(filename, String)
        splitdata = @view split(data, "\n")[1:end-1]
        itmdata = split(splitdata[1], " ")
        local_point_count = parse(Int,itmdata[3])
        ghost_point_count = parse(Int,itmdata[4])
        files_length[iter] = local_point_count
        actual_files_length[iter] = local_point_count + ghost_point_count
        for (idx, itm) in enumerate(splitdata)
            if idx == 1
                continue
            elseif idx <= local_point_count + 1
                itmdata = split(itm, " ")
                global_local_direct_index[parse(Int,itmdata[1])] = index_flag
                global_local_map_index[(parse(Float64,itmdata[2]), parse(Float64, itmdata[3]))] = index_flag
                index_flag += 1
            end
        end
    end
    return nothing
end

function readDistribuedFile(folder_name::String, defprimal, p, global_local_map_index)
    println("Reading multiple files")
    println(folder_name)
    iter = p - 1
    if iter - 1 < 10
        filename = folder_name * "/" * "partGrid000" * string(iter-1)
    elseif iter - 1 < 100
        filename = folder_name * "/" * "partGrid00" * string(iter-1)
    elseif iter - 1 < 1000
        filename = folder_name * "/" * "partGrid0" * string(iter-1)
    else
        filename = folder_name * "/" * "partGrid" * string(iter-1)
    end
    println(filename)
    data = read(filename, String)
    splitdata = @view split(data, "\n")[1:end-1]
    itmdata = split(splitdata[1], " ")
    local_point_count = parse(Int,itmdata[3])
    ghost_point_count = parse(Int,itmdata[4])
    local_points_holder = Array{Point,1}(undef, local_point_count)
    for (idx, itm) in enumerate(splitdata)
        if idx == 1
            continue
        elseif idx <= local_point_count + 1
            itmdata = split(itm, " ")
            globalID = global_local_map_index[(parse(Float64,itmdata[2]), parse(Float64, itmdata[3]))]
            local_points_holder[idx-1] = Point(parse(Int,itmdata[1]),
                parse(Float64,itmdata[2]),
                parse(Float64, itmdata[3]),
                parse(Int, itmdata[4]),
                parse(Int, itmdata[5]),
                parse(Int8,itmdata[6]),
                parse(Int8,itmdata[7]),
                parse(Float64,itmdata[8]),
                parse(Int8,itmdata[9]),
                parse.(Int, itmdata[10:end-1]),
                0.0,
                0.0,
                copy(defprimal),
                zeros(Float64, 4),
                Array{Array{Float64,1},1}(undef, 2), 0.0, 0, 0, 0, 0, Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0),
                Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4),
                globalID)
        end
    end
    return local_points_holder
end

function readDistribuedFileQ(folder_name::String, defprimal, p, global_local_map_index, files_length)
    println("Setting Up Q")
    # println(folder_name)
    iter = p - 1
    local_points_count = files_length[iter]

    local_points_holder = Array{TempQ,1}(undef, local_points_count)
    for idx in 1:local_points_count
        local_points_holder[idx] = TempQ(zeros(Float64, 4))
    end
    return local_points_holder
end

function readDistribuedFileDQ(folder_name::String, defprimal, p, global_local_map_index, files_length)
    println("Setting Up DQ")
    # println(folder_name)
    iter = p - 1
    local_points_count = files_length[iter]

    local_points_holder = Array{TempDQ,1}(undef, local_points_count)

    for idx in 1:local_points_count
        local_points_holder[idx] = TempDQ(Array{Array{Float64,1},1}(undef, 2))
    end

    return local_points_holder
end

function readGhostFile(folder_name::String, ghost_holder, global_local_map_index, dist_globaldata)
    # println(ghost_folder_name)
    for iter in 1:length(workers())
        if iter - 1 < 10
            filename = folder_name * "/" * "partGrid000" * string(iter-1)
        elseif iter - 1 < 100
            filename = folder_name * "/" * "partGrid00" * string(iter-1)
        elseif iter - 1 < 1000
            filename = folder_name * "/" * "partGrid0" * string(iter-1)
        else
            filename = folder_name * "/" * "partGrid" * string(iter-1)
        end
        println(filename)
        data = read(filename, String)
        splitdata = @view split(data, "\n")[1:end-1]
        ghost_holder[iter] = Dict{Int32,Point}()
        itmdata = split(splitdata[1], " ")
        local_point_count = parse(Int,itmdata[3])
        ghost_point_count = parse(Int,itmdata[4])

        for (idx, itm) in enumerate(splitdata)
            if idx > local_point_count + 1
                itmdata = split(itm, " ")
                ghost_holder[iter][idx-1] = dist_globaldata[global_local_map_index[(parse(Float64,itmdata[2]), parse(Float64, itmdata[3]))]]
            end
        end
    end
end

function readDistribuedFileMutables(folder_name::String, p, actual_files_length)
    println("Setting Up Mutables")
    # println(folder_name)
    iter = p - 1
    local_points_count = actual_files_length[iter]

    local_points_holder = zeros(Float64, 28, local_points_count)
    return local_points_holder
end

function writeToFile(loc_globaldata, numPoints)
    file = open("results/split/primvalscuda"* "_" * string(getConfig()["core"]["max_iters"]) *
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

# function readGhostFileMutables(folder_name::String, ghost_holder_mutable)
#     # println(ghost_folder_name)
#     for iter in 1:length(workers())
#         if iter - 1 < 10
#             filename = folder_name * "/" * "partGrid000" * string(iter-1)
#         elseif iter - 1 < 100
#             filename = folder_name * "/" * "partGrid00" * string(iter-1)
#         elseif iter - 1 < 1000
#             filename = folder_name * "/" * "partGrid0" * string(iter-1)
#         else
#             filename = folder_name * "/" * "partGrid" * string(iter-1)
#         end
#         println(filename)
#         data = read(filename, String)
#         splitdata = @view split(data, "\n")[1:end-1]
#         ghost_holder_mutable[iter] = Dict{Int32,Array{Float64,1}}()
#         itmdata = split(splitdata[1], " ")
#         local_point_count = parse(Int,itmdata[3])
#         ghost_point_count = parse(Int,itmdata[4])
#         for (idx, itm) in enumerate(splitdata)
#             if idx > local_point_count + 1
#                 ghost_holder_mutable[iter][idx-1] = zeros(Float64, 20)
#             end
#         end
#     end
# end