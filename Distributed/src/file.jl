function returnFileLength(file_name::String)
    data1 = read(file_name, String)
    splitdata = split(data1, "\n")
    return length(splitdata) - 2
end

function createGlobalLocalMapIndex(global_local_map_index, global_local_direct_index, folder_name::String)
    index_flag = 1
    # println("Reading multiple files")
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
        # println(filename)
        data = read(filename, String)
        splitdata = @view split(data, "\n")[1:end-1]
        itmdata = split(splitdata[1])
        local_point_count = parse(Int,itmdata[3])
        for (idx, itm) in enumerate(splitdata)
            if idx == 1
                continue
            elseif idx <= local_point_count + 1
                itmdata = split(itm)
                global_local_direct_index[parse(Int,itmdata[1])] = index_flag
                global_local_map_index[(parse(Float64,itmdata[2]), parse(Float64, itmdata[3]))] = index_flag
                index_flag += 1
            end
        end
    end
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
        # println(filename)
        data = read(filename, String)
        splitdata = @view split(data, "\n")[1:end-1]
        ghost_holder[iter] = Dict{Int64,Point}()
        itmdata = split(splitdata[1])
        local_point_count = parse(Int,itmdata[3])
        ghost_point_count = parse(Int,itmdata[4])

        for (idx, itm) in enumerate(splitdata)
            if idx > local_point_count + 1
                itmdata = split(itm)
                ghost_holder[iter][idx-1] = dist_globaldata[global_local_map_index[(parse(Float64,itmdata[2]), parse(Float64, itmdata[3]))]]
            end
        end
    end
end

function readDistribuedFile(folder_name::String, defprimal, p, global_local_map_index)
    # println("Reading multiple files")
    # println(folder_name)
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
    itmdata = split(splitdata[1])
    local_point_count = parse(Int,itmdata[3])
    ghost_point_count = parse(Int,itmdata[4])
    local_points_holder = Array{Point,1}(undef, local_point_count)
    for (idx, itm) in enumerate(splitdata)
        if idx == 1
            continue
        elseif idx <= local_point_count + 1
            itmdata = split(itm)
            globalID = global_local_map_index[(parse(Float64,itmdata[2]), parse(Float64, itmdata[3]))]
            local_points_holder[idx-1] = Point(parse(Int,itmdata[1]),
                parse(Float64, itmdata[2]),
                parse(Float64, itmdata[3]),
                parse(Int, itmdata[4]),
                parse(Int, itmdata[5]),
                parse(Int8,itmdata[6]),
                parse(Int8,itmdata[7]),
                parse(Float64,itmdata[8]),
                parse(Int8,itmdata[9]),
                parse.(Int, itmdata[10:end]),
                0.0,
                0.0,
                copy(defprimal),
                zeros(Float64, 4),
                zeros(Float64, 4),
                zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4), 0.0, 0, 0, 0, 0, Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0),
                Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4),
                globalID)
        end
    end
    return local_points_holder
end

function readDistribuedFileQuadtree(folder_name::String, defprimal, p, global_local_map_index)
    # println(folder_name)
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
    # println(filename)
    data = read(filename, String)
    splitdata = @view split(data, "\n")[1:end-1]
    itmdata = split(splitdata[1])
    local_point_count = parse(Int,itmdata[3])
    ghost_point_count = parse(Int,itmdata[4])
    local_points_holder = Array{Point,1}(undef, local_point_count)
    for (idx, itm) in enumerate(splitdata)
        if idx == 1
            continue
        elseif idx <= local_point_count + 1
            itmdata = split(itm)
            globalID = global_local_map_index[(parse(Float64,itmdata[2]), parse(Float64, itmdata[3]))]
            local_points_holder[idx-1] = Point(parse(Int,itmdata[1]),
                parse(Float64, itmdata[2]),
                parse(Float64, itmdata[3]),
                parse(Int, itmdata[4]),
                parse(Int, itmdata[5]),
                parse(Int8,itmdata[6]),
                parse(Int8,itmdata[7]),
                parse(Float64,itmdata[11]),
                parse(Int8,itmdata[12]),
                parse.(Int, itmdata[13:end]),
                parse(Float64, itmdata[7]),
                parse(Float64, itmdata[8]),
                copy(defprimal),
                SVector{4}([zero(Float64) for iter in 1:4]),
                zeros(Float64, 4),
                SVector{4}([zero(Float64) for iter in 1:4]), 
                SVector{4}([zero(Float64) for iter in 1:4]), 
                SVector{4}([zero(Float64) for iter in 1:4]), 
                SVector{4}([zero(Float64) for iter in 1:4]), 
                0.0, 0, 0, 0, 0, Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0),
                Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4),
                globalID)
        end
    end
    return local_points_holder
end

function readDistribuedFileQ(folder_name::String, defprimal, p, global_local_map_index)

    # println(folder_name)
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
    # println(filename)
    data = read(filename, String)
    splitdata = @view split(data, "\n")[1:end-1]
    itmdata = split(splitdata[1])
    local_point_count = parse(Int,itmdata[3])

    local_points_holder = Array{TempQ,1}(undef, local_point_count)
    for idx in 1:local_point_count + 1
        if idx == 1
            continue
        elseif idx <= local_point_count + 1
            # itmdata = split(itm)
            local_points_holder[idx-1] = TempQ(SVector{4}([zero(Float64) for iter in 1:4]))
        end
    end

    return local_points_holder
end

function readDistribuedFileQPack(folder_name::String, defprimal, p, global_local_map_index)

    # println(folder_name)
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
    # println(filename)
    data = read(filename, String)
    splitdata = @view split(data, "\n")[1:end-1]
    itmdata = split(splitdata[1])
    local_point_count = parse(Int,itmdata[3])

    local_points_holder = Array{TempQPack,1}(undef, local_point_count)
    for idx in 1:local_point_count + 1
        if idx == 1
            continue
        elseif idx <= local_point_count + 1
            # itmdata = split(itm)
            local_points_holder[idx-1] = TempQPack(SVector{4}([zero(Float64) for iter in 1:4]), 
                                            SVector{4}([zero(Float64) for iter in 1:4]), 
                                            SVector{4}([zero(Float64) for iter in 1:4]), 
                                            SVector{4}([zero(Float64) for iter in 1:4]), 
                                            SVector{4}([zero(Float64) for iter in 1:4]))
        end
    end
    return local_points_holder
end

function returnKeys(loc_ghost_holder)
    keysIter = keys(loc_ghost_holder[1])
    keysLength = length(keysIter)
    keyHolder = zeros(Int64, keysLength)
    idx = 1
    for item in keysIter
        keyHolder[idx] = item
        idx += 1
    end
    return keyHolder
end