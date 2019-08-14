function returnFileLength(file_name::String)
    data1 = read(file_name, String)
    splitdata = split(data1, "\n")
    return length(splitdata) - 1
end

function createGlobalLocalMapIndex(global_local_map_index, global_local_direct_index, folder_name::String)
    index_flag = 1
    for iter in 1:length(workers())
        filename = folder_name * "/" * "partGrid000" * string(iter-1)
        data = read(filename, String)
        splitdata = @view split(data, "\n")[1:end-1]
        itmdata = split(splitdata[1], " ")
        local_point_count = parse(Int,itmdata[3])
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
end

# function readFile(file_name::String, globaldata, table, defprimal, wallptsidx, outerptsidx, Interiorptsidx, shapeptsidx,
#         wallpts, Interiorpts, outerpts, shapepts, numPoints)
#     data1 = read(file_name, String)
#     splitdata = @view split(data1, "\n")[1:end-1]
#     # print(splitdata[1:3])
#     @showprogress 1 "Computing ReadFile" for (idx, itm) in enumerate(splitdata)
#         itmdata = split(itm, " ")
#         globaldata[idx] = Point(idx,
#                     parse(Float64,itmdata[1]),
#                     parse(Float64, itmdata[2]),
#                     parse(Int, itmdata[3]),
#                     parse(Int, itmdata[4]),
#                     parse(Int8,itmdata[5]),
#                     parse(Int8,itmdata[6]),
#                     parse(Float64,itmdata[7]),
#                     parse(Int8,itmdata[8]),
#                     parse.(Int, itmdata[9:end-1]),
#                     0.0,
#                     0.0,
#                     copy(defprimal),
#                     zeros(Float64, 4),
#                     zeros(Float64, 4),
#                     Array{Array{Float64,1},1}(undef, 2), 0.0, 0, 0, 0, 0, Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0),
#                     Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4))
#     end
#     return nothing
# end

function readGhostFile(folder_name::String, ghost_holder, global_local_map_index, dist_globaldata)
    # println(ghost_folder_name)
    for iter in 1:length(workers())
        filename = folder_name * "/" * "partGrid000" * string(iter-1)
        println(filename)
        data = read(filename, String)
        splitdata = @view split(data, "\n")[1:end-1]
        ghost_holder[iter] = Dict{Int64,Point}()
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

function readDistribuedFile(folder_name::String, local_points_holder, defprimal, p, global_local_map_index)
    println("Reading multiple files")
    println(folder_name)
    for iter in p-1:p-1
        filename = folder_name * "/" * "partGrid000" * string(iter-1)
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
                if iter == 1
                    local_points_holder[idx-1] = Point(parse(Int,itmdata[1]),
                        parse(Float64,itmdata[2]),
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
                        Array{Array{Float64,1},1}(undef, 2), 0.0, 0, 0, 0, 0, Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0),
                        Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4),
                        globalID)
                else
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
                        zeros(Float64, 4),
                        Array{Array{Float64,1},1}(undef, 2), 0.0, 0, 0, 0, 0, Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0),
                        Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4),
                        globalID)
                end
            end
        end
    end
    return local_points_holder
end

# function readBasicInfo