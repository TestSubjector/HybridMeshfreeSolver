function returnFileLength(file_name::String)
    data1 = read(file_name, String)
    splitdata = split(data1, "\n")
    return length(splitdata) - 1
end

function readFile(file_name::String, globaldata, table, defprimal, wallptsidx, outerptsidx, Interiorptsidx, shapeptsidx,
        wallpts, Interiorpts, outerpts, shapepts, numPoints)
    data1 = read(file_name, String)
    splitdata = @view split(data1, "\n")[1:end-1]
    # print(splitdata[1:3])
    @showprogress 1 "Computing ReadFile" for (idx, itm) in enumerate(splitdata)
        itmdata = split(itm, " ")
        globaldata[idx] = Point(idx,
                    parse(Float64,itmdata[1]),
                    parse(Float64, itmdata[2]),
                    parse(Int, itmdata[3]),
                    parse(Int, itmdata[4]),
                    parse(Int8,itmdata[5]),
                    parse(Int8,itmdata[6]),
                    parse(Float64,itmdata[7]),
                    parse(Int8,itmdata[8]),
                    parse.(Int, itmdata[9:end-1]),
                    0.0,
                    0.0,
                    copy(defprimal),
                    zeros(Float64, 4),
                    zeros(Float64, 4),
                    Array{Array{Float64,1},1}(undef, 2), 0.0, 0, 0, 0, 0, Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0),
                    Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4))
    end
    return nothing
end

function readGhostFile(ghost_folder_name::String, ghost_holder, globaldata)
    println(ghost_folder_name)
    for iter in 1:length(workers())
        filename = ghost_folder_name * "/" * "partGrid000" * string(iter-1)
        println(filename)
        data = read(filename, String)
        splitdata = @view split(data, "\n")[1:end-1]
        ghost_holder[iter] = Dict{Int64,Point}()
        for (idx, itm) in enumerate(splitdata)
            # print(idx)
            ghost_holder[iter][parse(Float64,itm)] = globaldata[idx]
        end
    end
end

function readDistribuedFile(folder_name::String, local_points_holder, defprimal, p)
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
        @showprogress 1 "Computing ReadFile" for (idx, itm) in enumerate(splitdata)
            if idx == 1
                continue
            elseif idx <= local_point_count + 1
                itmdata = split(itm, " ")
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
                        Array{Int32,1}(undef, 0), Array{Int32,1}(undef, 0), 0.0, zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4))
            end
        end
    end
    return local_points_holder
end

# function readBasicInfo