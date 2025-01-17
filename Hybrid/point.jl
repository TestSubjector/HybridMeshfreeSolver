
# point = Point(locaslID, x, y, left, right, flag_1, flag_2, nbhs, conn, nx, ny, prim, flux_res, q, dq, entropy, xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs, xpos_conn, xneg_conn, ypos_conn, yneg_conn, delta)

function setNormals(self::Point, n)
    self.nx = n[1]
    self.ny = n[2]
    return nothing
end

function getxy(self)
    return (self.x, self.y)
end

function setConnectivity(self, conn)
    self.xpos_conn = conn[1]
    self.xpos_nbhs = length(conn[1])
    self.xneg_conn = conn[2]
    self.xneg_nbhs = length(conn[2])
    self.ypos_conn = conn[3]
    self.ypos_nbhs = length(conn[3])
    self.yneg_conn = conn[4]
    self.yneg_nbhs = length(conn[4])
    return nothing
end

# function setSmallestPointDistance(self::Point, distance)
#     self.short_distance = distance
# end

function convertToFixedArray(targetArray1, originalStruct::Point, idx)
    targetArray1[idx] = FixedPoint(originalStruct.localID,
                                        originalStruct.x,
                                        originalStruct.y,
                                        originalStruct.left,
                                        originalStruct.right,
                                        originalStruct.flag_1,
                                        originalStruct.flag_2,
                                        originalStruct.short_distance,
                                        originalStruct.nbhs,
                                        originalStruct.nx,
                                        originalStruct.ny,
                                        0.0
                                    )
    return nothing
end



function convertToNeighbourArray(targetArray2, originalStruct::Point, idx)
    targetArray2[1] = originalStruct.xpos_nbhs
    targetArray2[2] = originalStruct.xneg_nbhs
    targetArray2[3] = originalStruct.ypos_nbhs
    targetArray2[4] = originalStruct.yneg_nbhs
    targetArray2[5:4 + originalStruct.nbhs, idx] = originalStruct.conn
    targetArray2[15:14 + originalStruct.xpos_nbhs, idx] = originalStruct.xpos_conn
    targetArray2[25:24 + originalStruct.xneg_nbhs, idx] = originalStruct.xneg_conn
    targetArray2[35:34 + originalStruct.ypos_nbhs, idx] = originalStruct.ypos_conn
    targetArray2[45:44 + originalStruct.yneg_nbhs, idx] = originalStruct.yneg_conn
    return nothing
end