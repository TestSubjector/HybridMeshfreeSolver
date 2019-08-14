function cal_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, configData)
	phi_i = zeros(Float64,4)
	phi_k = zeros(Float64,4)
	G_i = zeros(Float64,4)
    G_k = zeros(Float64,4)
	result = zeros(Float64,4)
	qtilde_i = zeros(Float64,4)
	qtilde_k = zeros(Float64,4)
	dist_length = length(loc_globaldata)

	# updateLocalGhost(loc_ghost_holder, globaldata)

	for (idx, _) in enumerate(loc_globaldata)
		if loc_globaldata[idx].flag_1 == 0
			wallindices_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, dist_length, configData, idx, phi_i, phi_k)
		elseif loc_globaldata[idx].flag_1 == 2
			outerindices_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, dist_length, configData, idx, phi_i, phi_k)
		elseif loc_globaldata[idx].flag_1 == 1
			interiorindices_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, dist_length, configData, idx, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k)
		end
	end
	return nothing
end

function wallindices_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, dist_length, configData, itm, phi_i, phi_k)
	# for itm in wallindices
		# println(itm)
		Gxp = wall_dGx_pos(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k)
		Gxn = wall_dGx_neg(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k)
		Gyn = wall_dGy_neg(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k)
		GTemp = @.((Gxp + Gxn + Gyn) * 2)
		loc_globaldata[itm].flux_res = GTemp
		# if itm == 3
		# 	println(IOContext(stdout, :compact => false), Gxp)
		# 	println(IOContext(stdout, :compact => false), Gxp + Gxn)
		# 	println(IOContext(stdout, :compact => false), Gxp + Gxn + Gyn)
		# end
	# end
	return nothing
end

function outerindices_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, dist_length, configData, itm, phi_i, phi_k)
	# for itm in outerindices
		Gxp = outer_dGx_pos(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k)
		Gxn = outer_dGx_neg(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k)
		Gyp = outer_dGy_pos(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k)
		GTemp = @.(Gxp + Gxn + Gyp)

		loc_globaldata[itm].flux_res = GTemp
	# end
	return nothing
end

function interiorindices_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, dist_length, configData, itm, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k)
	# for itm in interiorindices
		Gxp = interior_dGx_pos(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k)
		Gxn = interior_dGx_neg(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k)
		Gyp = interior_dGy_pos(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k)
		Gyn = interior_dGy_neg(loc_globaldata, globaldata, loc_ghost_holder, dist_length, itm, configData, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k)
		# if itm == 1
		# 	println("=======")
		# 	println(IOContext(stdout, :compact => false), Gxp)
		# 	println(IOContext(stdout, :compact => false), Gxn + Gxp)
		# 	println(IOContext(stdout, :compact => false), Gxn + Gxp + Gyp)
		# 	println(IOContext(stdout, :compact => false), Gxn + Gxp + Gyp + Gyn)
		# 	println()
		# end
		GTemp = @.(Gxp + Gxn + Gyp + Gyn)
		loc_globaldata[itm].flux_res = GTemp
	# end
	return nothing
end
