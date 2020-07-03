function cal_flux_residual(loc_globaldata, globaldata, loc_ghost_holder, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k,
    result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, main_store)

    dist_length = length(loc_globaldata)
    power = main_store[53]
    limiter_flag = main_store[55]
    vl_const = main_store[56]
    gamma = main_store[59]

	for idx in 1:dist_length
		if loc_globaldata[idx].flag_1 == 0
			wallindices_flux_residual(loc_globaldata, loc_ghost_holder, dist_length, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const )
		elseif loc_globaldata[idx].flag_1 == 2
			outerindices_flux_residual(loc_globaldata, loc_ghost_holder, dist_length, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const )
		elseif loc_globaldata[idx].flag_1 == 1
			interiorindices_flux_residual(loc_globaldata, loc_ghost_holder, dist_length, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const)
		end
	end
	return nothing
end

function wallindices_flux_residual(loc_globaldata, loc_ghost_holder, dist_length, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const )
		wall_dGx_pos(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gxp )
		wall_dGx_neg(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gxn )
		wall_dGy_neg(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gyn )
		loc_globaldata[idx].flux_res = SVector{4}((Gxp + Gxn + Gyn) * 2)

	return nothing
end

function outerindices_flux_residual(loc_globaldata, loc_ghost_holder, dist_length, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const )
	outer_dGx_pos(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gxp )
	outer_dGx_neg(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gxn )
	outer_dGy_pos(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gyp )
	loc_globaldata[idx].flux_res = SVector{4}(Gxp + Gxn + Gyp)
	return nothing
end

function interiorindices_flux_residual(loc_globaldata, loc_ghost_holder, dist_length, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const)
	interior_dGx_pos(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gxp)
	interior_dGx_neg(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gxn)
	interior_dGy_pos(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gyp)
	interior_dGy_neg(loc_globaldata, loc_ghost_holder, dist_length, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, ∑_Δx_Δf, ∑_Δy_Δf, power, limiter_flag, vl_const, Gyn)
	loc_globaldata[idx].flux_res = SVector{4}(Gxp + Gxn + Gyp + Gyn)
	return nothing
end
