include("read_func.jl")

#########################
### Scaling dimension ###
#########################
###
###  FUNCTION: Scaling Dimension (Not Susceptibility?)
###
@everywhere function h_lnpp_inv(conf::Conf_type, spt::Spt_type, i_field::Int, i::Int)
    local const N = size(conf.phi, 1)
    local Q = div(N, 2)
    local const L = Int(floor(sqrt(N)))
    local p = p_matrix2(N) * (L-1)^2
    local num_f = size(conf.phi, 2)
    local J = zeros(Complex{Float64}, 2N*num_f, 2N*num_f)
    for it = 1:num_f
        local n = (it-1) * 2N
        J[(1+n):(N+n),    (1+n):(N+n)]    = p
        J[(N+1+n):(2N+n), (N+1+n):(2N+n)] = - conj(p)
    end
    local tmp = zeros(Complex{Float64}, N)#
    local rs = zeros(Complex{Float64}, N)#
    for j = 1:size(conf.phi, 3)
        if (norm(conf.phi[:, i_field, j, i]) == 0)
            continue
        else
            local wall = superpotential_2(conf.phi[:, :, j, i], spt)
            local col = 1
            for i1 = 1:num_f
                local r = (i1-1) * 2N
                for i2 = i1:num_f
                    local c = (i2-1) * 2N
                    local w2 = fieldtomatrix(wall[:,col])
                    J[(N+1+r):(2N+r), (1+c):(N+c)] = w2
                    J[(1+r):(N+r), (N+1+c):(2N+c)] = w2'
                    if i1!=i2
                        J[(N+1+c):(2N+c), (1+r):(N+r)] = w2
                        J[(1+c):(N+c), (N+1+r):(2N+r)] = w2'
                    end
                    col += 1
                end
            end
            local Jinv = inv(J)
            local n_field = (i_field-1) * 2N
            local u = rs + diag(Jinv)[(1+n_field):(N+n_field)] * conf.sgn[j, i]#
            local s = tmp + u#
            local t = s - tmp#
            rs = u - t#
            tmp = s#
        end
    end
    return tmp
end
#  \langle \psi_1(p) \Bar\psi_\Dot{1}(-p)
#  with Standard Deviation
function h_lnpp(conf::Conf_type, spt::Spt_type, i_field::Int)
    local const N = size(conf.phi, 1)
    local Q = div(N, 2)
    local const L = Int(floor(sqrt(N)))
    local pp = zeros(Complex{Float64}, N, size(conf.phi, 4))
    pp = @parallel hcat for i = 1:size(conf.phi, 4)
        h_lnpp_inv(conf, spt, i_field, i)
    end
    local pp_sum, pp_dev
    pp_sum, pp_dev = wt_dev(pp)
    local Delta, Del_dev
    Delta, Del_dev = witten_ind(conf.sgn)
    local pm = abs(p_matrix(N))
    local res = zeros(Float64, N-1, 3)
    res[1:Q, 1]       = log(abs2(pm)[1:Q])
    res[(Q+1):(2Q),1] = log(abs2(pm)[(Q+2):N])
    res[1:Q, 2]       = log((L-1)^4 * abs(pp_sum[1:Q])     ./ pm[1:Q]     / Delta)
    res[(Q+1):(2Q),2] = log((L-1)^4 * abs(pp_sum[(Q+2):N]) ./ pm[(Q+2):N] / Delta)
    local tmp = (real(pp_sum).*real(pp_dev)).^2 + (imag(pp_sum).*imag(pp_dev)).^2
    res[1:Q, 3]       = sqrt(tmp[1:Q]    ./(abs2(pp_sum[1:Q]    )).^2 + (Del_dev/Delta)^2)
    res[(Q+1):(2Q),3] = sqrt(tmp[(Q+2):N]./(abs2(pp_sum[(Q+2):N])).^2 + (Del_dev/Delta)^2)
    return res
end
#function h_lnpp_old(conf::Conf_type, spt::Spt_type, i_field::Int)
#    local const N = size(conf.phi, 1)
#    local Q = div(N, 2)
#    local const L = Int(floor(sqrt(N)))
#    local p = p_matrix2(N) * (L-1)^2
#    local J = zeros(Complex{Float64}, 2N, 2N)
#    J[1:N, 1:N] = p
#    J[(N+1):2N, (N+1):2N] = - conj(p)
#    local pp = zeros(Complex{Float64}, N, size(conf.phi, 4))
#    for i = 1:size(conf.phi, 4)
#        local tmp = zeros(Complex{Float64}, N)#
#        local rs = zeros(Complex{Float64}, N)#
#        for j = 1:size(conf.phi, 3)
#            if (norm(conf.phi[:, i_field, j, i]) == 0)
#                continue
#            else
#                local wall = superpotential_2(conf.phi[:, :, j, i], spt)
#                local col = 0
#                for i2 = 1:size(conf.phi, 2)
#                    local r = 0
#                    for j2 = i2:size(conf.phi, 2)
#                        col += 1
#                        if (i2==i_field)&&(j2==i_field)
#                            r = 1
#                            break
#                        end
#                    end
#                    if r==1
#                        break
#                    end
#                end
#                local w2 = fieldtomatrix(wall[:, col])
#                J[(N+1):2N, 1:N] = w2
#                J[1:N, (N+1):2N] = w2'
#                local Jinv = inv(J)
#                #pp[:, i] += diag(Jinv)[1:N] * sgn[j, i]
#                local u = rs + diag(Jinv)[1:N] * conf.sgn[j, i]#
#                local s = tmp + u#
#                local t = s - tmp#
#                rs = u - t#
#                tmp = s#
#            end
#        end
#        pp[:, i] = tmp#
#    end
#    local pp_sum, pp_dev
#    pp_sum, pp_dev = wt_dev(pp)
#    local Delta, Del_dev
#    Delta, Del_dev = witten_ind(conf.sgn)
#    local pm = abs(p_matrix(N))
#    local res = zeros(Float64, N-1, 3)
#    res[1:Q, 1]       = log(abs2(pm)[1:Q])
#    res[(Q+1):(2Q),1] = log(abs2(pm)[(Q+2):N])
#    res[1:Q, 2]       = log((L-1)^4 * abs(pp_sum[1:Q])     ./ pm[1:Q]     / Delta)
#    res[(Q+1):(2Q),2] = log((L-1)^4 * abs(pp_sum[(Q+2):N]) ./ pm[(Q+2):N] / Delta)
#    local tmp = (real(pp_sum).*real(pp_dev)).^2 + (imag(pp_sum).*imag(pp_dev)).^2
#    res[1:Q, 3]       = sqrt(tmp[1:Q]    ./(abs2(pp_sum[1:Q]    )).^2 + (Del_dev/Delta)^2)
#    res[(Q+1):(2Q),3] = sqrt(tmp[(Q+2):N]./(abs2(pp_sum[(Q+2):N])).^2 + (Del_dev/Delta)^2)
#    return res
#end
#  Calc scaling dimension
function h_fit_range(chi::Array{Float64, 2}, r::Float64, ir_uv::Int)
    local N = size(chi, 1) + 1
    #local Q = div(N, 2)
    #local L = Int(floor(sqrt(N)))
    local res = zeros(Float64, N-1, size(chi, 2))
    local res_size = 0
    for i = 1:(N-1)
        if (ir_uv == 0 && chi[i, 1] >= 2log(r) && chi[i, 1] < 2log(pi))
            res_size += 1
            res[res_size, :] = chi[i, :]
        elseif (ir_uv == 1 && chi[i, 1] < 2log(r))
            res_size += 1
            res[res_size, :] = chi[i, :]
        end
    end
    if res_size == 0
        return zeros(Float64, 1, size(chi, 2))
    end
    return res[1:res_size, :]
end
function h_fit_range(chi::Array{Float64, 2}, r::Int)
    local N = size(chi, 1) + 1
    local Q = div(N, 2)
    #local L = Int(floor(sqrt(N)))
    local res = zeros(Float64, N-1, size(chi, 2))
    local res_size = 0
    local p = abs2(p_index(N))
    for i = 1:N
        if p[i] > 0 && p[i] <= r
            local j
            if i <= Q
                j = i
            elseif i > Q
                j = i-1
            end
            res_size += 1
            res[res_size, :] = chi[j, :]
        end
    end
    return res[1:res_size, :]
end
function h_fit_range(chi::Array{Float64, 2}, rmin::Int, rmax::Int)
    local N = size(chi, 1) + 1
    local Q = div(N, 2)
    #local L = Int(floor(sqrt(N)))
    local res = zeros(Float64, N-1, size(chi, 2))
    local res_size = 0
    local p = abs2(p_index(N))
        #(div(L, 2) - div(i-1, L))^2 + (div(L, 2) - (i-1)%L)^2
    for i = 1:N
        if p[i] >= rmin && p[i] < rmax
            local j
            if i <= Q
                j = i
            elseif i > Q
                j = i-1
            end
            res_size += 1
            res[res_size, :] = chi[j, :]
        end
    end
    return res[1:res_size, :]
end
function h_fit_range(chi::Array{Float64, 2}, rmin::Float64, rmax::Float64)
    local N = size(chi, 1) + 1
    #local Q = div(N, 2)
    #local L = Int(floor(sqrt(N)))
    local res = zeros(Float64, N-1, size(chi, 2))
    local res_size = 0
    for i = 1:(N-1)
        if chi[i, 1] >= 2log(rmin) && chi[i, 1] <= 2log(rmax)
            res_size += 1
            res[res_size, :] = chi[i, :]
        end
    end
    return res[1:res_size, :]
end
function h_fit(chi::Array{Float64, 2})
    local dy2, xdy2, x2dy2, ydy2, y2dy2, xydy2
    dy2   = sum(chi[:, 3].^(-2))
    xdy2  = sum(chi[:, 1] .* (chi[:, 3].^(-2)))
    x2dy2 = sum((chi[:, 1].^2) .* (chi[:, 3].^(-2)))
    ydy2  = sum(chi[:, 2] .* (chi[:, 3].^(-2)))
    y2dy2 = sum((chi[:, 2].^2) .* (chi[:, 3].^(-2)))
    xydy2 = sum(chi[:, 1] .* chi[:, 2] .* (chi[:, 3].^(-2)))
    ## parameters:
    ##  a  delta_a
    ##  b  delta_b
    ##  rsc(reduced chi-square) d.o.f.
    local dn
    local par = zeros(Float64, 3, 2)
    dn = dy2 * x2dy2 - xdy2^2
    par[1,1] = ( dy2 * xydy2 - xdy2 * ydy2 ) / dn
    par[1,2] = sqrt( abs(dy2 / dn ))
    par[2,1] = ( x2dy2 * ydy2 - xdy2 * xydy2 ) / dn
    par[2,2] = sqrt( abs(x2dy2 / dn ))
    par[3,1] = ( y2dy2 - par[1,1] * xydy2 - par[2,1] * ydy2 ) / (size(chi, 1) - 2)
    par[3,2] = size(chi, 1)
    par[1, 1] *= -1
    return par
end

###
###  FUNCTION: Output Correlation Functions
###
#  For scaling dimension
function h_output(chi::Array{Float64, 2}, huv::Array{Float64, 2}, h2::Array{Float64, 2}, h4::Array{Float64, 2}, hr::Array{Float64, 3}, spt::Spt_type, i_field::Int)
    local fout, ffit, fh
    local fn_spt = filename_spt(spt)
    local fn_tmp = string(FN_OUTPUTDIR, "h2/h2_", fn_spt)
    fout = string(fn_tmp, ".dat")
    local ft
    if i_field==1
        ft = "w"
    else
        ft = "a"
    end
    open(fout, ft) do fo
        writedlm(fo, chi)
        write(fo, "\n\n")
    end
    ffit = string(fn_tmp, "_fit.dat")
    open(ffit, ft) do fo
        write(fo, "# Scaling dimension (fitting: pi^2/2<=p^2<pi^2, dof = $(huv[3, 2]))\n")
        write(fo, "h = $(huv[1, 1]); dh = $(huv[1, 2]); b = $(huv[2, 1]); db = $(huv[2, 2]); c = $(huv[3, 1])\n")
        write(fo, "# Scaling dimension (fitting: p^2 = (2pi/L)^2 (1~2), dof = $(h2[3, 2]))\n")
        write(fo, "h2 = $(h2[1, 1]); dh2 = $(h2[1, 2]); b2 = $(h2[2, 1]); db2 = $(h2[2, 2]); c2 = $(h2[3, 1])\n")
        write(fo, "# Scaling dimension (fitting: p^2 = (2pi/L)^2 (1~4), dof = $(h4[3, 2]))\n")
        write(fo, "h4 = $(h4[1, 1]); dh4 = $(h4[1, 2]); b4 = $(h4[2, 1]); db4 = $(h4[2, 2]); c4 = $(h4[3, 1])\n\n\n")
    end
    fh = string(fn_tmp, "_fit")
    open(fh, ft) do fo
        write(fo, "# p h st-err(h) chi^2\n")
        for i = 1:size(hr, 3)
            write(fo, "$(2pi * i /spt.Li)  $(hr[1,1, i])  $(hr[1,2, i])  $(hr[3,1, i])\n")
        end
        write(fo, "\n\n")
    end
    return 0
end
function h_output(h2::Array{Float64, 2}, h4::Array{Float64, 2}, spt::Spt_type, i_field::Int, ft)
    local fh_spt = filename_spt2(spt)
    local fh = string(FN_OUTPUTDIR, "h2/h2_", fh_spt, "f$(i_field)_fit")
    open(fh, ft) do fo
        if ft == "w"
            write(fo, "# L 2pi/L h st-err(h) sy-err(h) st-err(h4) chi^2 chi^2(h4)\n")
        end
        write(fo, "$(spt.Li)  $(2pi/spt.Li)  $(h2[1,1])  $(h2[1,2])  $(h4[1,1])  $(h4[1,2])  $(h2[3,1])  $(h4[3,1])\n")
    end
    return 0
end

###
###  MAIN FUNCTION: Calculation process
###
#  Scaling dimension
function main_h2(spt::Spt_type, i_field::Int, ft)
    #  Read Configuraitons
    conf = read_conf(spt)
    if conf==0
        return 0
    end
    
    #  Scalar field correlation
    local chi, chiuv, chi2, chi4
    local h4  = zeros(Float64, 3, 2)
    local hr  = zeros(Float64, 3, 2, div(spt.Li, 2))
    local huv = zeros(Float64, 3, 2)
    local r4 = 4 #max of (index of p)^2
    local ruv = pi/sqrt(2) #min of ln(p^2)
    chi = h_lnpp(conf, spt, i_field)
    chi4  = h_fit_range(chi, r4)
    chiuv = h_fit_range(chi, ruv, 0)
    h4  = h_fit(chi4)
    huv = h_fit(chiuv)
    for r = 1:div(spt.Li, 2)
        local chi_r
        chi_r  = h_fit_range(chi, r^2, (r+1)^2)
        hr[:, :, r]  = h_fit(chi_r)
    end
    h_output(chi, huv, hr[:,:,1], h4, hr, spt, i_field)
    h_output(hr[:,:,1], h4, spt, i_field, ft)
    return 0
end
function main_h2(num::Int, strn::Int, endn::Int)
    local ft
    for i=strn:endn
        if i == 1
            ft = "w"
        else
            ft = "a"
        end
        local spt = set_spt(i,num)
        if spt==0
            return 0
        end
        local num_f = set_num_f(spt)
        for j=1:num_f
            main_h2(spt, j, ft)
        end
    end
    return 0
end
#########################



######################
### Central charge ###
######################
###
###  FUNCTION: Calc Central Charge (Energy-momentum tensor)
###
#  \psi_2(p) \Bar{\psi}_{\Dot{2}}(q)
@everywhere function c_pp(phi::Array{Complex{Float64}, 4}, spt::Spt_type, num::Int)
    local const N = size(phi, 1)
    local const L = Int(floor(sqrt(N)))
    local p = p_matrix2(N) * (L-1)^2
    local num_f = size(phi, 2)
    local J = zeros(Complex{Float64}, 2N*num_f, 2N*num_f)
    for i = 1:num_f
        local n = (i-1) * 2N
        J[(1+n):(N+n),    (1+n):(N+n)]    = p
        J[(N+1+n):(2N+n), (N+1+n):(2N+n)] = - conj(p)
    end
    local pp = zeros(Complex{Float64}, N, N, num_f, num_f, size(phi, 3))
    for j = 1:size(phi, 3)
        if (norm(phi[:, 1, j, num]) == 0)
            continue
        else
            local wall = superpotential_2(phi[:, :, j, num], spt)
            local col = 1
            for i1 = 1:num_f
                local r = (i1-1) * 2N
                for i2 = i1:num_f
                    local c = (i2-1) * 2N
                    local w2 = fieldtomatrix(wall[:,col])
                    J[(N+1+r):(2N+r), (1+c):(N+c)] = w2
                    J[(1+r):(N+r), (N+1+c):(2N+c)] = w2'
                    if i1!=i2
                        J[(N+1+c):(2N+c), (1+r):(N+r)] = w2
                        J[(1+c):(N+c), (N+1+r):(2N+r)] = w2'
                    end
                    col += 1
                end
            end
            local Jinv = inv(J)
            for i1 = 1:num_f
                local r = (i1-1) * 2N
                for i2 = 1:num_f
                    local c = (i2-1) * 2N
                    pp[:, :, i1, i2, j] =  Jinv[(N+1+r):(2N+r), (N+1+c):(2N+c)]
                end
            end
        end
    end
    return pp
end
#  \langle S^+(p) S^-(-p) \rangle
@everywhere function c_current(pm2::Array{Complex{Float64}, 2}, phi2::Array{Complex{Float64}, 4}, pp::Array{Complex{Float64}, 5}, p::Int)
    local const N = size(phi2, 1)
    local s = zeros(Complex{Float64}, size(phi2, 4))
    for n1 = 1:size(pp, 3)
        for n2 = 1:size(pp, 4)
            s += [( (pm2[p,:].*phi2[p,:,n1,i]).' * pp[:,:,n1,n2,i] * (pm2[p,:].*conj(phi2[p,:,n2,i])) )[1] for i=1:size(phi2,4)]
        end
    end
    return s
end
@everywhere function c_corr(conf::Conf_type, spt::Spt_type, pm2::Array{Complex{Float64}, 2}, num::Int)
    #(phi::Array{Complex{Float64}, 3}, pm2::Array{Complex{Float64}, 2}, sgn::Array{Int64, 2}, num::Int, k::Int, lm::Float64)
    local const N = size(conf.phi, 1)
    local pp = c_pp(conf.phi, spt, num)
    local phi2 = zeros(Complex{Float64}, N, N, size(conf.phi, 2), size(conf.phi, 3))
    for i = 1:size(conf.phi, 2)
        for j = 1:size(conf.phi, 3)
            phi2[:, :, i, j] = fieldtomatrix(conf.phi[:, i, j, num])
        end
    end
    local ss = zeros(Complex{Float64}, N)
    for p = 1:N
        ss[p] = (c_current(pm2, phi2, pp, p).' * conf.sgn[:,num])[1]
    end
    return ss
end
@everywhere function emt_corr(conf::Conf_type, spt::Spt_type, pm::Array{Complex{Float64}, 1}, num::Int)
    #(phi::Array{Complex{Float64}, 3}, pm::Array{Complex{Float64}, 1}, sgn::Array{Int64, 2}, num::Int, k::Int, lm::Float64)
    local pm2 = fieldtomatrix(pm)
    local ss = c_corr(conf, spt, pm2, num)
    local tt = - pm .* (ss[:] - ss[end:-1:1]) / 16
end
function emt_tt(conf::Conf_type, spt::Spt_type)
    #(phi::Array{Complex{Float64}, 3}, sgn::Array{Int64, 2}, k::Int, lm::Float64)
    local const N = size(conf.phi, 1)
    local const L = Int(floor(sqrt(N)))
    local pm = p_matrix(N)
    local tt = zeros(Complex{Float64}, N, size(conf.phi, 4))
    tt = @parallel hcat for j = 1:size(conf.phi, 4)
        emt_corr(conf, spt, pm, j)
    end
    local tt_ave, tt_dev
    tt_ave, tt_dev = wt_dev(tt)
    local Delta, Del_dev
    Delta, Del_dev = witten_ind(conf.sgn)
    local res = (2pi)^2 * tt_ave / Delta
    local res_dev = abs(imag(res))im .* sqrt(abs2(imag(tt_dev)./imag(tt_ave)) + abs2(Del_dev/Delta))
    res_dev += abs(real(res)) .* sqrt(abs2(real(tt_dev)./real(tt_ave)) + abs2(Del_dev/Delta))
    return res, res_dev
end
#  Calc cengral charge
function emt_fit_p(L::Int)
    local p = zeros(Float64, L)
    for i = 1:L
        p[i] = 2pi * (div(L, 2) - (i-1)) / (L-1)
    end
    local res = zeros(Complex{Float64}, L*L)
    for i = 0:(L*L-1)
        if (p[div(i,L)+1]^2 + p[i%L+1]^2)==0
            res[i+1] = 0.0im
        else
            res[i+1] += ( p[div(i,L)+1]^4 - 6p[div(i,L)+1]^2*p[i%L+1]^2 + p[i%L+1]^4 ) / ( p[div(i,L)+1]^2 + p[i%L+1]^2 )
            res[i+1] += 4p[div(i,L)+1]*p[i%L+1] * ( -p[div(i,L)+1]^2 + p[i%L+1]^2 )im / ( p[div(i,L)+1]^2 + p[i%L+1]^2 )
        end
    end
    return res
end
function emt_fit_range(tt::Array{Complex{Float64}, 1}, tt_dev::Array{Complex{Float64}, 1}, p::Array{Complex{Float64}, 1}, pran::Array{Float64, 1}, pmin::Int, pmax::Int)
    local const N = length(tt)
    local tt_r = zeros(Complex{Float64}, N)
    local tt_dev_r = zeros(Complex{Float64}, N)
    local p_r = zeros(Complex{Float64}, N)
    local count = 0
    for i = 1:N
        if pmin <= pran[i] && pran[i] < pmax
            count += 1
            tt_r[count] = tt[i]
            tt_dev_r[count] = tt_dev[i]
            p_r[count] = p[i]
        end
    end
    return tt_r[1:count], tt_dev_r[1:count], p_r[1:count]
end
function emt_fit_range(tt::Array{Complex{Float64}, 1}, tt_dev::Array{Complex{Float64}, 1}, p::Array{Complex{Float64}, 1}, pran::Array{Float64, 1}, pfit::Float64)
    local const N = length(tt)
    local tt_r = zeros(Complex{Float64}, N)
    local tt_dev_r = zeros(Complex{Float64}, N)
    local p_r = zeros(Complex{Float64}, N)
    local count = 0
    for i = 1:N
        if pfit == pran[i]
            count += 1
            tt_r[count] = tt[i]
            tt_dev_r[count] = tt_dev[i]
            p_r[count] = p[i]
        end
    end
    return tt_r[1:count], tt_dev_r[1:count], p_r[1:count]
end
function emt_fit_fit(tt::Array{Complex{Float64}, 1}, tt_dev::Array{Complex{Float64}, 1}, p::Array{Complex{Float64}, 1}, N::Int)
    local const L = Int(floor(sqrt(N)))
    local del = sum(real(p).^2 ./ real(tt_dev).^2 + imag(p).^2 ./ imag(tt_dev).^2)
    local c = sum(real(tt) .* real(p) ./ real(tt_dev).^2 + imag(tt) .* imag(p) ./ imag(tt_dev).^2) / del
    local chi = sum(real(tt).^2 ./ real(tt_dev).^2 + imag(tt).^2 ./ imag(tt_dev).^2) - c^2 * del
    return [c / ((L-1)^2 * pi / 48), 1 / (sqrt(del) * (L-1)^2 * pi / 48), chi / (Float64(2length(tt))-1), Float64(length(tt))]
end
function emt_fit(tt::Array{Complex{Float64}, 1}, tt_dev::Array{Complex{Float64}, 1})
    local const N = length(tt)
    local const L = Int(floor(sqrt(N)))
    local p3 = emt_fit_p(L)
    local pran = abs(p_index(N))
        #sqrt( (div(L, 2) - div(i-1, L))^2 + (div(L, 2) - (i-1)%L)^2 )
    local len_c = div(L, 2)-1 #Int(floor( div(L, 2)*sqrt(2) ))
    local c = zeros(Float64, len_c, 5)#|p|,central charge,deviation,reduced-chi^2,dof
    for i = 1:len_c
        local tt_i, tt_dev_i, p3_i
        tt_i, tt_dev_i, p3_i = emt_fit_range(tt, tt_dev, p3, pran, i, i+1)
        c[i, 1]  = 2pi * i / (L-1)
        c[i,2:5] = emt_fit_fit(tt_i, tt_dev_i, p3_i, N)
    end
    local c_sys = zeros(Float64, 2, 5)
    for i = 1:2
        local tt_i, tt_dev_i, p3_i
        tt_i, tt_dev_i, p3_i = emt_fit_range(tt, tt_dev, p3, pran, sqrt(i))
        c_sys[i, 1]  = 2pi * sqrt(i) / (L-1)
        c_sys[i,2:5] = emt_fit_fit(tt_i, tt_dev_i, p3_i, N)
    end
    return c, c_sys
end

###
###  FUNCTION: Output Correlation Functions
###
#  For central charge
function emt_output(tt::Array{Complex{Float64}, 1}, tt_dev::Array{Complex{Float64}, 1}, c::Array{Float64, 2}, c_sys::Array{Float64, 2}, spt::Spt_type, ft)
    local fn_spt = filename_spt(spt)
    local fn_tmp = string(FN_OUTPUTDIR, "emt/emt_", fn_spt)
    ffit = string(fn_tmp, "_fit")
    open(ffit, "w") do fo
        for i = 1:size(c_sys, 1)
            write(fo, "# ")
            writedlm(fo, c_sys[i, :].')
        end
        writedlm(fo, c)
    end
    local pm = p_matrix((spt.Li + 1)^2)
    local num_col = 5
    local corr = zeros(Float64, (spt.Li + 1)^2, num_col)
    corr[:, 1] = imag(pm)
    corr[:, 2] = real(tt)
    corr[:, 3] = real(tt_dev)
    corr[:, 4] = imag(tt)
    corr[:, 5] = imag(tt_dev)
    local fout0, fout1
    fout0 = string(fn_tmp, "_p1fix.dat")
    open(fout0, "w") do fo
        for p1 = 1:(spt.Li+1)
            writedlm(fo, [corr[(spt.Li+1) * i + p1, j] for i=0:spt.Li, j=1:num_col])
            write(fo, "\n\n")
        end
    end
    fout1 = string(fn_tmp, "_p0fix.dat")
    corr[:, 1]  = real(pm)
    open(fout1, "w") do fo
        for p0 = 0:spt.Li
            writedlm(fo, [corr[i + (spt.Li+1) * p0, j] for i=1:(spt.Li+1), j=1:num_col])
            write(fo, "\n\n")
        end
    end
    local fall
    local fn_spt2 = filename_spt2(spt)
    fall = string(FN_OUTPUTDIR, "emt/emt_", fn_spt2, "_fit")
    open(fall, ft) do fo
        if ft == "w"
            write(fo, "# L c st-err sys-err1 sys-err2 chi^2\n")
        end
        write(fo, "$(spt.Li)  $(c[1,2])  $(c[1,3])  $(abs(c[1,2]-c_sys[1,2]))  $(abs(c[1,2]-c_sys[2,2]))  $(c[1,4])\n")
    end
    return 0
end

###
###  MAIN FUNCTION: Calculation process
###
#  Central charge
function main_emt(spt::Spt_type, ft)
    #  Read Configuraitons
    conf = read_conf(spt)
    
    #  Central charge
    local tt, tt_dev, c, c_sys
    tt, tt_dev = emt_tt(conf, spt)
    c, c_sys = emt_fit(tt, tt_dev)
    emt_output(tt, tt_dev, c, c_sys, spt, ft)
    return 0
end
function main_emt(num::Int, strn::Int, endn::Int)
    local ft
    for i=strn:endn
        if i == 1
            ft = "w"
        else
            ft = "a"
        end
        local spt = set_spt(i, num)
        if spt==0
            return 0
        end
        main_emt(spt, ft)
    end
    return 0
end
######################



##############################
### C-function in progress ###
##############################
###
###  FUNCTION: Calc C-function (Energy-momentum tensor)
###
#  \psi_1(p) \Bar\psi_{\Dot{1}}(q)             & \psi_1(p) \psi_2(q)
#  \Bar\psi_{\Dot{2}}(p) \Bar\psi_{\Dot{1}}(q) & \Bar\psi_{\Dot{2}}(p) \psi_2(q)
@everywhere function cfunc_pp(phi::Array{Complex{Float64}, 4}, spt::Spt_type, num::Int)
    local const N = size(phi, 1)
    local const L = Int(floor(sqrt(N)))
    local p = p_matrix2(N) * (L-1)^2
    local num_f = size(phi, 2)
    local J = zeros(Complex{Float64}, 2N*num_f, 2N*num_f)
    for i = 1:num_f
        local n = (i-1) * 2N
        J[(1+n):(N+n),    (1+n):(N+n)]    = p
        J[(N+1+n):(2N+n), (N+1+n):(2N+n)] = - conj(p)
    end
    local pp = zeros(Complex{Float64}, 2N, 2N, num_f, num_f, size(phi, 3))
    for j = 1:size(phi, 3)
        if (norm(phi[:, 1, j, num]) == 0)
            continue
        else
            local wall = superpotential_2(phi[:, :, j, num], spt)
            local col = 1
            for i1 = 1:num_f
                local r = (i1-1) * 2N
                for i2 = i1:num_f
                    local c = (i2-1) * 2N
                    local w2 = fieldtomatrix(wall[:,col])
                    J[(N+1+r):(2N+r), (1+c):(N+c)] = w2
                    J[(1+r):(N+r), (N+1+c):(2N+c)] = w2'
                    if i1!=i2
                        J[(N+1+c):(2N+c), (1+r):(N+r)] = w2
                        J[(1+c):(N+c), (N+1+r):(2N+r)] = w2'
                    end
                    col += 1
                end
            end
            local Jinv = inv(J)
            for i1 = 1:num_f
                local r = (i1-1) * 2N
                for i2 = 1:num_f
                    local c = (i2-1) * 2N
                    pp[:, :, i1, i2, j] =  Jinv[(1+r):(2N+r), (1+c):(2N+c)]
                end
            end
        end
    end
    return pp
end
#  S_{\Bar{z}}^+(p) S_{\Bar{z}}^-(-p) & S_{\Bar{z}}^+(p) S_z^-(-p)
#  S_z^+(p) S_{\Bar{z}}^-(-p)         & S_z^+(p) S_z^-(-p)
@everywhere function cfunc_current(pm2::Array{Complex{Float64}, 2}, nic2::Array{Complex{Float64}, 3}, phi2::Array{Complex{Float64}, 4}, pp::Array{Complex{Float64}, 5}, p::Int)
    local const N = size(phi2, 1)
    local s11 = zeros(Complex{Float64}, size(phi2, 4))
    local s12 = zeros(Complex{Float64}, size(phi2, 4))
    local s21 = zeros(Complex{Float64}, size(phi2, 4))
    local s22 = zeros(Complex{Float64}, size(phi2, 4))
    for n1 = 1:size(pp, 3)
        for n2 = 1:size(pp, 4)
            s11 += [( (conj(nic2[p,:,n1])+pm2[p,:].*conj(phi2[p,:,n1,i])).' * pp[1:N,1:N,n1,n2,i] * (nic2[p,:,n1]-pm2[p,:].*phi2[p,:,n2,i]) )[1] for i=1:size(phi2,4)]
            s12 += [( (conj(nic2[p,:,n1])+pm2[p,:].*conj(phi2[p,:,n1,i])).' * pp[1:N,(N+1):2N,n1,n2,i] * (pm2[p,:].*conj(phi2[p,:,n2,i])) )[1] for i=1:size(phi2,4)]
            s21 += [( (pm2[p,:].*phi2[p,:,n1,i]).' * pp[(N+1):2N,1:N,n1,n2,i] * (nic2[p,:,n1]-pm2[p,:].*phi2[p,:,n2,i]) )[1] for i=1:size(phi2,4)]
            s22 += [( (pm2[p,:].*phi2[p,:,n1,i]).' * pp[(N+1):2N,(N+1):2N,n1,n2,i] * (pm2[p,:].*conj(phi2[p,:,n2,i])) )[1] for i=1:size(phi2,4)]
        end
    end
    return s11, s12, s21, s22
end
@everywhere function cfunc_corr(conf::Conf_type, spt::Spt_type, pm2::Array{Complex{Float64}, 2}, num::Int)
    local const N = size(conf.phi, 1)
    local pp = cfunc_pp(conf.phi, spt, num)
    local phi2 = zeros(Complex{Float64}, N, N, size(conf.phi, 2), size(conf.phi, 3))
    local nic2 = zeros(Complex{Float64}, N, N, size(conf.nic, 2))
    for i = 1:size(conf.phi, 2)
        for j = 1:size(conf.phi, 3)
            phi2[:, :, i, j] = fieldtomatrix(conf.phi[:, i, j, num])
        end
        nic2[:, :, i] = fieldtomatrix(conf.nic[:, i, num])
    end
    local ss11 = zeros(Complex{Float64}, N)
    local ss12 = zeros(Complex{Float64}, N)
    local ss21 = zeros(Complex{Float64}, N)
    local ss22 = zeros(Complex{Float64}, N)
    for p = 1:N
        local tmp11, tmp12, tmp21, tmp22
        tmp11,tmp12,tmp21,tmp22 = cfunc_current(pm2, nic2, phi2, pp, p)
        ss11[p] = (tmp11.' * conf.sgn[:,num])[1]
        ss12[p] = (tmp12.' * conf.sgn[:,num])[1]
        ss21[p] = (tmp21.' * conf.sgn[:,num])[1]
        ss22[p] = (tmp22.' * conf.sgn[:,num])[1]
    end
    return ss11, ss12, ss21, ss22
end
@everywhere function emtfunc_corr(conf::Conf_type, spt::Spt_type, pm::Array{Complex{Float64}, 1}, num::Int)
    local pm2 = fieldtomatrix(pm)
    local ss11, ss12, s21, s22
    ss11,ss12,ss21,ss22 = cfunc_corr(conf, spt, pm2, num)
    local N = size(ss11, 1)
    local tt = zeros(Complex{Float64}, 3N)
    tt[1:N]       = - pm .* (ss22[:] - ss22[end:-1:1]) / 16
    tt[(N+1):2N]  = - pm .* (ss21[:] - ss12[end:-1:1]) / 16
    tt[(2N+1):3N] = - pm .* (ss11[:] - ss11[end:-1:1]) / 16
    return tt
end
function emtfunc_tt(conf::Conf_type, spt::Spt_type)
    local const N = size(conf.phi, 1)
    local const L = Int(floor(sqrt(N)))
    local pm = p_matrix(N)
    local tt = zeros(Complex{Float64}, 3N, size(conf.phi, 4))
    println("      Parallell computing start")
    tt = @parallel hcat for j = 1:size(conf.phi, 4)
        emtfunc_corr(conf, spt, pm, j)
    end
    println("      Parallell computing end")
    local tt_ave = zeros(Complex{Float64}, N, 3)
    local tt_dev = zeros(Complex{Float64}, N, 3)
    for i=1:3
        tt_ave[:, i], tt_dev[:, i] = wt_dev(tt[((i-1)*N+1):i*N, :])
    end
    local Delta, Del_dev
    Delta, Del_dev = witten_ind(conf.sgn)
    local res = (2pi)^2 * tt_ave / Delta
    local res_dev = (2pi)^2 * sqrt.(abs2.(imag(tt_dev)) + abs2.(imag(tt_ave)) .* abs2.(Del_dev/Delta))im / Delta
    res_dev += (2pi)^2 * sqrt.(abs2.(real(tt_dev)) + abs2.(real(tt_ave)) .* abs2.(Del_dev/Delta)) / Delta
    return res, res_dev
end
#  <phi_1(p)phi_2(-p)> --> <phi_1(x_0,x_1=0)phi_2(0)>
#  x_0 = i * Li / x_steps, where i = 1:x_steps
function momentum2space(corr::Array{Complex{Float64}, 1}, corr_dev::Array{Complex{Float64}, 1}, stp::Int)
    local const N = size(corr, 1)
    local const Li = Int(floor(sqrt(N))) - 1
    local pm = imag(p_matrix(N))
    local res = zeros(Complex{Float64}, stp)
    local res_dev = zeros(Complex{Float64}, stp)
    for i=1:stp
        local x = Float64(i*Li) / stp
        res[i] = ( (exp.(im*x*pm)).' * corr )[] / Li^4
        local dev_imag = ( abs2.(real(exp.(im*x*pm))).' * abs2.(imag(corr_dev)) )[1] + ( abs2.(imag(exp.(im*x*pm))).' * abs2.(real(corr_dev)) )[1]
        local dev_real = ( abs2.(real(exp.(im*x*pm))).' * abs2.(real(corr_dev)) )[1] + ( abs2.(imag(exp.(im*x*pm))).' * abs2.(imag(corr_dev)) )[1]
        res_dev[i] =  sqrt(dev_real) / Li^4 + sqrt(dev_imag)im / Li^4
    end
    return res, res_dev
end
function momentum2space(corr::Array{Complex{Float64}, 2}, corr_dev::Array{Complex{Float64}, 2}, stp::Int)
    local res = zeros(Complex{Float64}, stp, size(corr,2))
    local res_dev = zeros(Complex{Float64}, stp, size(corr,2))
    for i=1:size(corr,2)
        res[:, i],res_dev[:, i] = momentum2space(corr[:,i], corr_dev[:,i], stp)
    end
    return res, res_dev
end
function cfunc_summation(corrx::Array{Complex{Float64}, 2}, corrx_dev::Array{Complex{Float64}, 2}, spt::Spt_type, xstp::Int)
    local res = zeros(Complex{Float64}, xstp)
    local res_dev = zeros(Complex{Float64}, xstp)
    for i=1:xstp
        local x = Float64(i * spt.Li) / xstp
        local f = x^4 * corrx[i, 1]
        local f_dev = x^4 * corrx_dev[i, 1]
        local g = 4 * x^4 * corrx[i, 2]
        local g_dev = 4 * x^4 * corrx_dev[i, 2]
        local h = x^4 * corrx[i, 3]
        local h_dev = x^4 * corrx[i, 3]
        res[i] = 2f - g - 3h/8
        res_dev[i] = sqrt(abs2(imag(2f_dev)) + abs2(imag(g_dev)) + abs2(imag(3h_dev/8)))im
        res_dev[i] += sqrt(abs2(real(2f_dev)) + abs2(real(g_dev)) + abs2(real(3h_dev/8)))
    end
    return res, res_dev
end

###
###  FUNCTION: Output Correlation Functions
###
#  For C-function
function cfunc_p_output(tt::Array{Complex{Float64}, 2}, tt_dev::Array{Complex{Float64}, 2}, spt::Spt_type)
    local fn_spt = filename_spt(spt)
    local fn_tmp = string(FN_OUTPUTDIR, "cfunc/cfunc_p_", fn_spt)
    local pm = p_matrix((spt.Li + 1)^2)
    local num_col = 5
    local corr = zeros(Float64, (spt.Li + 1)^2, num_col, 3)
    for i=1:3
        corr[:, 1, i] = imag(pm)
        corr[:, 2, i] = real(tt[:, i])
        corr[:, 3, i] = real(tt_dev[:, i])
        corr[:, 4, i] = imag(tt[:, i])
        corr[:, 5, i] = imag(tt_dev[:, i])
    end
    local fout0, fout1
    for k=1:3
        local corr_type = ""
        if (k==1)
            corr_type = "f"
        elseif (k==2)
            corr_type = "g"
        elseif (k==3)
            corr_type = "h"
        end
        fout0 = string(fn_tmp, corr_type, "_p1fix.dat")
        open(fout0, "w") do fo
            for p1 = 1:(spt.Li+1)
                writedlm(fo, [corr[(spt.Li+1)*i+p1,j,k] for i=0:spt.Li, j=1:num_col])
                write(fo, "\n\n")
            end
        end
        fout1 = string(fn_tmp, corr_type, "_p0fix.dat")
        corr[:, 1, k]  = real(pm)
        open(fout1, "w") do fo
            for p0 = 0:spt.Li
                writedlm(fo, [corr[i+(spt.Li+1)*p0,j,k] for i=1:(spt.Li+1), j=1:num_col])
                write(fo, "\n\n")
            end
        end
    end
    return 0
end
function cfunc_x_output(tt::Array{Complex{Float64}, 2}, tt_dev::Array{Complex{Float64}, 2}, spt::Spt_type, xsteps::Int)
    local fn_spt = filename_spt(spt)
    local fn_tmp = string(FN_OUTPUTDIR, "cfunc/cfunc_x_", fn_spt)
    local num_col = 5
    local cfunc_tmp, cfunc_dev_tmp
    cfunc_tmp,cfunc_dev_tmp = cfunc_summation(tt, tt_dev, spt, xsteps)
    local corr = zeros(Float64, size(tt, 1), num_col, 3)
    local cfunc = zeros(Float64, size(tt, 1), num_col)
    for i=1:size(cfunc, 1)
        cfunc[i, 1]   = Float64(i*spt.Li) / xsteps
    end
    cfunc[:, 2] = real(cfunc_tmp)
    cfunc[:, 3] = real(cfunc_dev_tmp)
    cfunc[:, 4] = imag(cfunc_tmp)
    cfunc[:, 5] = imag(cfunc_dev_tmp)
    for i=1:3
        for j=1:size(corr, 1)
            corr[j, 1, i] = Float64(j*spt.Li) / xsteps
        end
        corr[:, 2, i] = real(tt[:, i])
        corr[:, 3, i] = real(tt_dev[:, i])
        corr[:, 4, i] = imag(tt[:, i])
        corr[:, 5, i] = imag(tt_dev[:, i])
    end
    local fout = ""
    for k=1:3
        local corr_type = ""
        if (k==1)
            corr_type = "f_"
        elseif (k==2)
            corr_type = "g_"
        elseif (k==3)
            corr_type = "h_"
        end
        fout = string(fn_tmp, corr_type, "x$(xsteps).dat")
        open(fout, "w") do fo
            writedlm(fo, corr[:,:,k])
        end
    end
    local fout_cfunc = string(fn_tmp, "x$(xsteps).dat")
    open(fout_cfunc, "w") do fo
        writedlm(fo, cfunc)
    end
    return 0
end

###
###  MAIN FUNCTION: Calculation process
###
#  Central charge
function main_cfunc(spt::Spt_type, xsteps::Int, ft)
    println("  Calc start")
    #  Read Configuraitons
    println("    Read confs")
    conf = read_conf(spt)
    
    #  <T_zzT_zz>, <T_zzT_z\z>, <T_z\zTz\z>
    println("    EMT correlator")
    local tt, tt_dev
    tt, tt_dev = emtfunc_tt(conf, spt)
    cfunc_p_output(tt, tt_dev, spt)
    println("    EMT in x-space")
    local tt_x, tt_x_dev
    tt_x, tt_x_dev = momentum2space(tt, tt_dev, xsteps)
    cfunc_x_output(tt_x, tt_x_dev, spt, xsteps)
    println("  Calc end")
    return 0
end
function main_cfunc(num::Int, strn::Int, endn::Int, xsteps::Int)
    println("Main start")
    local ft
    for i=strn:endn
        println("Conf number: $i")
        if i == 1
            ft = "w"
        else
            ft = "a"
        end
        local spt = set_spt(i, num)
        if spt==0
            return 0
        end
        main_cfunc(spt, xsteps, ft)
    end
    println("Main end")
    return 0
end
##############################


#main_wt(1, 5,5)
#main_cls(1)

#main_h2(1, 5,5)
#main_emt(1, 5,5)

