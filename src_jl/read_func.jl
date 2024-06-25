###########################################################
##                                                       ##
##  Julia: Analyzing Wess--Zumino model with ADE-type W  ##
##                                                       ##
##    Cervo (landau, lifshitz2, lifshitz3)               ##
##    input-dir:                                         ##
FN_INPUTDIR  = "/home/morikawa/lg-cutoff_data2/confs/"   ##
FN_INPUTDIR2  = "/home/morikawa/lg-cutoff_data2/"        ##
##    output-dir:                                        ##
FN_OUTPUTDIR = "/home/morikawa/Dropbox/lg-cfunc/corr/"   ##
##                                                       ##
###########################################################

#################
### Structure ###
#################
@everywhere type Spt_type
    Li::Int
    spt::Int
    k::Array{Int, 1}
    lm::Array{Float64, 1}
    diglm::Array{Float64, 1}
end
@everywhere type Conf_type
    sgn::Array{Int64, 2}
    dmp::Array{Int64, 1}
    nic::Array{Complex{Float64}, 3}
    phi::Array{Complex{Float64}, 4}
    max_res::Float64
end

###
###  FUNCTION: Set some paramters
###
function set_spt(i_L::Int, num::Int)
    if num==1
        if (i_L == 1)
            return Spt_type(24, 0, [3], [0.2135], [4])
        elseif (i_L == 2)
            return Spt_type(28, 0, [3], [0.2538], [4])
        elseif (i_L == 3)
            return Spt_type(32, 0, [3], [0.3], [1])
        elseif (i_L == 4)
            return Spt_type(36, 0, [3], [0.342], [3])
        elseif (i_L == 5)
            return Spt_type(40, 0, [3], [0.3888], [4])
        elseif (i_L == 6)
            return Spt_type(44, 0, [3], [0.45], [2])
        elseif (i_L == 7)
            return Spt_type(48, 0, [3], [0.51], [2])
        elseif (i_L == 8)
            return Spt_type(52, 0, [3], [0.5705], [4])
        end
    else
        if (i_L == 1)
            return Spt_type(12, 0, [3], [0.2135], [4])
        end
    end
    return 0
end
function set_num_sol(spt::Spt_type)
    return 6
end
function set_num_nic(spt::Spt_type)
    if (spt.Li == 12)
        return 5120
    elseif (spt.Li <= 28)
        return 5120
    elseif (spt.Li == 32)
        return 2592
    elseif (spt.Li >= 36)
        return 2592
    end
    return 0
end
function set_num_nic_d(spt::Spt_type)
    if (spt.Li <= 28)
        return 2560
    elseif (spt.Li >= 32 && spt.Li <= 40)
        return 288
    elseif (spt.Li == 44)
        return 144
    elseif (spt.Li >= 48)
        return 36
    end
    return 0
end
function set_num_f(spt::Spt_type)
    return 1
end

#######################
### Basic FUNCTIONS ###
#######################
###
###  FUNCTION: From a Complex number in C++ to a Complex number in Julia
###
function cpptojulia_comp(s)
    s = replace(replace(replace(s, " ", ""), "(", ""), ")", ",")
    local sl = split(s[1:end-1], ',')
    local num = Int(floor(length(sl)/2))
    local res = zeros(Complex{Float64}, num)
    for i = 1:num
        res[i] = parse(Float64, sl[2i-1]) + parse(Float64, sl[2i])im
    end
    return res
end

###
###  FUNCTION: Read Configurations (Nicolai and Scalar)
###
#  Read a Nicolai Map Configuration
function read_nic(fn, sgn::Array{Int64, 2}, dmp::Array{Int64, 1}, nic::Array{Complex{Float64}, 3}, i_n::Int)
    local res = 0.0
    open(fn) do fi
        for line in eachline(fi)
            if line == "# Signs of Determinants\n"
                local signs = split(readline(fi))
                [sgn[i, i_n] = parse(Int, signs[i]) for i = 1:length(signs)]
            elseif line == "# Max norm of residue\n"
                res = parse(Float64, readline(fi))
            elseif line == "# Configuration of Gaussian random numbers\n"
                for i = 1:size(nic, 1)
                    nic[i, :, i_n] = cpptojulia_comp(readline(fi))
                end
            elseif line == "# Number of Solutions and Dumped-Solutions\n"
                local items = split(readline(fi), ",")
                local tmp   = split(items[2], "/")
                dmp[i_n] = parse(Int, tmp[1])
            end
        end
    end
    return res
end
#  Read a Scalar Field Configuration
function read_phi(fp, phi::Array{Complex{Float64}, 4}, ip::Int, i_n::Int)
    local sign = 1
    local res  = 1.0
    open(fp) do fi
        for line in eachline(fi)
            if line == "# Sign of Determinant"
                sign == parse(Int, readline(fi))
            elseif line == "# Norm of residue\n"
                res = parse(Float64, readline(fi))
            elseif line == "# Configuration of a scalar field solution\n"
                for i = 1:size(phi, 1)
                    phi[i, :, ip, i_n] = cpptojulia_comp(readline(fi))
                end
            end
        end
    end
    return sign * res
end
#  File name template
function filename_spt(spt::Spt_type)
    #(Li::Int, spt_type::Int, k::Array{Int, 1}, lm::Array{Float64, 1}, diglm::Array{Int, 1})
    local fn_spt = ""
    local fn_lm  = ""
    local fn_k   = ""
    local zs = Array(String, length(spt.diglm))
    for i = 1:length(zs)
        zs[i] = ""
        for digit = 1:(6-spt.diglm[i])
            zs[i] = string(zs[i], "0")
        end
    end
    if spt.spt==0
        #fn_spt = "spt-a_"
        fn_spt = ""
    elseif spt.spt==1
        fn_spt = "spt-d_"
    elseif spt.spt==2
        fn_spt = "spt-e_"
    elseif spt.spt==3
        fn_spt = "spt-c_"
    end
    if length(spt.lm)==1
        fn_lm = string("lm$(spt.lm[1])", zs[1])
    else
        for i = 1:length(spt.lm)
            fn_lm = string(fn_lm, "lm$(i-1)-$(spt.lm[i])", zs[i])
        end
    end
    if length(spt.k)==1
        fn_k = "k$(spt.k[1])"
    else
        for i = 1:length(spt.k)
            fn_k = string(fn_k, "k$(i-1)-$(spt.k[i])")
        end
    end
    return string(fn_spt, fn_lm, fn_k, "L$(spt.Li)")
end
function filename_spt2(spt::Spt_type)
    #(Li::Int, spt_type::Int, k::Array{Int, 1}, lm::Array{Float64, 1}, diglm::Array{Int, 1})
    local fn_spt = ""
    local fn_lm  = ""
    local fn_k   = ""
    local zs = Array(String, length(spt.diglm))
    for i = 1:length(zs)
        zs[i] = ""
        for digit = 1:(6-spt.diglm[i])
            zs[i] = string(zs[i], "0")
        end
    end
    if spt.spt==0
        #fn_spt = "spt-a_"
        fn_spt = ""
    elseif spt.spt==1
        fn_spt = "spt-d_"
    elseif spt.spt==2
        fn_spt = "spt-e_"
    elseif spt.spt==3
        fn_spt = "spt-c_"
    end
    if length(spt.lm)==1
        fn_lm = string("lm$(spt.lm[1])", zs[1])
    else
        for i = 1:length(spt.lm)
            fn_lm = string(fn_lm, "lm$(i-1)-$(spt.lm[i])", zs[i])
        end
    end
    if length(spt.k)==1
        fn_k = "k$(spt.k[1])"
    else
        for i = 1:length(spt.k)
            fn_k = string(fn_k, "k$(i-1)-$(spt.k[i])")
        end
    end
    return string(fn_spt, fn_lm, fn_k)
end
function filename_spt3(spt::Spt_type)
    #(Li::Int, spt_type::Int, k::Array{Int, 1}, lm::Array{Float64, 1}, diglm::Array{Int, 1})
    local fn_spt = ""
    local fn_lm  = ""
    local fn_k   = ""
    if spt.spt==0
        #fn_spt = "spt-a_"
        fn_spt = ""
    elseif spt.spt==1
        fn_spt = "spt-d_"
    elseif spt.spt==2
        fn_spt = "spt-e_"
    elseif spt.spt==3
        fn_spt = "spt-c_"
    end
    if length(spt.lm)==1
        fn_lm = string("lm$(spt.lm[1])")
    else
        for i = 1:length(spt.lm)
            fn_lm = string(fn_lm, "lm$(i-1)-$(spt.lm[i])")
        end
    end
    if length(spt.k)==1
        fn_k = "k$(spt.k[1])"
    else
        for i = 1:length(spt.k)
            fn_k = string(fn_k, "k$(i-1)-$(spt.k[i])")
        end
    end
    return string(fn_spt, fn_lm, fn_k, "L$(spt.Li)")
end

#  Read Configurations
function read_conf(spt::Spt_type)
    local num_sol = set_num_sol(spt)
    local num_nic = set_num_nic(spt)
    local num_f   = set_num_f(spt)
    if num_sol == 0
        return 0
    end
    local sgn = zeros(Int, num_sol, num_nic)
    local dmp = zeros(Int, num_nic)
    local nic = zeros(Complex{Float64}, (spt.Li + 1)^2, num_f, num_nic )
    local phi = zeros(Complex{Float64}, (spt.Li + 1)^2, num_f, num_sol, num_nic )
    
    local max_res = 0.0
    local fn, fp
    local fn_spt = filename_spt(spt)
    for i_n = 1:size(sgn, 2)
        fn = string(FN_INPUTDIR, "nicolai_", fn_spt, "n$(i_n-1).dat")
        local res = read_nic(fn, sgn, dmp, nic, i_n)
        if res > max_res
            max_res = res
        end
        
        local ip = 1
        while ip <= size(sgn, 1) && sgn[ip, i_n] != 0
            fp = string(FN_INPUTDIR, "/phi_", fn_spt, "n$(i_n-1)_$(ip-1).dat")
            read_phi(fp, phi, ip, i_n)
            ip += 1
        end
    end
    return Conf_type(sgn, dmp, nic, phi, max_res)
end


###
### Average & Standard Deviation
###
function wt_dev(v::Array{Complex{Float64}, 2})
    local v_sum = mean(v, 2)[:]
    local v_dff = v - v_sum * ones(1, size(v, 2))
    local v_var = mean(abs2(real(v_dff)), 2)[:] + mean(abs2(imag(v_dff)), 2)[:]im
    local v_dev = sqrt(real(v_var)) + sqrt(imag(v_var))im
    return v_sum, v_dev / sqrt(size(v, 2)-1)
end
function wt_dev(v::Array{Float64, 2})
    local v_sum = mean(v, 2)[:]
    local v_dff = v - v_sum * ones(1, size(v, 2))
    local v_var = mean(abs2(v_dff), 2)[:]
    local v_dev = sqrt(v_var)
    return v_sum, v_dev / sqrt(size(v, 2)-1)
end

###
### Witten index
###
function witten_ind(sgn::Array{Int64, 2})
    local sgn_sum = sum(sgn, 1)[:]
    local wit = mean(sgn_sum)
    local wit_dev = sqrt(mean(abs2(sgn_sum - wit)) / (size(sgn,2)-1))
    return wit, wit_dev
end
function witten_ind(sgn::Array{Int64, 1})
    local wit = mean(sgn)
    local wit_dev = sqrt(mean(abs2(sgn - wit)) / (length(sgn)-1))
    return wit, wit_dev
end
#######################


##############################################
### Basic FUNCTIONS for Parallel computing ###
##############################################
###
###  FUNCTION: Calculate Matrices (Fields) as C++
###
#  V(p - q) -> M(p, q)
@everywhere function fieldtomatrix(v::Array{Complex{Float64}, 1})
    local const N = length(v)
    local const L = Int(floor(sqrt(N)))
    local m = zeros(Complex{Float64}, N, N)
    for i = 1:N
        for j = 1:N
            local ind_0 = div(L, 2) + div(i-1, L) - div(j-1, L)
            local ind_1 = div(L, 2) + (i-1)%L     - (j-1)%L
            if (0 <= ind_0 && ind_0 < L && 0 <= ind_1 && ind_1 < L)
                m[i, j] = v[L*ind_0 + ind_1 + 1]
            end
        end
    end
    return m
end

#  Convolution
##  v1 * v2
@everywhere function convolution(v1::Array{Complex{Float64}, 1}, v2::Array{Complex{Float64}, 1})
    local const N = length(v1)
    local const Li = Int(floor(sqrt(N))) - 1
    local m1 = fieldtomatrix(v1)
    return m1 * v2 / Li^2
end
## (v * v * ... * v) where #v = pw
@everywhere function convolution(v::Array{Complex{Float64}, 1}, pw::Int)
    local const N = length(v)
    local const Li = Int(floor(sqrt(N))) - 1
    if pw<0
        return zeros(Complex{Float64}, N)
    elseif pw==0
        local res = zeros(Complex{Float64}, N)
        res[Int(floor(N/2))+1] = Complex{Float64}(Li^2)
        return res
    elseif pw==1
        return v
    elseif pw==2
        return convolution(v, v)
    else
        local m  = fieldtomatrix(v)
        local mp = eye(N, N)
        for i = 1:(pw-2)
            mp = m * mp
        end
        return m * mp * v / (Li^2)^(pw-1)
    end
end
#  p_mat : 2 i p_z
@everywhere function p_matrix(N::Int)
    local const L = Int(floor(sqrt(N)))
    local p = zeros(Complex{Float64}, N)
    for i = 1:N
        p[i] = 2pi * (div(L, 2)im - div(i-1, L)im + div(L, 2) - (i-1)%L) / (L-1)
    end
    return p
end
@everywhere function p_matrix2(N::Int)
    local const L = Int(floor(sqrt(N)))
    local p = zeros(Complex{Float64}, N, N)
    for i = 1:N
        p[i, i] = 2pi * (div(L, 2)im - div(i-1, L)im + div(L, 2) - (i-1)%L) / (L-1)
    end
    return p
end
@everywhere function p_index(N::Int)
    local const L = Int(floor(sqrt(N)))
    local p = zeros(Complex{Int}, N)
    for i = 1:N
        p[i] = div(L, 2)im - div(i-1, L)im + div(L, 2) - (i-1)%L
    end
    return p
end

#  Superpotential and its Derivatives
@everywhere function superpotential_2(phi::Array{Complex{Float64}, 2}, spt::Spt_type)
    if spt.spt==0
        return superpotential_typealgebraA(phi, spt)
    elseif spt.spt==1
        return superpotential_typealgebraD(phi, spt)
    elseif spt.spt==2
        return superpotential_typealgebraE(phi, spt)
    elseif spt.spt==3
        return superpotential_typecustom(phi, spt)
    end
    return 0
end
@everywhere function superpotential_typealgebraA(phi::Array{Complex{Float64}, 2}, spt::Spt_type)
    local w = convolution(phi[:,1], spt.k[1]-2) * (spt.k[1]-1) * spt.lm[1]
    return hcat(w)
end
@everywhere function superpotential_typealgebraD(phi::Array{Complex{Float64}, 2}, spt::Spt_type)
    local xk2 = convolution(phi[:,1], spt.k[1]-2) * (spt.k[1]-1) * spt.lm[1]
    local x1 = phi[:, 1] * spt.lm[2]
    local y1 = phi[:, 2] * spt.lm[2]
    return hcat(xk2, y1, x1)
end
@everywhere function superpotential_typealgebraE(phi::Array{Complex{Float64}, 2}, spt::Spt_type)
    local x1 = phi[:, 1] * 2.0 * spt.lm[1]
    local xy = convolution(phi[:,1], phi[:,2]) * 2.0 * spt.lm[2]
    local y2 = convolution(phi[:,2], 2) * spt.lm[2]
    return hcat(x1, y2, xy)
end
@everywhere function superpotential_typecustom(phi::Array{Complex{Float64}, 2}, spt::Spt_type)
    local x1 = phi[:,1] * 2.0 * spt.lm[1]
    local y1 = phi[:,2] * 2.0 * spt.lm[2]
    local z1 = phi[:,3] * 2.0 * spt.lm[3]
    local xlm = phi[:,1] * spt.lm[4]
    local ylm = phi[:,2] * spt.lm[4]
    local zlm = phi[:,3] * spt.lm[4]
    return hcat(x1, zlm, ylm, y1, xlm, z1)
end
##############################################

##################################################
### FUCNTIONS for Ward--Takahashi identitities ###
##################################################
##
##  FUNCTION: Read ...
##
#  Classification of confs
function read_cls(sgn::Array{Int64, 2})
    local count = zeros(Int, 3, size(sgn, 2))
    local count_size = 0
    for i = 1:size(sgn, 2)
        local c1 = sum(abs(sgn[:, i]))
        local c2 = sum(    sgn[:, i] )
        local a  = 0
        for j = 1:count_size
            if count[1, j] == c1 && count[2, j] == c2
                count[3, j] += 1
                a = 1
                break
            end
        end
        if a == 0
            count_size += 1
            count[1, count_size] = c1
            count[2, count_size] = c2
            count[3, count_size] += 1
        end
    end
    return count[:, 1:count_size]
end
#  Read Output_file to get time
function read_out(spt::Spt_type, nic_str::Int, nic_end::Int)
    #local fn = "outputs/lm$(lm)k$(k)L$(Li)n$(nic_str)-$(nic_end)_output"
    local fn_spt = filename_spt3(spt)
    local fn = string(FN_INPUTDIR2,"outputs/",fn_spt,"n$(nic_str)-$(nic_end)_output")
    local tms = zeros(Int, 3, 2)
    open(fn) do fo
        for line in eachline(fo)
            local t, tm
            if length(line) > 3 && line[1:3] == "rea"
                t = split(line[6:end-2], "m")
                if parse(Float64, t[2]) > 30
                    tm = parse(Int, t[1])+1
                else
                    tm = parse(Int, t[1])
                end
                tms[1, 1] = div(tm, 60)
                tms[1, 2] = tm%60
            elseif length(line) > 3 && line[1:3] == "use"
                t = split(line[6:end-2], "m")
                if parse(Float64, t[2]) > 30
                    tm = parse(Int, t[1])+1
                else
                    tm = parse(Int, t[1])
                end
                tms[2, 1] = div(tm, 60)
                tms[2, 2] = tm%60
            elseif length(line) > 3 && line[1:3] == "sys"
                t = split(line[5:end-2], "m")
                if parse(Float64, t[2]) > 30
                    tm = parse(Int, t[1])+1
                else
                    tm = parse(Int, t[1])
                end
                tms[3, 1] = div(tm, 60)
                tms[3, 2] = tm%60
            end
        end
    end
    return tms
end

###
###  FUNCTION: SUSY Ward--Takahashi Identities (One-point)
###
#  \delta defined by Eq.~(5.8) in [Kamata--Suzuki]
#  with Standard Deviation
function wt_onepoint_delta(nic::Array{Complex{Float64}, 3}, sgn::Array{Int64, 2})
    local const N = size(nic, 1)
    local const L = Int(floor(sqrt(N)))
    local sb_nm = [sum(abs2(nic[:, i_f, i])) * sum(sgn[:, i]) / (L-1)^2 for i = 1:size(nic, 3), i_f = 1:size(nic, 2)]
    local sb_nv = sum(sb_nm, 2)
    local sb_n = mean(sb_nv)
    local Delta, Del_dev
    Delta, Del_dev = witten_ind(sgn)
    local sb_exp = sb_n / Delta
    local sb_err = abs(sb_exp) * sqrt(mean(abs2(sb_nv - sb_n)) / abs2(sb_n) + abs2(Del_dev/Delta))
    return sb_exp / (N*size(nic, 2)) - 1, ( sb_err / sqrt(size(nic, 3)-1) ) / (N*size(nic,2))
end

###
###  FUNCTION: SUSY Ward--Takahashi Identities (Two-point)
###
#  \langle \psi_1(p) \Bar{\psi}_{\Dot{1}}(-p) \rangle
#  with Standard Deviation
function wt_pp(phi::Array{Complex{Float64}, 3}, sgn::Array{Int64, 2}, k::Int, lm::Float64)
    local const N = size(phi, 1)
    local const L = Int(floor(sqrt(N)))
    local p = p_matrix2(N) * (L-1)^2
    local J = zeros(Complex{Float64}, 2N, 2N)
    J[1:N, 1:N] = p
    J[(N+1):2N, (N+1):2N] = - conj(p)
    local pp = zeros(Complex{Float64}, N, size(phi, 3))
    for i = 1:size(phi, 3)
        local tmp = zeros(Complex{Float64}, N)#
        local rs = zeros(Complex{Float64}, N)#
        for j = 1:size(phi, 2)
            if (norm(phi[:, j, i]) == 0)
                continue
            else
                local w2 = fieldtomatrix(superpotential_2(phi[:, j, i], k, lm))
                J[(N+1):2N, 1:N] = w2
                J[1:N, (N+1):2N] = w2'
                local Jinv = inv(J)
                #pp[:, i] += diag(Jinv)[1:N] * sgn[j, i]
                local u = rs + diag(Jinv)[1:N] * sgn[j, i]#
                local s = tmp + u#
                local t = s - tmp#
                rs = u - t#
                tmp = s#
            end
        end
        pp[:, i] = tmp#
    end
    local pp_sum, pp_dev
    pp_sum, pp_dev = wt_dev(pp)
    #local Delta = sum(sgn) / size(phi, 3)
    local Delta, Del_dev
    Delta, Del_dev = witten_ind(sgn)
    local res = (L-1)^4 * pp_sum / Delta
    local res_dev = abs(imag(res))im .* sqrt(abs2(imag(pp_dev)./imag(pp_sum)) + abs2(Del_dev/Delta))
    res_dev += abs(real(res)) .* sqrt(abs2(real(pp_dev)./real(pp_sum)) + abs2(Del_dev/Delta))
    return res, res_dev
end
#  \langle A^*(-p) A(p) \rangle
#  with Standard Deviation
function wt_aa(phi::Array{Complex{Float64}, 3}, sgn::Array{Int64, 2})
    local aa = zeros(Float64, size(phi, 1), size(phi, 3))
    for i = 1:size(phi, 3)
        aa[:, i] = abs2(phi[:, :, i]) * sgn[:, i]
    end
    local aa_sum, aa_dev
    aa_sum, aa_dev = wt_dev(aa)
    #local Delta = sum(sgn) / size(phi, 3)
    local Delta, Del_dev
    Delta, Del_dev = witten_ind(sgn)
    local res = aa_sum / Delta
    return res, abs(res) .* sqrt(abs2(aa_dev./aa_sum) + abs2(Del_dev/Delta))
end
#  \langle F^*(-p) F(p) \rangle
#  with Standard Deviation
function wt_ff(phi::Array{Complex{Float64}, 3}, nic::Array{Complex{Float64}, 2}, sgn::Array{Int64, 2})
    local const N = size(phi, 1)
    local const L = Int(floor(sqrt(N)))
    local pz = p_matrix(N)
    local ff = zeros(Float64, size(phi, 1), size(phi, 3))
    for i = 1:size(phi, 3)
        ff[:, i] = [abs2( nic[p, i] - pz[p]*phi[p, j, i] ) for p=1:N, j=1:size(phi, 2)] * sgn[:, i]
    end
    local ff_sum, ff_dev
    #local Delta = sum(sgn) / size(phi, 3)
    local Delta, Del_dev
    Delta, Del_dev = witten_ind(sgn)
    ff_sum, ff_dev = wt_dev(ff - Delta * (L-1)^2)
    local res = ff_sum / Delta
    return res, abs(res) .* sqrt(abs2(ff_dev./ff_sum) + abs2(Del_dev/Delta))
end

###
###  FUNCTION: Output Correlation Functions
###
#  For SUSY Ward--Takahashi identities
#  table: p_0, real(psi-psi&err), imag(~), real(2ip(psi-psi)&err), imag(~), real(2ip(A-A)&err), imag(~), F-F&err
#  blocks with fixed p_1
function wt_output(pp::Array{Complex{Float64}, 1}, pp_dev::Array{Complex{Float64}, 1}, aa::Array{Float64, 1}, aa_dev::Array{Float64, 1}, ff::Array{Float64, 1}, ff_dev::Array{Float64, 1}, Li::Int, k::Int, lm::Float64)
    local pm = p_matrix((Li + 1)^2)
    local p_pp, p_pp_dev, p_aa, p_aa_dev, p2_aa, p2_aa_dev
    p_pp = pm .* pp
    p_pp_dev = sqrt( abs2(real(pm).*real(pp_dev)) + abs2(imag(pm).*imag(pp_dev)) ) + sqrt( abs2(real(pm).*imag(pp_dev)) + abs2(imag(pm).*real(pp_dev)) )im
    p_aa = -conj(pm) .* aa
    p_aa_dev = (real(pm) .* aa_dev) + (imag(pm) .* aa_dev)im
    p2_aa = abs2(pm) .* aa
    p2_aa_dev = abs2(pm) .* aa_dev
    local num_col = 17
    local corr = zeros(Float64, (Li + 1)^2, num_col)
    corr[:, 1]  = imag(pm)
    corr[:, 2]  = real(pp)
    corr[:, 3]  = real(pp_dev)
    corr[:, 4]  = imag(pp)
    corr[:, 5]  = imag(pp_dev)
    corr[:, 6]  = real(p_pp)
    corr[:, 7]  = abs(real(p_pp_dev))
    corr[:, 8]  = imag(p_pp)
    corr[:, 9]  = abs(imag(p_pp_dev))
    corr[:, 10] = real(p_aa)
    corr[:, 11] = abs(real(p_aa_dev))
    corr[:, 12] = imag(p_aa)
    corr[:, 13] = abs(imag(p_aa_dev))
    corr[:, 14] = ff
    corr[:, 15] = ff_dev
    corr[:, 16] = p2_aa
    corr[:, 17] = p2_aa_dev
    fout = string(FN_OUTPUTDIR, "wt/wt_lm$(lm)k$(k)L$(Li)_p1fix.dat")
    open(fout, "w") do fo
        for p1 = 1:(Li+1)
            writedlm(fo, [corr[(Li+1) * i + p1, j] for i=0:Li, j=1:num_col])
            write(fo, "\n\n")
        end
    end
    fout0 = string(FN_OUTPUTDIR, "wt/wt_lm$(lm)k$(k)L$(Li)_p0fix.dat")
    corr[:, 1]  = real(pm)
    open(fout0, "w") do fo
        for p0 = 0:Li
            writedlm(fo, [corr[i + (Li+1) * p0, j] for i=1:(Li+1), j=1:num_col])
            write(fo, "\n\n")
        end
    end
    return 0
end

###
###  MAIN FUNCTION: Calculation process
###
#  Ward--Takahashi identities
function main_wt(spt::Spt_type, ft)
    #  Read Configuraitons
    local conf = read_conf(spt)
    if conf==0
        return 0
    end
    local num_nic = size(conf.sgn, 2)
    local num_sol = size(conf.sgn, 1)
    local cls = read_cls(conf.sgn)
    
    #  Witten index
    local wit, wit_dev
    wit, wit_dev = witten_ind(conf.sgn)
    
    #  SUSY WT Idnetities
    ##  Onepoint
    local sb, sb_dev
    sb, sb_dev = wt_onepoint_delta(conf.nic, conf.sgn)
    
    ##  Twopoint
    ###  Calc \langle hoge(p) hoge(-p) \rangle
    #local pp, pp_dev, aa, aa_dev, ff, ff_dev
    #pp, pp_dev = wt_pp(phis, sgns, k, lm)
    #aa, aa_dev = wt_aa(phis, sgns)
    #ff, ff_dev = wt_ff(phis, nics, sgns)
    ###  Output correlation functions
    #wt_output(pp, pp_dev, aa, aa_dev, ff, ff_dev, Li, k, lm)
    
    #  Time
    local tm = zeros(Int, 3, 2)
    
    local num_nic_d = set_num_nic_d(spt)
    for i in 0:num_nic_d:(num_nic-num_nic_d)
        tm += read_out(spt, i, i+num_nic_d)
    end
    for i = 1:3
        if tm[i, 2] >= 60
            local tm_h = div(tm[i, 2], 60)
            tm[i, 1] += tm_h
            tm[i, 2] -= tm_h * 60
        end
    end
    
    local fn_spt = filename_spt2(spt)
    fout = string(FN_OUTPUTDIR, "wt_", fn_spt, "_output")
    open(fout, ft) do fo
        if spt.spt==0
            write(fo, "# SuperPotentialType\nAlgebra A\n")
        elseif spt.spt==1
            write(fo, "# SuperPotentialType\nAlgebra D\n")
        elseif spt.spt==2
            write(fo, "# SuperPotentialType\nAlgebra E\n")
        elseif spt.spt==3
            write(fo, "# SuperPotentialType\nCustom\n")
        end
        write(fo, "# Physical Box Size, Power, Coupling\n")
        write(fo, "$(spt.Li), $(spt.k), $(spt.lm)\n")
        write(fo, "# Set Number of N-Configurations and its Scalar Solutions\n")
        write(fo, "$(num_nic), $(num_sol)\n")
        
        write(fo, "# Max Residue\n")
        write(fo, "$(conf.max_res)\n")
        write(fo, "# Number of Dumping (Max, Average)\n")
        write(fo, "$(maximum(conf.dmp)), $(mean(conf.dmp))\n")
        write(fo, "# Classification of confs\n")
        for i = 1:size(cls, 2)
            write(fo, "(")
            for j = 1:div(cls[1, i] + cls[2, i], 2)
                write(fo, "{+}")
            end
            for j = 1:div(cls[1, i] - cls[2, i], 2)
                write(fo, "{-}")
            end
            write(fo, ")_$(cls[2, i]) :  $(cls[3, i])\n")
        end
        
        write(fo, "# Witten index (Delta)\n")
        write(fo, "$(wit)\n$(wit_dev)\n")
        
        write(fo, "# WT Onepoint (delta)\n")
        write(fo, "$(sb)\n$(sb_dev)\n")
        
        #write(fo, "# WT Twopoint\n")
        #write(fo, "# Max Standard Deviation of (psi-psi, A-A, F-F) / norm\n")
        #write(fo, "$(maxabs(pp_dev) / norm(pp))\n")
        #write(fo, "$(maximum(aa_dev) / norm(aa))\n")
        #write(fo, "$(maximum(ff_dev) / norm(ff))\n")
        
        write(fo, "# Time[h:m] (real, user, sys)\n")
        for i = 1:3
            if tm[i, 2] < 10
                write(fo, "$(tm[i, 1]):0$(tm[i, 2])\n")
            else
                write(fo, "$(tm[i, 1]):$(tm[i, 2])\n")
            end
        end
        write(fo, "\n")
    end
    return 0
end
function main_wt(num::Int, strn::Int, endn::Int)
    local ft
    for i=strn:endn
        if i==1
            ft = "w"
        else
            ft = "a"
        end
        local spt = set_spt(i,num)
        if spt==0
            return 0
        end
        main_wt(spt, ft)
    end
    return 0
end
##################################################

#####################################################################
### FUNCTION for Tables (TeX) of Classification of Configurations ###
#####################################################################
function read_wtout(spt::Spt_type)
    local num_Ls = 20
    local tbl = Array{String, 2}(7, num_Ls+1)
    local row_init = 15
    local row = 0
    local cls     = Array{String, 1}(row_init)
    local cls_num = zeros(Int, row_init, num_Ls+1)
    local Lindex = 0
    local fn_spt = filename_spt2(spt)
    local fn = string(FN_OUTPUTDIR, "wt_", fn_spt, "_output")
    open(fn) do fi
        for line in eachline(fi)
            if line=="# Physical Box Size, Power, Coupling\n"
                local par = replace(readline(fi), " ", "")
                local pars = split(par, ",")
                Lindex += 1
                tbl[1,Lindex] = pars[1]
            elseif line=="# Max Residue\n"
                local res = readline(fi)
                local resn = round(parse(Float64, res)*(10^15), 6)
                tbl[2, Lindex] =  string(resn)
            elseif line=="# Classification of confs\n"
                while(true)
                    local tmp = replace(readline(fi), " ", "")
                    if tmp=="#Wittenindex(Delta)\n"
                        break
                    end
                    local tmps = split(tmp, ":")
                    for i = 1:(row+1)
                        if i == (row+1)
                            cls[i] = tmps[1]
                            cls_num[i, Lindex] = parse(Int, tmps[2])
                            row += 1
                        elseif cls[i] == tmps[1]
                            cls_num[i, Lindex] = parse(Int, tmps[2])
                            break
                        end
                    end
                end
            #elseif line=="# Witten index (Delta)\n"
                local wtt     = parse(Float64, readline(fi))
                local wtt_dev = parse(Float64, readline(fi))
                local ktt = 3
                local t
                t  = string(round(wtt, ktt))
                for i = 1:(2+ktt-length(t))
                    t = string(t, "0")
                end
                if wtt_dev == 0
                    tbl[3, Lindex] = "$(t)"
                else
                    local td
                    td = string(Int(round(wtt_dev * 10^ktt)))
                    tbl[3, Lindex] = "$(t)($(td))"
                end
            elseif line=="# WT Onepoint (delta)\n"
                local wt     = parse(Float64, readline(fi))
                local wt_dev = parse(Float64, readline(fi))
                local wt2 = abs(wt)
                local sign = (wt>0)?"":"\$-\$"
                local kt = 5
                local w  = string(round(wt2, kt))
                for i = 1:(2+kt-length(w))
                    w = string(w, "0")
                end
                local wd = string(Int(round(wt_dev * 10^kt)))
                tbl[4, Lindex] = "$(sign)$(w)($(wd))"
            elseif line=="# Time[h:m] (real, user, sys)\n"
                local tr = readline(fi)
                local tu = readline(fi)
                local ts = readline(fi)
                tbl[5, Lindex] = tr[1:end-1]
                tbl[6, Lindex] = tu[1:end-1]
                tbl[7, Lindex] = ts[1:end-1]
            end
        end
    end
    return tbl[:,1:Lindex], cls[1:row], cls_num[1:row, :]
end
function cls_output(tbl::Array{String, 2}, cls::Array{String, 1}, cls_num::Array{Int, 2}, Lmin::Int, Lmax::Int, spt::Spt_type)
    local fn_spt = filename_spt(spt)
    local fn = string("../notes/cls_", fn_spt, ".tex")
    local ft
    if (parse(Int, tbl[1, 1]) == Lmin)
        ft = "w"
    else
        ft = "a"
    end
    open(fn, ft) do fo
        write(fo, "\\begin{table}[ht]\n")
        write(fo, " \\centering\n")
        write(fo, " \\begin{tabular}{l")
        for i = 1:size(tbl, 2)
            write(fo, "r")
        end
        write(fo, "}\\toprule\n")
        write(fo, "  \$L\$         ")
        for i = 1:size(tbl, 2)
            write(fo, "& $(tbl[1, i]) ")
        end
        write(fo, "\\\\\\midrule\n")
        for i = 1:length(cls)
            write(fo, "  \$$(cls[i])\$ ")
            for j = 1:size(tbl, 2)
                write(fo, "& $(cls_num[i, j]) ")
            end
            write(fo, "\\\\\n")
        end
        write(fo, "  \$\\Delta\$ ")
        for i = 1:size(tbl, 2)
            write(fo, "& $(tbl[3, i]) ")
        end
        write(fo, "\\\\\n")
        write(fo, "  \$\\delta\$ ")
        for i = 1:size(tbl, 2)
            write(fo, "& $(tbl[4, i]) ")
        end
        write(fo, "\\\\\n  \\midrule\n")
        write(fo, "  %Residue(e-15) ")
        for i = 1:size(tbl, 2)
            write(fo, "& $(tbl[2, i]) ")
        end
        write(fo, "\\\\\n")
        write(fo, "  %real [h:m] ")
        for i = 1:size(tbl, 2)
            local rtime = parse(Int, tbl[5, i][1:end-3])+parse(Int, tbl[5,i][end-1:end])/60
            rtime = round(rtime*10^2)/10^2
            write(fo, "& $(rtime) ")
        end
        write(fo, "\\\\\n")
        write(fo, "  core\$\\cdot\$hour [h] ")
        for i = 1:size(tbl, 2)
            local utime = parse(Int, tbl[6, i][1:end-3])+parse(Int, tbl[6,i][end-1:end])/60
            utime = round(utime*10^2)/10^2
            write(fo, "& $(utime) ")
        end
        write(fo, "\\\\\n")
        write(fo, "  %sys [h:m] ")
        for i = 1:size(tbl, 2)
            local stime = parse(Int, tbl[7, i][1:end-3])+parse(Int, tbl[7,i][end-1:end])/60
            stime = round(stime*10^2)/10^2
            write(fo, "& $(stime) ")
        end
        write(fo, "\\\\\n  \\bottomrule\n")
        write(fo, " \\end{tabular}\n")
        write(fo, " \\caption{Classification of configurations")
        if spt.spt==0
            write(fo, " for~\$A_{$(spt.k[1]-1)}\$")
        elseif spt.spt==1
            write(fo, " for~\$D_{$(spt.k[1]+1)}\$")
        elseif spt.spt==2
            write(fo, " for~\$E_{7}\$")
        elseif spt.spt==2
            write(fo, "")
        end
        if (ft=="w")
            write(fo, "}\n")
        else
            write(fo, " (continued)}\n")
        end
        write(fo, " \\label{tab:class.")
        write(fo, fn_spt)
        if     (ft == "w")
            write(fo, "_1}\n")
        elseif (parse(Int, tbl[1, end]) == Lmax)
            write(fo, "_2}\n")
        end
        write(fo, "\\end{table}%\n")
    end
    return 0
end
function main_cls(num::Int)
    local spt = set_spt(1, num)
    local tbl, cls, cls_num
    tbl, cls, cls_num = read_wtout(spt)
    local Lmin = parse(Int, tbl[1, 1])
    local Lmax = parse(Int, tbl[1,end])
    local stp = 4
    for i = 1:stp:size(tbl,2)
        local j = i + stp - 1
        if j > size(tbl,2)
            j = size(tbl,2)
        end
        cls_output(tbl[:, i:j], cls, cls_num[:, i:j], Lmin, Lmax, spt)
    end
    return 0
end
#####################################################################



#spt = Spt_type(8, 3, [4], [0.3;0.4;0.5;0.6], [1;1;1;1])
#conf = read_conf(spt)

