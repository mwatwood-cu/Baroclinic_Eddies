function take_timestep_array(rhs::Function, L, u, dt)
    #u should be an array containing the current values of the desired quantities
    #rhs should be a function to operate on the explicit part of the differential equation
    #L should be a matrix representation of the implicit operator used
    #dt is the step size to take
    
    
    #Take a time step using RK stages for ARk[4](3)6L[2]SA 
    #see citation [FIX_ME] in the decription in the main notebook
    ae, b, ai, be  =get_rk_constants()
    
    
    M = 1 ./ (1 .- (0.25*dt).*L)
    
    
    #k is the derivative estimate of the explicit piece
    k = zeros(ComplexF64,size(u)[1], size(u)[2], size(u)[3], 6);
    #l is the derivative estimate of the implicit piece
    l= zeros(ComplexF64, size(u)[1], size(u)[2], size(u)[3], 6);
    
    #Stage 1
    k[:,:,:,1] = rhs(u)
    l[:,:,:,1] = L .* u
    
    #Stage 2
    u_tmp = M .* (u+dt*(ae[2,1]*k[:,:,:,1] + ai[2,1]*l[:,:,:,1]))
    k[:,:,:,2] = rhs(u_tmp)
    l[:,:,:,2] = L .* u_tmp
    
    #Stage 3
    u_tmp = M .* (u+dt*(ae[3,1]*k[:,:,:,1] + ae[3,2]*k[:,:,:,2] 
                       + ai[3,1]*l[:,:,:,2] + ai[3,2]*l[:,:,:,2]))
    k[:,:,:,3] = rhs(u_tmp) 
    l[:,:,:,3] = L .* u_tmp
    
    #Stage 4
    u_tmp = M .* (u+dt*(ae[4,1]*k[:,:,:,1] + ae[4,2]*k[:,:,:,2] + ae[4,3]*k[:,:,:,3] 
                       + ai[4,1]*l[:,:,:,1] + ai[4,2]*l[:,:,:,2] + ai[4,3]*l[:,:,:,3]))
    k[:,:,:,4] = rhs(u_tmp) 
    l[:,:,:,4] = L .* u_tmp
    
    #Stage 5
    u_tmp = M .* (u+dt*(ae[5,1]*k[:,:,:,1] + ae[5,2]*k[:,:,:,2]  + ae[5,3]*k[:,:,:,3] + ae[5,4]*k[:,:,:,4]
                       + ai[5,1]*l[:,:,:,1] + ai[5,2]*l[:,:,:,2]  + ai[5,3]*l[:,:,:,3] + ai[5,4]*l[:,:,:,4]))
    
    k[:,:,:,5] = rhs(u_tmp) 
    l[:,:,:,5] = L .* u_tmp
    
    #Stage 6
    u_tmp = M .*(u+dt*(ae[6,1]*k[:,:,:,1] + ae[6,2]*k[:,:,:,2]  + ae[6,3]*k[:,:,:,3] + ae[6,4]*k[:,:,:,4] + ae[6,5]*k[:,:,:,5]
                       + ai[6,1]*l[:,:,:,1] + ai[6,2]*l[:,:,:,2]  + ai[6,3]*l[:,:,:,3] + ai[6,4]*l[:,:,:,4] + ai[6,5]*l[:,:,:,5]))
    k[:,:,:,6] = rhs(u_tmp) 
    l[:,:,:,6] = L .* u_tmp
    
    #Error Check
    fourier_error_sum = zeros(ComplexF64,size(u)[1], size(u)[2], size(u)[3]) 
    for i=1:6
        fourier_error_sum = fourier_error_sum + be[i]*(k[:,:,:,i]+l[:,:,:,i])
    end
    real_error_sum = FFTW.irfft(fourier_error_sum, Init.parameters.N_points);
    error_max = dt*maximum(real_error_sum)
    
    #println(Init.isAdaptive)
    if(Init.isAdaptive && error_max > Init.tol )
        return take_timestep_array(rhs, L, u, 0.75*dt)
    end
    
    #On a successful round
    u_next = u
    for i=1:6
        u_next = u_next + dt*(b[i]*k[:,:,:,i]+l[:,:,:,i])
    end
 
    if(Init.isAdaptive)
        r0 = Init.time_parameters.r0;
        next_dt = dt*((.75*Init.tol/error_max)^(.3/4))*((Init.r0/error_max)^(.4/4))
        Init.time_parameters.r0 = error_max
    else
        next_dt = dt
    end
    
    #Test case 1 same forcing for both top and bottom
    #stocastic = Create_Stoch_Pert(Int(Init.parameters.N_points/2+1), Init.parameters.N_points);
    #u_next[:,:,1] = u_next[:,:,1] + sqrt(next_dt).*stocastic;
    #u_next[:,:,2] = u_next[:,:,2] + sqrt(next_dt).*stocastic;
    
    
    return u_next, next_dt
    
end

function get_rk_constants()
    a_explicit = zeros(Float64, 6,6)
    a_explicit[2,1] = 0.5
    a_explicit[3,1] = 13861/62500
    a_explicit[3,2] = 6889/62500
    a_explicit[4,1] = -116923316275/2393684061468
    a_explicit[4,2] = -2731218467317/15368042101831
    a_explicit[4,3] = 9408046702089/11113171139209
    a_explicit[5,1] = -451086348788/2902428689909
    a_explicit[5,2] = -2682348792572/7519795681897
    a_explicit[5,3] = 12662868775082/11960479115383
    a_explicit[5,4] = 3355817975965/11060851509271
    a_explicit[6,1] = 647845179188/3216320057751
    a_explicit[6,2] = 73281519250/8382639484533
    a_explicit[6,3] = 552539513391/3454668386233
    a_explicit[6,4] = 3354512671639/8306763924573
    a_explicit[6,5] = 4040/17871

    b = zeros(Float64, 6)
    b[1] = 82889/524892
    b[2] = 0
    b[3] = 15625/83664
    b[4] = 69875/102672
    b[5] = -2260/8211
    b[6] = 0.25

    a_implicit = zeros(Float64, 6,6)
    a_implicit[2,1] = 0.25
    a_implicit[2,2] = 0.25 
    a_implicit[3,3] = 0.25
    a_implicit[4,4] = 0.25
    a_implicit[5,5] = 0.25
    # All diagonal elements except a_implicit[1,1] are 1/4
    
    a_implicit[3,1] = 8611/62500
    a_implicit[3,2] = -1743/31250
    a_implicit[4,1] = 5012029/34652500
    a_implicit[4,2] = -654441/2922500
    a_implicit[4,3] = 174375/388108
    a_implicit[5,1] = 15267082809/155376265600
    a_implicit[5,2] = -71443401/120774400
    a_implicit[5,3] = 730878875/902184768
    a_implicit[5,4] = 2285395/8070912
    a_implicit[6,:] = b

    b_error = zeros(Float64, 6)
    b_error[1] = 31666707/9881966720
    b_error[2] = 0
    b_error[3] = -256875/105007616
    b_error[4] = -2768025/128864768
    b_error[5] = 169839/3864644
    b_error[6] = -5247/225920
    
    return a_explicit, b, a_implicit, b_error

end
