function Create_Stoch_Pert(nx_length, ny_length)
    pert_spectrum = zeros(ComplexF64, nx_length, ny_length)
    for j = 1:Int(ny_length)
        for i = 1:Int(nx_length)
            if Init.stoch_spec[i,j] >= 0
                real_value = rand(Normal(),1);
                complex_value = rand(Normal(),1);
                ## TODO Look back at square root of spectrum to create a flat spectrum
                pert_spectrum[i,j] = (Init.stoch_spec[i,j]*(real_value + complex_value*im))[1];
            end
        end
    end
    # Inverse FFT to move to physical space
    pert_phys = FFTW.irfft(pert_spectrum, Int(Init.parameters.N_points))
    
    #Scaling physically
    phys = pert_phys.*spatial_localize()
    
    phys_spec = FFTW.rfft(phys)
    
    scaled_spec = -(Init.k2).*phys_spec
    scaled_spec[1,1] = 0 +0im
    
    return scaled_spec
    
end

function spatial_localize()
    MR = ones(Init.parameters.N_points, Int(Init.parameters.N_points))
    for i = 1:Int(Init.parameters.N_points)
        for j = 1:Init.parameters.N_points
            if Init.y[i,j] < -Init.Ly/4
                MR[i,j] = 1 - exp(0.1/((2*Init.y[i,j]/Init.Ly+1)^2 - 0.25))
            end
            if Init.y[i,j] > Init.Ly/4  
                MR[i,j] = 1 - exp(0.1/((2*Init.y[i,j]/Init.Ly-1)^2 - 0.25))
            end
        end
    end
    MR = Init.MR_scale*MR
    return MR
end