function Setup_Psi()  
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:Init.parameters.N_points/2); collect(-Init.parameters.N_points/2+1:-1)]
    
    if(Init.parameters.model_type == 1)
        KX = Make_2d_Array(1, k);
        KY = Make_2d_Array(2, k);
        Laplacian = -(KX .^ 2 + KY .^ 2);
        Init.psi_param.InvBT = 1 ./ Laplacian; 
        Init.psi_param.InvBT[1,1] = 0;
        Init.psi_param.InvBC = 1 ./ (Laplacian .- Init.parameters.k_deform^2);
        Init.psi_param.InvBC[1,1] = 0;
    else
        KX = Make_2d_Array(1, k);
        KY = Make_2d_Array(2, k);
        K2 = KX .^ 2 + KY .^ 2;
        KS = sqrt.(K2);
        Init.psi_param.InvMat = zeros(Init.parameters.N_points, Init.parameters.N_points, 2, 2);
        Init.psi_param.InvMat[:,:,1,1] = 1 ./ (KS .* tanh.(KS));
        Init.psi_param.InvMat[:,:,1,2] = -1 ./ (KS .* sinh.(KS));
        Init.psi_param.InvMat[:,:,2,1] = 1 ./ (KS .* sinh.(KS));
        Init.psi_param.InvMat[:,:,2,2] = -1 ./ (KS .* tanh.(KS));
        Init.psi_param.InvMat[1,1,:,:] = 0;
    end
    
end

function Get_Psi(q_hat)
    if(Init.parameters.model_type == 1)
        # Invert for psi
        q_bt = .5*(q_hat[:,:,1] + q_hat[:,:,2]);
        q_bc = .5*(q_hat[:,:,1] - q_hat[:,:,2]);
        psi_bt = InvBT .* q_bt;
        psi_bc = InvBC .* q_bc;
        psi_hat[:,:,2] = psi_bt-psi_bc;
        psi_hat[:,:,1] = psi_bt+psi_bc;
    else
        #Invert for psi
        psi_hat[:,:,1] = InvMat[:,:,1,1] .* q_hat[:,:,1] + InvMat[:,:,1,2] .* q_hat[:,:,2];
        psi_hat[:,:,2] = InvMat[:,:,2,1] .* q_hat[:,:,1] + InvMat[:,:,2,2] .* q_hat[:,:,2];
    end
    return psi_hat
end

function Make_2d_Array(axis::Integer, reference_array)
    # Axis = 1 is x changing
    # Axis != 1 is y changing 
    a_size = size(reference_array)[1]
    return_array = zeros(a_size, a_size)
    for i in 1:a_size
        if axis == 1
            return_array[i,:] = reference_array
        else
            return_array[:,i] = transpose(reference_array)
        end
    end
    return return_array
end
            
function Make_2d_Array(axis::Integer, reference_array, size)
    # Axis = 1 is x changing
    # Axis != 1 is y changing 
    return_array = zeros(Int64(size), Int64(size))
    for i in 1:Int64(size)
        if axis == 1
            return_array[i,:] = reference_array
        else
            return_array[:,i] = transpose(reference_array)
        end
    end
    return return_array
end           
            
            
            