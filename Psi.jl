function Setup_Psi()  
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:Init.parameters.N_points/2); collect(-Init.parameters.N_points/2+1:-1)]
    KX = Make_2d_Array(1, k[Int(1):div(Init.parameters.N_points,2)+1], Init.parameters.N_points);
    KY = Make_2d_Array(2, k, div(Init.parameters.N_points,2)+1);
    if(Init.parameters.model_type == 1)
        Laplacian = -(KX .^ 2 + KY .^ 2);
        Init.psi_param.InvBaroTropic = 1 ./ Laplacian; 
        Init.psi_param.InvBaroTropic[1,1] = 0;
        Init.psi_param.InvBaroClinic = 1 ./ (Laplacian .- Init.parameters.k_deform^2);
        Init.psi_param.InvBaroClinic[1,1] = 0;
    else
        K2 = KX .^ 2 + KY .^ 2;
        KS = sqrt.(K2);
        Init.psi_param.InvMat = zeros(Init.parameters.N_points/2+1, Init.parameters.N_points/2+1, 2, 2);
        Init.psi_param.InvMat[:,:,1,1] = 1 ./ (KS .* tanh.(KS));
        Init.psi_param.InvMat[:,:,1,2] = -1 ./ (KS .* sinh.(KS));
        Init.psi_param.InvMat[:,:,2,1] = 1 ./ (KS .* sinh.(KS));
        Init.psi_param.InvMat[:,:,2,2] = -1 ./ (KS .* tanh.(KS));
        Init.psi_param.InvMat[1,1,:,:] = 0;
    end
    
end

function Get_Psi(q_hat)
    psi_hat = zeros(ComplexF64, size(q_hat)[1],size(q_hat)[2],2)
    if(Init.parameters.model_type == 1)
        # Invert for psi
        q_bt = 0.5*(q_hat[:,:,1] .+ q_hat[:,:,2]);
        q_bc = 0.5*(q_hat[:,:,1] - q_hat[:,:,2]);
        psi_bt = Init.psi_param.InvBaroTropic .* q_bt;
        psi_bc = Init.psi_param.InvBaroClinic .* q_bc;
        psi_hat[:,:,2] = psi_bt .- psi_bc;
        psi_hat[:,:,1] = psi_bt .+ psi_bc;
    else
        #Invert for psi
        psi_hat[:,:,1] = Init.psi_param.InvMat[:,:,1,1] .* q_hat[:,:,1] + Init.psi_param.InvMat[:,:,1,2] .* q_hat[:,:,2];
        psi_hat[:,:,2] = Init.psi_param.InvMat[:,:,2,1] .* q_hat[:,:,1] + Init.psi_param.InvMat[:,:,2,2] .* q_hat[:,:,2];
    end
    return psi_hat
end

function Make_2d_Array(axis::Integer, reference_array, off_axis_length)
    # Axis = 1 is x 
    # Axis != 1 is y 
    a_size = size(reference_array)[1]
    if axis == 1
        return_array = zeros(a_size, off_axis_length)
        for i in 1:off_axis_length
            return_array[:,i] = transpose(reference_array);
        end
        return return_array
    else
        return_array = zeros(off_axis_length, a_size)
        for i in 1:off_axis_length
            return_array[i,:] = reference_array
        end
        return return_array
    end
end         
            
            
            