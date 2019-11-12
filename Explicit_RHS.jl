# Function to setup a struct of parameters used in the QG RHS function
function setup_QG_RHS()
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:Init.parameters.N_points/2); collect(-Init.parameters.N_points/2+1:-1)]
    dX = 1im*Make_2d_Array(1, k[Int(1):div(Init.parameters.N_points,2)+1], Init.parameters.N_points)
    dY = 1im*Make_2d_Array(2, k, div(Init.parameters.N_points,2)+1)
    Init.laplace_param.Laplacian = dX.^2+dY.^2;
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:Init.parameters.N_points/2-1); 0; collect(-Init.parameters.N_points/2+1:-1)];
    Init.laplace_param.dX = 1im*Make_2d_Array(1, k[Int(1):div(Init.parameters.N_points,2)+1], Init.parameters.N_points);
    Init.laplace_param.dY = 1im*Make_2d_Array(2, k, div(Init.parameters.N_points,2)+1);
    # For the dealiased jacobian:
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:.75*Init.parameters.N_points-1); 0; collect(-.75*Init.parameters.N_points+1:-1)];
    DX = 1im*Make_2d_Array(1, k[1:Int(div(1.5*Init.parameters.N_points,2)+1)], Int(1.5*Init.parameters.N_points));
    Init.laplace_param.DX = zeros(ComplexF64, size(DX)[1],size(DX)[2],2);
    Init.laplace_param.DX[:,:,1] = DX;
    Init.laplace_param.DX[:,:,2] = DX;
    DY = 1im*Make_2d_Array(2, k, Int(div(1.5*Init.parameters.N_points,2)+1));
    Init.laplace_param.DY = zeros(ComplexF64, size(DY)[1],size(DY)[2],2);
    Init.laplace_param.DY[:,:,1] = DY;
    Init.laplace_param.DY[:,:,2] = DY;
end

function QG_RHS(q_hat)
    
l_params = Init.laplace_param;
params = Init.parameters;    
# Function takes Fourier coefficients of PV (q_hat) and struct containing
# parameters (p) and evaluates RHS of 2-layer QG equations except for
# high-k dissipation. Returns Fourier coefficients of RHS.
# Jacobian is dealiased via 3/2 rule.
#########################################################################

# Invert for psi
psi_hat = Get_Psi(q_hat);

RHS = zeros(ComplexF64, div(params.N_points,2)+1, params.N_points, 2);
if( params.model_type==1  && params.k_beta != 0)
    # beta plus mean shear
    RHS[:,:,1] =  params.k_beta^2 * dX .* psi_hat[:,:,1];
    RHS[:,:,2] = -params.k_beta^2 * dX .* psi_hat[:,:,2];
else
    # mean buoyancy gradient plus mean shear
    RHS[:,:,1] =-l_params.dX .* q_hat[:,:,1] + 2 * l_params.dX .* psi_hat[:,:,1];
    RHS[:,:,2] = l_params.dX .* q_hat[:,:,2] + 2 * l_params.dX .* psi_hat[:,:,2];
end

# Bottom drag
#u = -FFTW.irfft(dY.*psi_hat[:,:,2]);
#v = FFTW.irfft(dX.*psi_hat[:,:,2]);
#U = sqrt(u.^2+v.^2);
#RHS[:,:,2] = RHS[:,:,2] - params.cd*l_params.dX.*FFTW.rfft(U.*v)+params.cd*l_params.dY.*FFTW.rfft(U.*u);
RHS[:,:,2] = RHS[:,:,2] - params.r_ekman * l_params.Laplacian .* psi_hat[:,:,2];

# The following code computes the Jacobian J[psi,q], dealiased using a 3/2-rule.
# I.e., if you set N=256 then the code uses 3*N/2=384 Fourier modes (in each
# direction) to compute the Jacobian.
# physical space, 3/2 grid; factor of (9/4) scales fft 
# Note: The FFTW uses only 1/2 the size of the Kx direction so this is not a square matrix setup
Psi_hat = zeros(ComplexF64, Int(div(1.5*params.N_points,2)+1), Int(1.5*params.N_points), 2);
    
# Make upper left piece
Psi_hat[1:div(params.N_points,2)+1,1:div(params.N_points,2)+1,:] = 
    (9/4) * psi_hat[1:div(params.N_points,2)+1, 1:div(params.N_points,2)+1, :];
    
# Make the upper right piece
Psi_hat[1:div(params.N_points,2)+1,params.N_points+2:Int(1.5*params.N_points),:] = (9/4)*psi_hat[1:div(params.N_points,2)+1,div(params.N_points,2)+2:params.N_points,:];
    
## Make lower left piece
#Psi_hat[params.N_points+2:Int(div(1.5*params.N_points,2)+1),1:div(params.N_points,2)+1,:] = #(9/4)*psi_hat[div(params.N_points,2)+2:params.N_points,1:div(params.N_points,2)+1,:];
    
# Make lower right piece
#Psi_hat[params.N_points+2:1.5*params.N_points,params.N+2:1.5*params.N_points,:] = (9/4)*psi_hat[div(params.N_points,2)+2:params.N_points,div(params.N_points,2)+2:params.N_points,:];
    
Q_hat = zeros(ComplexF64, Int(div(1.5*params.N_points,2)+1), Int(1.5*params.N_points), 2);
    
Q_hat[1:div(params.N_points,2)+1,1:div(params.N_points,2)+1,:] = (9/4)*q_hat[1:div(params.N_points,2)+1,1:div(params.N_points,2)+1,:];
    
Q_hat[1:div(params.N_points,2)+1,params.N_points+2:Int(1.5*params.N_points),:] = (9/4)*q_hat[1:div(params.N_points,2)+1,div(params.N_points,2)+2:params.N_points,:];
    
#Q_hat[params.N_points+2:1.5*params.N_points,1:div(params.N_points,2)+1,:] = (9/4)*q_hat[div(params.N_points,2)+2:params.N_points,1:div(params.N_points,2)+1,:];
    
#Q_hat[params.N_points+2:1.5*params.N_points,params.N+2:1.5*params.N_points,:] = (9/4)*q_hat[div(params.N_points,2)+2:params.N_points,div(params.N_points,2)+2:params.N_points,:];
    
# calculate u.gradq on 3/2 grid
u = FFTW.irfft( -l_params.DY.*Psi_hat, Int(params.N_points*1.5));
v = FFTW.irfft(  l_params.DX.*Psi_hat, Int(params.N_points*1.5));
qx= FFTW.irfft(  l_params.DX.*Q_hat, Int(params.N_points*1.5));
qy= FFTW.irfft(  l_params.DY.*Q_hat, Int(params.N_points*1.5));
jaco_real = u.*qx+v.*qy;
# fft, 3/2 grid; factor of (4/9) scales fft
Jaco_hat = (4/9)*FFTW.rfft(jaco_real);
# reduce to normal grid
jaco_hat = zeros(ComplexF64,div(params.N_points,2)+1, params.N_points, 2);
#Make top left
jaco_hat[1:div(params.N_points,2)+1,1:div(params.N_points,2)+1,:] = Jaco_hat[1:div(params.N_points,2)+1,1:div(params.N_points,2)+1,:];
    
# Make top right
jaco_hat[1:div(params.N_points,2)+1,div(params.N_points,2)+2:params.N_points,:] = Jaco_hat[1:div(params.N_points,2)+1,params.N_points+2:Int(1.5*params.N_points),:];

#Make bottom left
#jaco_hat[div(params.N_points,2)+2:params.N_points,1:div(params.N_points,2)+1,:] = Jaco_hat[params.N_points+2:Int(1.5*params.N_points),1:div(params.N_points,2)+1,:];
    
# Make bottom right
#jaco_hat[div(params.N_points,2)+2:params.N_points,div(params.N_points,2)+2:params.N_points,:] = Jaco_hat[params.N+2:Int(1.5*params.N_points),params.N_points+2:Int(1.5*params.N_points),:];
    
# Put it all together
RHS = RHS - jaco_hat;
        
end