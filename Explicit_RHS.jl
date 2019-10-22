# Function to setup a struct of parameters used in the QG RHS function
function setup_QG_RHS()
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:Init.parameters.N_points/2); collect(-Init.parameters.N_points/2+1:-1)]
    dX = 1im*Make_2d_Array(1, k)
    dY = 1im*Make_2d_Array(2, k)
    Laplacian = dX.^2+dY.^2;
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:Init.parameters.N_points/2-1); 0; collect(-Init.parameters.N_points/2+1:-1)];
    Init.lap_param.dX = 1im*Make_2d_Array(1, k);
    Init.lap_param.dY = 1im*Make_2d_Array(2, k);
    # For the dealiased jacobian:
    k = (2*pi/Init.parameters.domain_wid)*[collect(0:.75*Init.parameters.N_points-1); 0; collect(-.75*Init.parameters.N_points+1:-1)];
    DX = 1im*Make_2d_Array(1, k, 1.5*Init.parameters.N_points);
    Init.lap_param.DX = zeros(ComplexF64, size(DX)[1],size(DX)[1],2);
    Init.lap_param.DX[:,:,1] = DX;
    Init.lap_param.DX[:,:,2] = DX;
    DY = 1im*Make_2d_Array(2, k, 1.5*Init.parameters.N_points);
    Init.lap_param.DY = zeros(ComplexF64, size(DY)[1],size(DY)[1],2);
    Init.lap_param.DY[:,:,1] = DY;
    Init.lap_param.DY[:,:,2] = DY;
end

function QG_RHS(q_hat, params, l_params)
    
# Function takes Fourier coefficients of PV (q_hat) and struct containing
# parameters (p) and evaluates RHS of 2-layer QG equations except for
# high-k dissipation. Returns Fourier coefficients of RHS.
# Jacobian is dealiased via 3/2 rule.
#########################################################################

# Invert for psi
psi_hat = GetPsi(q_hat, params);

RHS = zeros([params.N params.N 2]);
if( params.model==1 )
    # beta plus mean shear
    RHS[:,:,1] =-dX.*q_hat(:,:,1)-(params.kb^2 + params.kd^2)*dX.*psi_hat(:,:,1);
    RHS[:,:,2] = dX.*q_hat(:,:,2)-(params.kb^2 - params.kd^2)*dX.*psi_hat(:,:,2);
else
    # mean buoyancy gradient plus mean shear
    RHS[:,:,1] =-dX.*q_hat[:,:,1]+2*dX.*psi_hat[:,:,1];
    RHS[:,:,2] = dX.*q_hat[:,:,2]+2*dX.*psi_hat[:,:,2];
end

# Bottom drag
u = -FFTW.irfft(dY.*psi_hat[:,:,2]);
v = FFTW.irfft(dX.*psi_hat[:,:,2]);
U = sqrt(u.^2+v.^2);
RHS[:,:,2] = RHS[:,:,2] - params.cd*dX.*FFTW.fft(U.*v)+params.cd*dY.*FFTW.fft(U.*u);
RHS[:,:,2] = RHS[:,:,2] - params.r*Laplacian.*psi_hat(:,:,2);

# The following code computes the Jacobian J[psi,q], dealiased using a 3/2-rule.
# I.e., if you set N=256 then the code uses 3*N/2=384 Fourier modes (in each
# direction) to compute the Jacobian.
    # physical space, 3/2 grid; factor of (9/4) scales fft
    Psi_hat = zeros([1.5*params.N 1.5*params.N 2]);
    Psi_hat[1:params.N/2+1,1:params.N/2+1,:] = (9/4)*psi_hat(1:params.N/2+1,1:params.N/2+1,:);
    Psi_hat[1:params.N/2+1,params.N+2:1.5*params.N,:] = (9/4)*psi_hat(1:params.N/2+1,params.N/2+2:params.N,:);
    Psi_hat[params.N+2:1.5*params.N,1:params.N/2+1,:] = (9/4)*psi_hat(params.N/2+2:params.N,1:params.N/2+1,:);
    Psi_hat[params.N+2:1.5*params.N,params.N+2:1.5*params.N,:] = (9/4)*psi_hat(params.N/2+2:params.N,params.N/2+2:params.N,:);
    Q_hat = zeros([1.5*params.N 1.5*params.N 2]);
    Q_hat[1:params.N/2+1,1:params.N/2+1,:] = (9/4)*q_hat(1:params.N/2+1,1:params.N/2+1,:);
    Q_hat[1:params.N/2+1,params.N+2:1.5*params.N,:] = (9/4)*q_hat(1:params.N/2+1,params.N/2+2:params.N,:);
    Q_hat[params.N+2:1.5*params.N,1:params.N/2+1,:] = (9/4)*q_hat(params.N/2+2:params.N,1:params.N/2+1,:);
    Q_hat[params.N+2:1.5*params.N,params.N+2:1.5*params.N,:] = (9/4)*q_hat(params.N/2+2:params.N,params.N/2+2:params.N,:);
    # calculate u.gradq on 3/2 grid
    u = FFTW.irfft(-DY.*Psi_hat);
    v = FFTW.irfft( DX.*Psi_hat);
    qx= FFTW.irfft( DX.*Q_hat);
    qy= FFTW.irfft( DY.*Q_hat);
    jaco_real = u.*qx+v.*qy;
    # fft, 3/2 grid; factor of (4/9) scales fft
    Jaco_hat = (4/9)*FFTW.rfft(jaco_real);
    # reduce to normal grid
    jaco_hat = zeros([params.N params.N 2]);
    jaco_hat[1:params.N/2+1,1:params.N/2+1,:] = Jaco_hat[1:params.N/2+1,1:params.N/2+1,:];
    jaco_hat[1:params.N/2+1,params.N/2+2:params.N,:] = Jaco_hat[1:params.N/2+1,params.N+2:1.5*params.N,:];
    jaco_hat[params.N/2+2:params.N,1:params.N/2+1,:] = Jaco_hat[params.N+2:1.5*params.N,1:params.N/2+1,:];
    jaco_hat[params.N/2+2:params.N,params.N/2+2:params.N,:] = Jaco_hat[params.N+2:1.5*params.N,params.N+2:1.5*params.N,:];
    
# Put it all together
RHS = RHS - jaco_hat;
        
end