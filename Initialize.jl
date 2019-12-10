module Init

    mutable struct Phys_Params
        k_deform::Float64
        domain_wid::Float64
        k_beta::Float64
        r_ekman::Float64
        drag_coeff::Float64
        nu_buoy::Float64
        model_type::Integer
        N_points::Integer
    end

    mutable struct Psi_Parameters
        InvBaroTropic 
        InvBaroClinic 
        InvMat
    end

    mutable struct Laplacian_Parameters
        dX
        dY
        DX
        DY
        Laplacian
    end

    mutable struct Time_Step_Parameters
        r0
    end

    psi_param = Psi_Parameters(0, 0, 0)
    laplace_param = Laplacian_Parameters(0, 0, 0, 0, 0)

    export psi_param, lap_param
    
    using FFTW

    #   Set simulation parameters
    model = 1;
    # Number of grid points in each direction
    N=128
    # Compute diagnostics every countDiag steps
    countDiag = 5;
    # Is the time step adaptive
    isAdaptive = true;
    # Starting time step; will change if adaptive
    dt = .01;
    # Number of Time steps
    Nt = 1e3;
    # Blowup threshold on q
    qlim = 1e3;

    export countDiag, isAdaptive, dt, Nt, qlim

    #  Set physical parameters
    # Nondimensional deformation wavenumber. 2-surface assumes kd=1.
    kd = 1; 
    # Nondimensional domain width
    LX = 2*pi
    Ly = 2*pi
    # Nondimensional beta wavenumber; assumed to be 0 for 2-surface
    kb = 0;
    # Nondimensional Ekman friction coefficient
    r_ekman = 0.5;
    # Nondimensional quadratic drag coefficient
    c_d = 0.1;
    # Coefficient of hyperviscous PV/buoyancy diffusion
    nu = 1E-11;

    
    parameters = Phys_Params(kd, LX, kb, r_ekman, c_d, nu, model, N);

    # Set up hyperviscous PV dissipation
    # wavenumbers
    k = (2*pi/LX)*vcat(0:N/2, -N/2+1:-1);

    L = zeros(Int(div(N,2)+1), N, 2);
    for i=1:Int(N/2+1)
        for j=1:N
            kr = sqrt(k[i]^2+k[j]^2);
            L[i,j,1] = -nu*kr^2;
            L[i,j,2] = -nu*kr^2;
        end
    end

    export k, L, parameters

    # Adaptive stepping
    tol= 1E-2;
    r0 = .8*tol;
    time_parameters = Time_Step_Parameters(r0)
    export tol, time_parameters

    # Initialize
    t = 0;
    # qp is the "physical" potential vorticity, i.e. values on the grid
    qp = randn(N, N, 2); 
    # q is the Fourier coefficients of potential vorticity
    q = FFTW.rfft(qp, 1:2);
    # If model = 2 then qp is surface buoyancy on the grid

    # Initialize Tracer Diagnostics
    T = zeros(1, convert(Int64,Nt/countDiag));

    export t, qp, q, T

    #Stocastic bounds
    kMin = 19*pi
    kMax = 21*pi

    #Creating Stochastic Forcing Spectrum
    kx = zeros(Int(N/2+1), N)
    ky = zeros(Int(N/2+1), N)
    k2 = zeros(Int(N/2+1), N)

    #Initialize wavenumbers
    for i=1:Int(N/2+1)
        kx[i,:] = ones(N).*(2*pi/LX) * (i-1)
    end

    for j=1:Int(N/2+1)
        ky[:,j] = ones(Int(N/2+1)).*(2*pi/Ly) * (j-1)
    end
    for j=Int(N/2+2):N
        ky[:,j] = ones(Int(N/2+1)).*(2*pi/Ly) * (j-N-1)
    end

    k2 = kx.^2+ky.^2
    stoch_spec = zeros(Int(N/2+1), N)
    for j=1:N
       for i=1:Int(N/2+1)
           if ( (kMin^2<=k2[i,j]) && (k2[i,j]<=kMax^2) ) 
               stoch_spec[i,j] = 1/sqrt.(k2[i,j])
           end
       end
   end
   stoch_spec[1,1] = 0

   dy = Ly/N

   y = zeros(N,N) 
   
   for j=1:N
        y[:,j] = ones(N).*(-(Ly+dy)/2 + dy*j) 
   end

   export stoch_spec, k2, y, Ly
end