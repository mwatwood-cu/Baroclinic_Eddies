module Init

    struct Phys_Params
        k_deform
        domain_wid
        k_beta
        r_eckman
        drag_coeff
        nu_buoy
        model_type
        N_points
    end
    using FFTW

    #   Set simulation parameters
    model = 1;
    # Number of grid points in each direction
    N=256
    # Compute diagnostics every countDiag steps
    countDiag = 5;
    # Is the time step adaptive
    isAdaptive = false;
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
    LX = 2*pi;
    # Nondimensional beta wavenumber; assumed to be 0 for 2-surface
    kb = 0;
    # Nondimensional Ekman friction coefficient
    r = 0.1;
    # Nondimensional quadratic drag coefficient
    c_d = 0.1;
    # Coefficient of hyperviscous PV/buoyancy diffusion
    nu = 1E-11;

    
    params = Phys_Params(kd, LX, kb, r, c_d, nu, model, N);

    # Set up hyperviscous PV dissipation
    # wavenumbers
    k = (2*pi/LX)*vcat(0:N/2, -N/2+1:-1);

    L = zeros(N, N);
    for j=1:N
        for i=1:N
            kr = sqrt(k[i]^2+k[j]^2);
            L[i,j] = -nu*kr^2;
        end
    end

    export k, L, params

    # Adaptive stepping
    tol= 1E-2;
    r0 = .8*tol;

    export tol, r0

    # Initialize
    t = 0;
    # qp is the "physical" potential vorticity, i.e. values on the grid
    qp = randn(N, N); 
    # q is the Fourier coefficients of potential vorticity
    q = fft(qp);
    # If model = 2 then qp is surface buoyancy on the grid

    # Initialize Tracer Diagnostics
    T = zeros(1, convert(Int64,Nt/countDiag));

    export t, qp, q, T

end