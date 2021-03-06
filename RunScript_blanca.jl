using FFTW, PyPlot, FileIO, DelimitedFiles
using Random, Distributions

include("RK_Functions.jl")
include("Initialize.jl")
include("Psi.jl");
include("Explicit_RHS.jl")
include("Stocastic.jl")

function RunScript()
    println("Starting Script")
    Setup_Psi()
    setup_QG_RHS();

    io = open(string("/projects/mawa7160/Tests/state_4500_9000.txt"), "r")
    start  = readdlm(io, '\t', Float64, '\n')
    close(io)
    q = FFTW.rfft(start, 1:2);
    u_next = zeros(ComplexF64,641,1280,2)
    u_next[:,:,1] = q
    u_next[:,:,2] = q
    println("Starting time steps")
    flush(stdout)
    dt_next = 0.0008
    for i in 4500:9000
        u_next, dt_next = take_timestep_array(QG_RHS, Init.L, u_next, dt_next)
        if(i%50 ==0)
            println(string(i, "\t", dt_next, "\n"))
	    flush(stdout)
        end
        if(i%250 == 0)
            qp = FFTW.irfft(u_next, Init.parameters.N_points)
            io = open(string("/rc_scratch/mawa7160/Tests/state_",string(i),"_9000.txt"), "w+")
            writedlm(io, qp[:,:,1])
            close(io)
        end
    end
    qp = FFTW.irfft(u_next, Init.parameters.N_points)
    io = open("/rc_scratch/mawa7160/Tests/final_state9000.txt", "w+")
    writedlm(io, qp[:,:,1])
    close(io)
end

RunScript()
