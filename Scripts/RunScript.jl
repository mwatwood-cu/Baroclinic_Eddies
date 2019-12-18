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

    println("Starting Ouput")
    io = open("/scratch/summit/mawa7160/Tests/starting500.txt","w+")
    writedlm(io,Init.qp[:,:,1])
    close(io)
    println("Starting time steps")
    flush(stdout)
    u_next = Init.q
    dt_next = 0.000001
    for i in 1:500
        u_next, dt_next = take_timestep_array(QG_RHS, Init.L, u_next, dt_next)
        if(i%50 ==0)
            println(string(i, "\t", dt_next, "\n")
            flush(stdout)
        end
        if(i%250 == 0)
            qp = FFTW.irfft(u_next, Init.parameters.N_points)
            io = open(string("/scratch/summit/mawa7160/Tests/state_",string(i),"_500.txt"), "w+")
            writedlm(io, qp[:,:,1])
            close(io)
        end
    end
    qp = FFTW.irfft(u_next, Init.parameters.N_points)
    io = open("/scratch/summit/mawa7160/Tests/final_state500.txt", "w+")
    writedlm(io, qp[:,:,1])
    close(io)
end

RunScript()
