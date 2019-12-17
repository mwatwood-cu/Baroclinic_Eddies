using FFTW, PyPlot, FileIO, DelimitedFiles
using Random, Distributions

include("RK_Functions.jl")
include("Initialize.jl")
include("Psi.jl");
include("Explicit_RHS.jl")
include("Stocastic.jl")

function RunScript()
    println("Strating Script")
    Setup_Psi()
    setup_QG_RHS();

    println("Starting Ouput")
    io = open("/scratch/summit/mawa7160/Tests/starting2.txt","w")
    writedlm(io,Init.qp[:,:,1])
    close(io)
    println("Starting time steps")
    u_next = Init.q
    dt_next = 0.000001
    for i in 1:4000
        u_next, dt_next = take_timestep_array(QG_RHS, Init.L, u_next, dt_next)
        if(i%50 ==0)
            println(string(i, "\t", dt_next, "\n"))
        end
        if(i%250 == 0)
            qp = FFTW.irfft(u_next, Init.parameters.N_points)
            io = open(string("/projects/mawa7160/Tests/state2_",string(i),".txt"), "w")
            writedlm(io, qp[:,:,1])
            close(io)
            #imshow(qp[:,:,1])
        end
    end
    qp = FFTW.irfft(u_next, Init.parameters.N_points)
    io = open("/scratch/summit/mawa7160/Tests/final_state2.txt", "w")
    writedlm(io, qp[:,:,1])
    close(io)
end

function test()
    println("Testing")
end
RunScript()
#test()
