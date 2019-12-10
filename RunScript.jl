using FFTW, PyPlot, FileIO, DelimitedFiles
using Random, Distributions

include("RK_Functions.jl")
include("Initialize.jl")
include("Psi.jl");
include("Explicit_RHS.jl")
include("Stocastic.jl")

function RunScript()
    Setup_Psi()
    setup_QG_RHS();


    io = open("\projects\mawa7160\Tests\starting.txt","w")
    writedlm(io,Init.qp[:,:,1])
    close(io)
    io = open("\projects\mawa7160\Tests\results.txt","w")
    print(io,"Iter\t Timestep\n")
    u_next = Init.q
    dt_next = 0.01
    for i in 1:50
        u_next, dt_next = take_timestep_array(QG_RHS, Init.L, u_next, dt_next)
        if(i%5 ==0)
            print(io, string(i, "\t", dt_next, "\n"))
            #qp = FFTW.irfft(u_next, Init.parameters.N_points)
            #imshow(qp[:,:,1])
        end
    end
    close(io)
    qp = FFTW.irfft(u_next, Init.parameters.N_points)
    io = open("\projects\mawa7160\Tests\final_state.txt", "w")
    writedlm(io, qp[:,:,1])
    close(io)
end

RunScript()
