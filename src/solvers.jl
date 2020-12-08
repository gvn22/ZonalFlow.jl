## NL
function nl(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,
    νn::Float64=0.0;jw::Float64=0.1,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,kwargs...)
    A = acoeffs(ly,ny,Ξ,τ)
    B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
    Cp,Cm = ccoeffs(lx,ly,nx,ny)
    p = [nx,ny,A,B,Cp,Cm]
    tspan = (0.0,t_end)
    prob = ODEProblem(nl_eqs!,ic,tspan,p)
    @info "Solving NL equations on $(nx-1)x$(ny-1) grid"
    @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"
    solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false)
end

## GQL
function gql(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,νn::Float64=0.0;
    sfparams,jw::Float64=0.1,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,kwargs...)

    A = acoeffs(ly,ny,Ξ,τ)
    B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
    # Cp,Cm = ccoeffs(nx,ny)
    Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)

    tprev = 0.0
    trenew,mmin,mmax,amp = sfparams
    η = fill!(similar(ic),0.0 + 0.0im)
    η̂ = fill!(similar(ic),0.0 + 0.0im)
    fcoeffs!(nx,ny,mmin,mmax,amp,η̂)
    # η̂ = fcoeffs(nx,ny,mmin,mmax,amp)
    F = (tprev,trenew,η,η̂)

    p = [nx,ny,Λ,A,B,Cp,Cm,F]
    tspan = (0.0,t_end)
    prob = ODEProblem(gql_eqs!,ic,tspan,p)
    @info "Solving GQL equations on $(nx-1)x$(ny-1) grid with Λ = $Λ"
    @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"

    if trenew > 0.0

        integrator = init(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
        for i in integrator

            step!(integrator)

            tprev,trenew,η,η̂ = integrator.p[8]
            Δt = integrator.t - tprev
            onebytr = trenew > 0.0 ? 1.0/trenew : 0.0
            R = (1.0 - Δt*onebytr)/(1.0 + Δt*onebytr)
            η .= R*η .+ sqrt((1.0 - R*R)*onebytr)*η̂

            tprev = integrator.t
            # η̂ = fcoeffs(nx,ny,mmin,mmax,amp)
            fcoeffs!(nx,ny,mmin,mmax,amp,η̂)
            F = (tprev,trenew,η,η̂)
            integrator.p = [nx,ny,Λ,A,B,Cp,Cm,F]

        end
        return integrator.sol
    else
        solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
    end
end

## GCE2
function gce2(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,
    νn::Float64=0.0;jw::Float64=0.1,icnl::Bool=true,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,saveinfo::Bool=false,saveinfofreq::Int=50,poscheck::Bool=false,poscheckfreq::Int=20,kwargs...)
    A = acoeffs(ly,ny,Ξ,τ)
    B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
    Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)
    @info "Solving GCE2 equations on $(nx-1)x$(ny-1) grid with Λ = $Λ"
    @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"
    tspan = (0.0,t_end)
    u0 = icnl == true ? ic_cumulants(nx,ny,Λ,ic) : ic_cumulants(nx,ny,Λ,1e-6,ic)
    p = [nx,ny,Λ,A,B,Cp,Cm,fill!(similar(u0.x[1]),0),fill!(similar(u0.x[2]),0),fill!(similar(u0.x[2]),0)]
    prob = ODEProblem(gce2_eqs!,u0,tspan,p)
    if poscheck && Λ < nx - 1
        poschecktimes = [tt for tt=1.0:poscheckfreq:t_end]
        condition(u,t,integrator) = t ∈ poschecktimes && !ispositive(u.x[2],nx,ny,Λ)
        affect!(integrator) = positivity!(integrator.u.x[2],nx,ny,Λ)
        cbp = DiscreteCallback(condition,affect!,save_positions=(false,false))
        # return solve(prob,RK4(),callback=cbp,tstops=poschecktimes,dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
    end
    if saveinfo && Λ < nx - 1
        saved_values = SavedValues(Float64, Tuple{Float64,Array{Float64,1}})
        save_func(u,t,integrator) = rankis(integrator.u.x[2],nx,ny,Λ)
        saveat_array = [i for i=0.0:saveinfofreq:t_end]
        cbs = SavingCallback(save_func,saved_values,save_start=true,save_everystep=false,saveat=saveat_array)
    end
    if saveinfo && !poscheck
        @info "Saving twopoint correlation rank info"
        sol = solve(prob,RK4(),dt=dt,callback=cbs,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
        return sol,saved_values
    elseif saveinfo && poscheck
        @info "Running with positivity condition and saving twopoint correlation rank info"
        cb = CallbackSet(cbp,cbs)
        sol = solve(prob,RK4(),dt=dt,callback=cb,tstops=poschecktimes,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
        return sol,saved_values
    else
        return solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
    end
end

function ce2(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,
    νn::Float64=0.0;jw::Float64=0.1,icnl::Bool=true,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,saveinfo::Bool=false,saveinfofreq::Int=50,poscheck::Bool=false,poscheckfreq::Int=20,kwargs...)
    A = acoeffs(ly,ny,Ξ,τ)
    B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
    Cp,Cm = ccoeffs(lx,ly,nx,ny,0)
    @info "Solving CE2 equations on $(nx-1)x$(ny-1) grid"
    @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"
    tspan = (0.0,t_end)
    u0 = icnl == true ? ic_cumulants(nx,ny,ic) : ic_cumulants(nx,ny,1e-6,ic)
    p = [nx,ny,A,B,Cp,Cm,fill!(similar(u0.x[1]),0),fill!(similar(u0.x[2]),0),fill!(similar(u0.x[2]),0)]
    prob = ODEProblem(ce2_eqs!,u0,tspan,p)
    if saveinfo
        saved_values = SavedValues(Float64, Tuple{Int64,Array{Float64,1},Array{Int64,1},Array{Float64,2}})
        save_func(u,t,integrator) = rankis(integrator.u.x[2],nx,ny)
        saveat_array = [i for i=0.0:saveinfofreq:t_end]
        cbs = SavingCallback(save_func,saved_values,save_start=true,save_everystep=false,saveat=saveat_array)
        sol = solve(prob,RK4(),dt=dt,callback=cbs,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
        return sol,saved_values
    end
    solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
end
