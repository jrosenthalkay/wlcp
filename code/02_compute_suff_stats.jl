"""

Code for the computation of sufficient statistics for model counterfactuals

"The Winners and Losers of Climate Policies: A Sufficient Statistics Approach"
by Thomas Bourany and Jordan Rosenthal-Kay

Date created: February 2025
Last updated: September 2025

ran on the following environment: 
    Julia Version 1.11.1
    Commit 8f5b7ca12ad (2024-10-16 10:53 UTC)
    Build Info:
    Official https://julialang.org/ release
    Platform Info:
    OS: macOS (arm64-apple-darwin22.4.0)
    CPU: 10 × Apple M2 Pro
    WORD_SIZE: 64
    LLVM: libLLVM-16.0.6 (ORCJIT, apple-m2)
    Threads: 1 default, 0 interactive, 1 GC (on 6 virtual cores)
    Environment:
    JULIA_EDITOR = code
    JULIA_VSCODE_REPL = 1
    
"""

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --- set environment ------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# >>>>>>>> SET YOUR ROOT FOLDER HERE <<<<<<<<
cd("C:/Users/path/to/root/wlcp") # <-- Change this to the repo root on your machine

# packages
using Pkg, Base.Threads # Pkg must be installed, so as to install other dependencies

function require(pkg::String)
    # check if package exists. if not, install
    if !haskey(Pkg.project().dependencies,pkg)
        println("Installing $pkg...")
        Pkg.add(pkg)
    end
    @eval using $(Symbol(pkg))
end

Pkg_list = ["Parameters",
    "LinearAlgebra",
    "Plots",
    "BasisMatrices",
    "SparseArrays",
    "QuantEcon",
    "Arpack",
    "Roots",
    "ForwardDiff",
    "LaTeXStrings",
    "CSV",
    "DataFrames",
    "Optim"
]

# load / install packages
for pkg in Pkg_list 
    require(pkg)
end

# pin NL solve for julia 1.11 or greater issues 
# NLsolve / forwarddiff compatability 
Pkg.pin(Pkg.PackageSpec(name="NLsolve", version="4.5.1"))
require("NLsolve")

# create a blank .txt file that records statistics that do not appear in figures / tables.
outfile = "output/"
write(outfile*"SCC_ests.txt", "")

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --- load and clean data --------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# read in sufficient statistics csv
data3_raw = CSV.read("data/model_input/baseline_csv_suff_stat.csv", DataFrame)
data3 = zeros(193,19) ; 

# rescale gdp and population numbers
data3_raw.population = data3_raw.population./1e6
data3_raw.gdp_per_capita = data3_raw.gdp_per_capita/1000

# fill data matrix
data3[:,1] = data3_raw.population 
data3[:,2] = data3_raw.gdp_per_capita
data3[:,3] = data3_raw.primary_energy_consumption ./data3_raw.population # from TWh/million ppl = MWh/capita
data3[:,4] = data3_raw.fossil_share
data3[:,5] = data3_raw.coal_share
data3[:,6] = data3_raw.renewable_share
data3[:,7] = data3_raw.e_x ./data3_raw.population
data3[:,8] = data3_raw.consumption_share; 
data3[:,9] = data3_raw.fossil_rent_share
data3[:,10] = data3_raw.coal_rent_share
data3[:,11] = 1 ./ data3_raw.inv_nu_f_trunc
data3[:,12] = 1 ./ data3_raw.inv_nu_c_trunc
data3[:,13] = 0.0181.*data3_raw.e_x  ./ (data3_raw.gdp_per_capita .* data3_raw.population) # vex, rescale carbon intensity
data3[:,14] = 0.0181.*data3_raw.primary_energy_consumption .*data3_raw.fossil_share ./ (data3_raw.gdp_per_capita .* data3_raw.population) # vef
data3[:,15] = data3_raw.nrg_export_share
data3[:,16] = data3_raw.nrg_import_share
data3[:,17] = data3_raw.temp
data3[:,18] = data3_raw.co2_em_mt
data3[:,19] = data3_raw.carbon_intensity

# flag important iso3 codes for plotting and use for counterfactuals
iso3 = ["USA", "CAN", "CHN", "DEU", "ESP", "FRA", "GBR", "ITA", "IND", "PAK", "NGA", "ZAF", "EGY", "IRN", "SAU", "TUR", "RUS", "AUS", "JPN", "KOR", "IDN", "THA", "ARG", "BRA", "MEX"]

iso3_l = data3_raw.iso3 ;
iso3_ls = iso3[1:25] ;
ind_ls = zeros(Int,25); 

for i = 1:25
    ind_ls[i] = Int( sum((1:193) .*(iso3_l.==iso3_ls[i])) )
end

iso3_25 = iso3[1:25]
iso3_25w = [iso3_25 ; "world"]

# Identify blocs of countries
# Current EU and ASEAN member states as of March 2025 (ISO3 codes)
eu_members = ["AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE"]
asean_members = ["BRN", "KHM", "IDN", "LAO", "MYS", "MMR", "PHL", "SGP", "THA", "VNM"]

# Create a vector of 0s and 1s (1 if in EU, 0 if not; ditto for ASEAN)
eu_membership = [iso_3 in eu_members ? 1 : 0 for iso_3 in iso3_l]
asean_membership = [iso_3 in asean_members ? 1 : 0 for iso_3 in iso3_l]

is_eu = Bool.(eu_membership)
is_asean = Bool.(asean_membership)

# Clean trade matrix 
# unique iso3 list for trade matrix data 
unique_data1 = unique(data3_raw.iso3)

# ingest trade matrix data
data3_trade_raw = CSV.read("data/model_input/trade_shares.csv", DataFrame)

# origins and destinations
unique_data2 = unique(data3_trade_raw.iso3_o) ; 
unique_data2_ = unique(data3_trade_raw.iso3_d) ; 

# sanity check
setdiff(unique_data1, unique_data2) 
setdiff(unique_data1, unique_data2_) 
setdiff(unique_data2,unique_data1) 

# filter trade data
data3_trade_raw_filter_o = filter(row -> row.iso3_o in data3_raw.iso3, data3_trade_raw)
data3_trade_raw_filter_od = filter(row -> row.iso3_d in data3_raw.iso3, data3_trade_raw_filter_o)

# ffrom to 38809 to 33352
# 193^2 = 37249; 
dftrade = data3_trade_raw_filter_od;

# Get the unique sorted list of iso3 codes
iso3_list = sort(union(dftrade.iso3_o, dftrade.iso3_d))

@assert all(iso3_list .== iso3_l) # ensure matrix rows line up.

# Map iso3 codes to indices
iso3_indices = Dict(code => idx for (idx, code) in enumerate(iso3_list))

# Initialize the matrix with zeros
datatrademat = zeros(Float64, 193, 193)

# Fill the matrix with s_od values
for row in eachrow(dftrade)
    i = iso3_indices[row.iso3_o]
    j = iso3_indices[row.iso3_d]
    datatrademat[j,i] = row.trade_shares
end

# check 
sum(datatrademat,dims=2) # not quite numerically 1 -- not quite; will handle after checking for 0s on diagonal

# check for zeros on the diagonal
for i=1:193
    if datatrademat[i,i] <0.2
        println("iso3: ", iso3_list[i], " dom trade: ",round(datatrademat[i,i], digits=5), "        share wld pop (pc) ",round(100*data3[i,1]./sum(data3[:,1]), digits=5) , "          share wld gdp ", round(100*data3[i,1].*data3[i,2]./sum(data3[:,1].*data3[:,2]), digits=5))
    end 
end 
# there are none; good.

# flag domestic trade = 0 (sanity check; this does not occur)
datatradematnew = copy(datatrademat)

datatradematnew = datatradematnew./sum(datatradematnew,dims=2) # ensure row sums are one

# now, plot the trade matrix 

# sanity check: view full matrix: many zeros
heatmap(1:193,1:193,datatradematnew, xrotation = 30, xtickfont = font(6, "Arial"),ytickfont = font(6, "Arial"), xticks = (1:193,iso3_list), yticks = (1:193,iso3_list))

# small trade matrix of major players
datatrademat_sm = datatradematnew[ind_ls,ind_ls]

ps=heatmap(1:25,1:25, 
    reverse(datatrademat_sm.^0.32,dims=1), 
    color = cgrad(:Blues_5, rev = false),
    xticks=(1:25, iso3_ls), xrotation = 40,
    yticks=(1:25, reverse(iso3_ls)), 
    legend=:none, xtickfontsize =7, ytickfontsize = 7,
    title = title = "Trade share matrix, " * string(L"\mathbf{S}"),
    title_location = :left)

l = @layout [a{0.95w} b{0.05w}]

p = plot(ps, heatmap((0:0.01:1).*ones(101,1), 
    legend=:none, xticks=:none, 
    yticks=(1:10:101, vcat(string.([0,0.3,1.8]),string.(Int.(round.(100*(0.3:0.1:1).^3.125))))),
    color = cgrad(:Blues_5, rev = false)), 
    layout=l) 

savefig("output/figures/data_tradeshare_matrix_v5_25.png")

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --- define model objects -------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# structure for model parameters
@with_kw mutable struct Params_IAMHAT_suffstats_o

    # handle data
    data = data
    datatrade = datatrade 
    Ntot = size(data,1)

    pp = data[:,1]
    PP = sum(pp)
    p̂p̂ = pp ./ PP ;
    tt =  data[:,17]
    yy =  data[:,2]
    yy_pp = yy .* pp
    ŷŷ = yy_pp ./ sum(yy_pp) ;

    #1 TWh = 85,984 Toe
    ee      =  data[:,3]
    sh_ef   =  data[:,4]
    sh_ec   =  data[:,5]
    sh_er   =  data[:,6]
    ϵϵ      =  data[:,18]

    exex    = data[:,7].* pp ./ sum( pp.* data3[:,7]) 
    ex      = data[:,7].* pp 
    ef      = pp.*ee.*sh_ef


    ΠΠex = data[:,9]
    ΠΠec = data[:,10]

    #Energy data (MTOE)
    #EE = 13961.33 ;  # total energy in MTOE 
    #Share emissions/energy
    EE = sum(ee .* pp) 
    # =132708 TWh = 11410 MTOE

    #Emission data tons of carbon (MT)
    #EmEm::Float64 = 9945.75 ; 
    #EmCO2::Float64 = 36467.76
    #Em_ghg_others::Float64 = 4360
    EmEm = sum(data3[:,18])
    # = 30454 TCO2 
    S0 = 484000 * 44/12

    sh_Ef = 0.56104
    sh_Ec = 0.26935
    sh_Er = 0.169599

    # carbon intensity 
    # convert kg of C02 emitted per dollar of GDP to metric ton / $
    em_intensity = data3[:,19] ./ 1000

    # calibrated parameters 
    θ = 5.0
    sh_E = 0.11; 

    ηef = ΠΠex ; 
    ηec = ΠΠec ;  ; 
    ηer = 0.005 .*ones(Ntot) ; 
    ηy = 1 .- ηef .-  ηec .- ηer ;
    ηc = 0.95 .* ones(Ntot); 
    
    # carbon intensities 
    #2.76 = (32.192 .* 3.068 .+ 23.913.* 2.349 ) ./ (32.192 + 23.913)
    ### ξf = 2.761  / (44/12)    
    ### ξc = 3.961 / (44/12) 
    ξf = 0.238   
    ξc = 0.341 

    # energy supply parameters
    σ = 0.6
    σe = 2.0 
    νf = data[:,11] 
    νc = data[:,12] 
    νr = 2.7*ones(Ntot) 
    ν̄f =  1 ./(sum(exex ./ νf)) ;

    #other sufficient stats: 
    λef = pp .* ee .* sh_ef ./sum(pp .* ee .* sh_ef)
    λec = pp .* ee .* sh_ec ./sum(pp .* ee .* sh_ec)
   
    s_fϵ = EE.*sh_Ef.*ξf ./ (EE.*sh_Ef.*ξf .+ EE.*sh_Ec.*ξc)
    s_cϵ = EE.*sh_Ec.*ξc ./ (EE.*sh_Ef.*ξf .+ EE.*sh_Ec.*ξc)

    T_horizon = 92
    s_ES = T_horizon .* EmEm ./ (S0 .+ T_horizon.*EmEm) 
    
    tstar = 14.016
    t̄ = 0.5 .* tstar .+ 0.5 * tt
    # κ+ = 0.00311, and κ− = 0.00456. Krusell Smith 2022
    # κ+ = κ− = 0.002388    Kotlikoff et al climate policy and welfare
    #γ2₊::Float64 = 1.25*0.0028388  #1.25*0.0028388 
    #γ2₋::Float64 = 0.3*0.0028388    #0.002337157

    γdam = 0.012
    γ̄ = 2 .* γdam .* (tt - t̄) .* tt .* s_ES

    # Price TOE 
    # 1 barrel = 0.1364 TOE
    # 1 barrel of oil: 70$
    # 1 TOE of oil: 513$

    # 1000 cubic feet = 0.024 TOE
    # 1000 cubic feet: 4-7$ (=5)
    # 1 TOE of gas: 208$

    # TOE of oil and gas: 370 
    # (200*24+500*32)/(32+24)
    # = 371.

    #ex= EE.* sh_Ef.*exex
    #ef = pp.*ee.*sh_ef
    #yp = (1 .- ΠΠex).* yy

    q̄f = 0.0319
    vy = data[:,8] # yp/v

    vex = data[:,15]
    vef = data[:,16]
    vne = data[:,15] .- data[:,16]
    v = (ηy .* yy + vne .* yy).*pp

    Smat = datatrade ;
    Tmat_v = (Smat' .* ((v)' ./ (v))); 
    Tmat_y_n = (Smat' .* ((yy_pp)' ./ (yy_pp))); 

    Tmat_y = Tmat_y_n ./ sum(Tmat_y_n,dims=2); 

    # utility paramaters
    ηcrra = 1.5 
    ω_u = pp .* yy.^(1 - ηcrra) ./ sum(pp .* yy.^(1 - ηcrra))
    ω_n = ŷŷ

end

# output container (?)
@with_kw mutable struct Model_suffstats_o

    par::Params_IAMHAT_suffstats_o

    Ntot = par.Ntot 
    dϵ = 0.0;
    dϵ_eq = 0.0;

    dtϵ = zeros(Ntot)
    dsϵ = zeros(Ntot)
    dtb = zeros(Ntot,Ntot)
    𝕁 = zeros(Ntot)

    dD = zeros(Ntot)
    dp = zeros(Ntot)
    dP = zeros(Ntot)
    dP_r = zeros(Ntot)
    dy = zeros(Ntot)
    dqf = 0;
    dqc = zeros(Ntot)
    dqr = zeros(Ntot)
    def = zeros(Ntot)
    dec = zeros(Ntot)
    der = zeros(Ntot)
    dπf = zeros(Ntot)
    dπc = zeros(Ntot)
    dπr = zeros(Ntot)
    dex = zeros(Ntot)
    dEf = 0
    dEc = 0 

    dW_Dy   = zeros(Ntot)
    dW_Du   = zeros(Ntot)
    dW_p    = zeros(Ntot)
    dW_P    = zeros(Ntot)
    dW_ef   = zeros(Ntot)
    dW_ec   = zeros(Ntot)
    dW_er   = zeros(Ntot)
    dW_πf   = zeros(Ntot)
    dW_πc   = zeros(Ntot)
    dW_πr   = zeros(Ntot)
    dWtot   = zeros(Ntot)
    dWdec   = zeros(Ntot,5)
    dW      = zeros(Ntot,11)

    dW_wld_u = 0 
    dW_Dy_wld_u = 0 
    dW_p_wld_u = 0 
    dW_e_wld_u = 0 
    dW_π_wld_u = 0  

    dW_wld_n   = 0 
    dW_Dy_wld_n= 0 
    dW_p_wld_n= 0 
    dW_e_wld_n= 0 
    dW_π_wld_n = 0

    A = zeros(Ntot,Ntot)

    dx =hcat(dtϵ,dsϵ,𝕁,dD ,dp ,dP ,dy ,dqc,dqr,def,dec,der,dπf,dπc,dπr,dex,dW_Dy,dW_Du,dW_p ,dW_P ,dW_ef,dW_ec,dW_er,dW_πf,dW_πc,dW_πr,dWtot)
    dX = hcat(dϵ,dϵ_eq,dqf,dEf,dEc, dW_wld_u,dW_Dy_wld_u,dW_p_wld_u,dW_e_wld_u,dW_π_wld_u ,  dW_wld_n, dW_Dy_wld_n,dW_p_wld_n,dW_e_wld_n,dW_π_wld_n )
    name_varx = ["dtϵ","dsϵ","J","dD ","dp ","dP ","dy ","dqc","dqr","def","dec","der","dπf","dπc","dπr","dex","dW_Dy","dW_Du","dW_p ","dW_P ","dW_ef","dW_ec","dW_er","dW_πf","dW_πc","dW_πr","dWtot"]
    name_varX = ["dϵ","dϵ_eq","dqf","dEf","dEc", "dW_wld_u"," dW_Dy_wld_u", " dW_p_wld_u", " dW_e_wld_u", " dW_π_wld_u ", "   dW_wld_n", "  dW_Dy_wld_n", " dW_p_wld_n", " dW_e_wld_n", " dW_π_wld_n"]
end 

# set up model parameters
parss = Params_IAMHAT_suffstats_o(data = data3,datatrade=datatradematnew ) ; 

# print income share matrix
ps=heatmap(1:25,1:25, 
    reverse(parss.Tmat_y[ind_ls,ind_ls].^0.32,dims=1), 
    color = cgrad(:Blues_5, rev = false),xticks=(1:25, iso3_ls),
    xrotation = 40,yticks=(1:25, reverse(iso3_ls)), 
    legend=:none, xtickfontsize =7, ytickfontsize = 7,
    title = title = "Income flow matrix, " * string(L"\mathbf{T}"),
title_location = :left)
l = @layout [a{0.95w} b{0.05w}]
p = plot(ps, heatmap((0:0.01:1).*ones(101,1), 
    legend=:none, xticks=:none, 
    yticks=(1:10:101, vcat(string.([0,0.3,1.8]),string.(Int.(round.(100*(0.3:0.1:1).^3.125))))),
    color = cgrad(:Blues_5, rev = false)), 
    layout=l) # Plot them set y values of color bar accordingly

savefig("output/figures/data_incomeshare_matrix_v5_25.png")

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --- set up model functions -----------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# compute energy block effects
function energy_qe_cr(par::Params_IAMHAT_suffstats_o)
    
    @unpack Ntot = par
    @unpack pp,yy,tt,yy_pp,ϵϵ,exex, EE, EmEm, sh_Ef,sh_Ec, sh_Er, sh_E, ΠΠex, sh_ef, sh_ec, sh_er = par
    @unpack ξf,ξc,σ,σe,νf,νc,νr,ν̄f,λef,λec,s_fϵ,s_cϵ, γ̄,ex,ef = par
    @unpack Smat,Tmat_y, θ = par

    e_c = zeros(Ntot,7)
    e_r = zeros(Ntot,7)
    q_c = zeros(Ntot,5)
    q_r = zeros(Ntot,5)

    for i = 1:Ntot
        B = [ -((σ/(1 - sh_E))*sh_ec[i]+(1-sh_ec[i])*σe)   (σe - (σ/(1 - sh_E)))*sh_er[i] ; (σe - (σ/(1 - sh_E)))*sh_ec[i]    -((σ/(1 - sh_E))*sh_er[i]+(1-sh_er[i])*σe)   ]
        νmat = [νc[i] 0 ; 0  νr[i]]

        Bν_inv =  inv( I(2) .- B*νmat )

        # 1. Climate Damage
        e_c[i,1] = ((σ/(1 - sh_E)) * Bν_inv * ones(2) )[1]
        e_r[i,1] = ((σ/(1 - sh_E)) * Bν_inv * ones(2) )[2]

        # 2. Good prices
        e_c[i,2] = (Bν_inv * (B .+ σ/(1 - sh_E)) * ones(2) )[1]
        e_r[i,2] = (Bν_inv * (B .+ σ/(1 - sh_E)) * ones(2) )[2]

        # 3. Oil-gas price 
        e_c[i,3] = ((σe - (σ/(1 - sh_E)))*sh_ef[i] *  Bν_inv * ones(2) )[1]
        e_r[i,3] = ((σe - (σ/(1 - sh_E)))*sh_ef[i] *  Bν_inv * ones(2) )[2]

        # 4. Carbon tax 
        e_c[i,4] = ((σe - (σ/(1 - sh_E)))*sh_ef[i] *  Bν_inv * ones(2) .* ξf )[1] .+ ((σ/(1 - sh_E)) * Bν_inv * B * [ξc ; 0] )[1]
        e_r[i,4] = ((σe - (σ/(1 - sh_E)))*sh_ef[i] *  Bν_inv * ones(2) .* ξf )[2] .+ ((σ/(1 - sh_E)) * Bν_inv * B * [ξc ; 0] )[2]

        e_c[i,6] = ((σe - (σ/(1 - sh_E)))*sh_ef[i] *  Bν_inv * ones(2) .* ξf )[1]
        e_c[i,7] =  ((σ/(1 - sh_E)) * Bν_inv * B * [ξc ; 0] )[1]

        e_r[i,6] = ((σe - (σ/(1 - sh_E)))*sh_ef[i] *  Bν_inv * ones(2) .* ξf )[2] 
        e_r[i,7] =  ((σ/(1 - sh_E)) * Bν_inv * B * [ξc ; 0] )[2]


        # 5. Renewable subsidy 
        e_c[i,5] =  ((σ/(1 - sh_E)) .* Bν_inv * B * [0 ; -1] )[1]
        e_r[i,5] =  ((σ/(1 - sh_E)) .* Bν_inv * B * [0 ; -1] )[2]

        # energy prices
        q_c[i,:] = νc[i] .* e_c[i,1:5] 
        q_c[i,2] = q_c[i,2] .+ 1

        q_r[i,:] = νr[i] .* e_r[i,1:5]
        q_r[i,2] = q_r[i,2] .+ 1

    end 
    return e_c,e_r,q_c,q_r
end

# compute energy market effects, section 3.6
function energy_ef(par::Params_IAMHAT_suffstats_o, qer,qec)
    
    @unpack Ntot = par
    @unpack pp,yy,tt,yy_pp,ϵϵ,exex, EE, EmEm, sh_Ef,sh_Ec, sh_Er, sh_E, ΠΠex, sh_ef, sh_ec, sh_er = par
    @unpack ξf,ξc,σ,σe,νf,νc,νr,ν̄f,λef,λec,s_fϵ,s_cϵ, γ̄,ex,ef = par
    @unpack Smat,Tmat_y, θ = par

    e_f = zeros(Ntot,5)

    e_f[:,1] = (σ/(1 - sh_E)) .*ones(Ntot)
    e_f[:,2] = (σ/(1 - sh_E)) .*ones(Ntot)
    e_f[:,3] = -(σe - (σ/(1 - sh_E)))*sh_ef
    e_f[:,4] = -(σe - (σ/(1 - sh_E)))*sh_ef .* ξf .+ (σe - (σ/(1 - sh_E)))*sh_ec .* ξc
    e_f[:,5] = -(σe - (σ/(1 - sh_E)))*sh_er

    e_f_d = e_f 
    
    e_f_i = e_f_d .+ (σe - (σ/(1 - sh_E))).* ( sh_ec .* qec .+ sh_er .* qer) 

    return e_f_i, e_f_d 
end 

# compute output changes, section 3.6
function  output_y(par::Params_IAMHAT_suffstats_o, qec,qer)
    
    @unpack Ntot = par
    @unpack pp,yy,tt,yy_pp,ϵϵ,exex, EE, EmEm, sh_Ef,sh_Ec, sh_Er, sh_E, ΠΠex, sh_ef, sh_ec, sh_er = par
    @unpack ξf,ξc,σ,σe,νf,νc,νr,ν̄f,λef,λec,s_fϵ,s_cϵ, γ̄,ex,ef = par
    @unpack Smat,Tmat_y, θ = par

    y_d = zeros(Ntot,5)
    y_i = zeros(Ntot,5)

    αy = (σ .* sh_E/(1 - sh_E)); 
    y_d[:,1] = (1 .+ αy) .*ones(Ntot)
    y_d[:,2] = αy .*ones(Ntot)
    y_d[:,3] = -αy .*sh_ef
    y_d[:,4] = -αy .*sh_ef.*ξf .- αy*sh_ec .* ξc
    y_d[:,5] = αy .*sh_er

    y_i = y_d .-  αy .*sh_ec .* qec .- αy .*sh_er .* qer

    return y_i, y_d
end 

# compute total energy market effects (e.g., equation 24)
function totalenergy_qef(par::Params_IAMHAT_suffstats_o, ∂ef, ∂ec, dϵ̄, dp, dtϵ, dsϵ)
    
    @unpack Ntot = par
    @unpack pp,yy,tt,yy_pp,ϵϵ,exex, EE, EmEm, sh_Ef,sh_Ec, sh_Er, sh_E, ΠΠex, sh_ef, sh_ec, sh_er = par
    @unpack ξf,ξc,σ,σe,νf,νc,νr,ν̄f,λef,λec,s_fϵ,s_cϵ, γ̄,ex,ef = par
    @unpack Smat,Tmat_y, θ = par

    dD =  .- γ̄ .* dϵ̄

    λqef = sum(λef .* ∂ef[:,3])

    def = ∂ef[:,1].*dD .+  ∂ef[:,2].*dp .+  ∂ef[:,4].*dtϵ .+ ∂ef[:,5].*dsϵ
    dec = ∂ec[:,1].*dD .+  ∂ec[:,2].*dp .+  ∂ec[:,4].*dtϵ .+ ∂ec[:,5].*dsϵ
    
    dqf = (1 ./ (1 .- ν̄f.*λqef)) .* ( ν̄f .* sum(λef .* def) .+ sum(exex.*(ν̄f./νf).*dp)  ) 

    dEf = sum(λef .* def)
    dEc = sum(λec .* dec)

    dϵ̄2 = s_fϵ .*dEf .+  s_cϵ .*dEc 

    #dlnexi =[(1/ν)dlnqf −dlnpi]

    return dqf, dEf,dEc, dϵ̄2
end 

# compute trade block, equation (25)
function tradeblock_p(par::Params_IAMHAT_suffstats_o, ∂y, ∂ef, dϵ̄, dqf, dtϵ, dsϵ, dtb, 𝕁)
    
    @unpack Ntot = par
    @unpack pp,yy,tt,yy_pp,ϵϵ,exex, EE, EmEm, sh_Ef,sh_Ec, sh_Er, sh_E,ΠΠex, sh_ef, sh_ec, sh_er = par
    @unpack ξf,ξc,σ,σe,νf,νc,νr,ν̄f,λef,λec,s_fϵ,s_cϵ, γ̄,ex,ef = par
    @unpack Smat,Tmat_y, θ = par
    @unpack vy, vef, vex, vne = par
    @unpack em_intensity = par

    ∂ex = zeros(Ntot,5)
    ∂ex[:,3] = (1 ./ νf)  
    ∂ex[:,2] = .- ones(Ntot)
    dex = ∂ex[:,3].*dqf


    # trade tariff of non-members
    ind𝕁 = findall(==(1),𝕁) ;
    j_in = zeros(Ntot, Ntot); 
    j_in[ind𝕁,:] .= 1.0 ; 
     
    #em_intensity = (ϵϵ ./(pp.* yy))./1000    # carbon intensity (in ton CO2 /$) for one unit of output
    𝐭bc_ij = ones(Ntot, Ntot) .* (j_in.*(1 .-j_in')).*em_intensity'.*dtb  ;

    Jtb = 𝐭bc_ij

    dD =  .- γ̄ .* dϵ̄

    def = ∂ef[:,1].*dD .+  ∂ef[:,3].*dqf .+  ∂ef[:,4].*dtϵ .+ ∂ef[:,5].*dsϵ
    dy  = ∂y[:,1].*dD .+ ∂y[:,3].*dqf .+  ∂y[:,4].*dtϵ .+ ∂y[:,5].*dsϵ

    C = Tmat_y.* vy' .- I(Ntot)

    A = θ.*I(Ntot) .- Tmat_y.* vy' .- (θ-1)*(Tmat_y*Smat) .-  Tmat_y .* (vex .* ∂ex[:,2] .-vef.*∂ef[:,2])' .- C .*(∂y[:,2]')

    RHS_output = C * dy 
    RHS_energy = Tmat_y * (vex .* dex .- vef .* def .+ vne .* dqf )  
    RHS_tariffs = θ* ( Tmat_y*(Smat .*𝐭bc_ij)*ones(Ntot) - diag( Tmat_y*((ones(Ntot,Ntot) .+ Smat').* 𝐭bc_ij')) )

    # solve
    Ainv =inv(A) 

    dp_output = Ainv* (RHS_output )
    dp_energy = Ainv * (RHS_energy)
    dp_tariffs = Ainv * (RHS_tariffs)
    dp = dp_output .+ dp_energy .+ dp_tariffs

    return dp, dp_output, dp_energy, dp_tariffs, 𝐭bc_ij

end 

# GE function wrapper
function GE_fct(par::Params_IAMHAT_suffstats_o, dp0, dϵ̄, dtϵ, dsϵ, dtb, 𝕁, resid)

    # resid: GE in dϵ̄
    
    ∂ec,∂er,∂qc,∂qr = energy_qe_cr(par)
    ∂ef,∂ef_d = energy_ef(par,∂qc, ∂qr)
    ∂y, ∂y_d = output_y(par, ∂qc, ∂qr)

    dqf, dEf,dEc, dϵ̄1 = totalenergy_qef(par, ∂ef, ∂ec, dϵ̄, dp0, dtϵ, dsϵ)

    dp, dp_y,dp_e,dp_t, 𝐭bc_ij = tradeblock_p(par, ∂y, ∂ef, dϵ̄, dqf, dtϵ, dsϵ, dtb, 𝕁)

    resid_p = dp .- dp0
    resid_ϵ̄ = dϵ̄1 .- dϵ̄

    if resid
        return vcat(resid_p,resid_ϵ̄)
    else
        return resid_p
    end 
end 


# compute GE solution 
function GE_solution(par::Params_IAMHAT_suffstats_o, dp, dϵ̄, dtϵ, dsϵ, dtb, 𝕁)
        
    @unpack Ntot = par
    @unpack pp,yy,tt,yy_pp,ϵϵ,exex, EE, EmEm, sh_Ef,sh_Ec, sh_Er, sh_E, ΠΠex,ΠΠec, sh_ef, sh_ec, sh_er = par
    @unpack ξf,ξc,σ,σe,νf,νc,νr,ν̄f,λef,λec,s_fϵ,s_cϵ, γ̄,ex,ef = par
    @unpack Smat,Tmat_y, θ = par
    @unpack vy, vef, vex, vne = par
    @unpack ηef, ηec, ηer, ηy, ηc, ω_n, ω_u = par

    ∂ec,∂er,∂qc,∂qr = energy_qe_cr(par)
    ∂ef,∂ef_d = energy_ef(par,∂qc, ∂qr)
    ∂y, ∂y_d = output_y(par, ∂qc, ∂qr)

    dqf, dEf,dEc, dϵ̄1 = totalenergy_qef(par, ∂ef, ∂ec, dϵ̄, dp, dtϵ, dsϵ)

    dp, dp_y,dp_e,dp_t, 𝐭bc_ij = tradeblock_p(par, ∂y, ∂ef, dϵ̄, dqf, dtϵ, dsϵ, dtb, 𝕁)

    
    dD =  .- γ̄ .* dϵ̄
    der = ∂er[:,1].*dD .+ ∂er[:,2].*dp .+ ∂er[:,3].*dqf .+ ∂er[:,4].*dtϵ .+ ∂ec[:,5].*dsϵ
    dec = ∂ec[:,1].*dD .+ ∂ec[:,2].*dp .+ ∂ec[:,3].*dqf .+ ∂ec[:,4].*dtϵ .+ ∂ec[:,5].*dsϵ
    def = ∂ef[:,1].*dD .+ ∂ef[:,2].*dp .+ ∂ef[:,3].*dqf .+ ∂ef[:,4].*dtϵ .+ ∂ef[:,5].*dsϵ
    dy  = ∂y[:,1].*dD  .+ ∂y[:,2].*dp  .+ ∂y[:,3].*dqf .+ ∂y[:,4].*dtϵ .+ ∂y[:,5].*dsϵ

    dex = (1 ./νf).*dqf .- dp 
    ∂ex_p = .- ones(Ntot)

    dqc = ∂qc[:,1].*dD .+ ∂qc[:,2].*dp .+ ∂qc[:,3].*dqf .+ ∂qc[:,4].*dtϵ .+ ∂qc[:,5].*dsϵ
    dqr = ∂qr[:,1].*dD .+ ∂qr[:,2].*dp .+ ∂qr[:,3].*dqf .+ ∂qr[:,4].*dtϵ .+ ∂qr[:,5].*dsϵ

    dP = Smat * dp +  (Smat .* 𝐭bc_ij) * ones(Ntot) 
    dP_r = Smat * dp 

    @unpack νf, νc, νr = par
    dπf = (1 .+ 1 ./νf).*dqf .- (1 ./νf).*dp
    dπc = (1 .+ 1 ./νc).*dqc .- (1 ./νc).*dp
    dπr = (1 .+ 1 ./νr).*dqr .- (1 ./νr).*dp
    
    dW_Dy = (ηy ./ ηc).* dD
    dW_p = (ηy ./ ηc).* dp
    dW_ef = .- (ηy ./ ηc).*sh_E .* sh_ef .* dqf
    dW_ec = .- (ηy ./ ηc).*sh_E .* sh_ec .* dqc
    dW_er = .- (ηy ./ ηc).*sh_E .* sh_er .* dqr
    dW_πf = (ηef./ ηc) .* dπf
    dW_πc = (ηec./ ηc) .* dπc
    dW_πr = (ηer./ ηc) .* dπr
    dW_P  = .- dP_r 
    dW_Du  = dD 

    dWtot = dW_Dy + dW_p + dW_P + dW_ef+dW_ec+dW_er+dW_πf+dW_πc+dW_πr

    dW = hcat(dWtot,dW_Dy, dW_Du, dW_p, dW_P, dW_ef, dW_ec,dW_er, dW_πf, dW_πc, dW_πr);

    dWdec = hcat(dWtot, dW_Dy, dW_p + dW_P, dW_ef+dW_ec+dW_er, dW_πf+dW_πc+dW_πr);

    dW_wld_u = ω_u' * dWtot
    dW_wld_n = ω_n' * dWtot

    # World Welfare decomposition 
    dW_Dy_wld_u = ω_u' * dW_Dy
    dW_Dy_wld_n = ω_n' * dW_Dy

    dW_p_wld_u = ω_u' * (dW_p + dW_P)
    dW_p_wld_n = ω_n' * (dW_p + dW_P)

    dW_e_wld_u = ω_u' * (dW_ef+dW_ec+dW_er)
    dW_e_wld_n = ω_n' * (dW_ef+dW_ec+dW_er)

    dW_π_wld_u = ω_u' * (dW_πf+dW_πc+dW_πr)
    dW_π_wld_n = ω_n' * (dW_πf+dW_πc+dW_πr)

    C = Tmat_y.* vy' .- I(Ntot)

    Ap = θ.*I(Ntot) .- Tmat_y.* vy' .- (θ-1)*(Tmat_y*Smat) .-  Tmat_y .* (vex .* ∂ex_p .-vef.*∂ef[:,2])' .- C .*(∂y[:,2]')

    A =inv(Ap) 

    m = Model_suffstats_o(par = par) ; 

    m.dϵ = dϵ̄
    m.dϵ_eq = dϵ̄1
    m.dtϵ = dtϵ
    m.dsϵ = dsϵ
    m.dtb = 𝐭bc_ij
    m.𝕁 = 𝕁
    m.dD = dD
    m.dp = dp
    m.dP = dP
    m.dP_r = dP_r
    m.dy = dy
    m.dqf = dqf
    m.dqc = dqc
    m.dqr = dqr
    m.def = def
    m.dec = dec
    m.der = der
    m.dπf = dπf
    m.dπc = dπc
    m.dπr = dπr
    m.dex = dex
    m.dEf = dEf
    m.dEc = dEc  

    m.A = A 

    # welfare
    m.dW_Dy = dW_Dy
    m.dW_Du = dW_Du
    m.dW_p = dW_p
    m.dW_P = dW_P
    m.dW_ef = dW_ef 
    m.dW_ec = dW_ec
    m.dW_er = dW_er
    m.dW_πf = dW_πf
    m.dW_πc = dW_πc
    m.dW_πr = dW_πr

    m.dWtot = dWtot
    m.dWdec = dWdec
    m.dW = dW
    m.dW_wld_u = dW_wld_u
    m.dW_wld_n = dW_wld_n
    
    dx =hcat(dtϵ,dsϵ,𝕁,dD ,dp ,dP ,dy ,dqc,dqr,def,dec,der,dπf,dπc,dπr,dex,dW_Dy,dW_Du,dW_p ,dW_P ,dW_ef,dW_ec,dW_er,dW_πf,dW_πc,dW_πr,dWtot)
    dX = hcat(dϵ̄,dϵ̄1,dqf,dEf,dEc,   dW_wld_u,dW_Dy_wld_u,dW_p_wld_u,dW_e_wld_u,dW_π_wld_u ,  dW_wld_n, dW_Dy_wld_n,dW_p_wld_n,dW_e_wld_n,dW_π_wld_n )

    m.dx = dx 
    m.dX = dX
    return m
end 

# solution wrapper
function solve_GE(par::Params_IAMHAT_suffstats_o, dϵ̄, dtϵ, dsϵ, dtb, 𝕁, resid)

    Ntot = par.Ntot
    function optim_GE!(F,x)
        if resid
            px = x[1:Ntot]
            ϵx = x[Ntot+1]
        else 
            px = x[1:Ntot]
            ϵx = dϵ̄
        end 
        F .= GE_fct(par,px, ϵx, dtϵ, dsϵ, dtb, 𝕁, resid) ; 

    end 
    if resid
        eq_init = zeros(Ntot+1) 
    else 
        eq_init = zeros(Ntot) 
    end 

    ss_eq = nlsolve(optim_GE!, eq_init, ftol = 1e-8,  autodiff = :forward)
    
    if resid
        dp_eq = ss_eq.zero[1:Ntot]
        dϵ_eq = ss_eq.zero[Ntot+1]
    else 
        dp_eq = ss_eq.zero
        dϵ_eq = dϵ̄; 
    end 

    print("Solution: ",  ss_eq.f_converged, " Iterations: ",ss_eq.iterations, "\n") ; 

    return dp_eq,dϵ_eq 
end

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --- counterfactuals ------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

####################################
###### Climate change shock ########
####################################

# set up number of countries
Ntot = 193; 

# generate carbon impulse
dϵ̄1 = 1.00 ; 
dp0 = zeros(Ntot);
dtϵ0 = zeros(Ntot);

# set policy instruments to zero 
dsϵ0 = zeros(Ntot);
dtb0 = 0.0;

# get GE price / carbon responses 
dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot), false) ; 

# solve for counterfactual
m1 = GE_solution(parss,dp1, dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot));

# sanity check output 
bar([m1.dWtot[ind_ls]; m1.dW_wld_u], xticks = (1:26,iso3_25w), xrotation=30,title="Total welfare")

# welfare output 
dWdec1 = m1.dWdec  
dWdec1w = m1.dX[:,6:10]
dWdec1w_n = m1.dX[:,11:15]

# save output for 1% shock
df1 = DataFrame([dWdec1 ;dWdec1w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df1.iso3 =  [iso3_l;"WLD"] 

CSV.write("data/model_output/welfare_climate_l_1pc.csv", df1)

# save output for 12% carbon impulse 
df12 = DataFrame(12*[dWdec1 ;dWdec1w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df12.iso3 = [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_climate_l_12pc.csv", df12)

df12_x = DataFrame(12*m1.dx,:auto)
rename!(df12_x, m1.name_varx)
df12_x.iso3 = iso3_l  

CSV.write("data/model_output/welfare_climate_l_12pc_dx.csv", df12_x)

df12_X = DataFrame(12*m1.dX,:auto)
rename!(df12_X, m1.name_varX)

CSV.write("data/model_output/welfare_climate_l_12pc_dXagg.csv", df12_X)


#
# Compute the social cost of carbon #
#
# How to get to 3•c = 
# 0.7 (time horizon) *0.12 (12% impulse) * 20ºc (world temperature (pop weighted) in 2010, i.e. at 1.3ºc) 
# corresponds to 75 * EmEm * 12/100 
#(1+0.7*0.12)*20 - (20-1.3) = 3•c

### UTILITARIAN-weighted

# first calculation: from world welfare cost:
SCC = - sum(1e9 .* parss.pp .* (parss.yy).^(1-parss.ηcrra) ) .* dWdec1w[1] .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

# second calculation: from country welfare cost -- numerically the same
SCC_sum∆U = - sum(1e9 * parss.pp .* (parss.yy).^(1-parss.ηcrra) .* dWdec1[:,1] ) .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

# third calculation: from country LCC_i 

## LCC  = −ci ∆E Ui -- aligns!
LCC_i = - 1e6 .* parss.pp .*  1e3 .* (parss.yy) .* dWdec1[:,1] .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

SCC_sumLCC = sum((parss.yy).^(-parss.ηcrra) .* LCC_i)

println("utilitarian SCC: $SCC_sumLCC")

open(outfile * "SCC_ests.txt", "a") do f
    write(f, "Utilitarian SCC: $SCC_sumLCC\n")
end

### NEGISHI-weighted

𝝎_n =  parss.yy.^(parss.ηcrra) ./ sum(parss.p̂p̂ .* parss.yy.^(parss.ηcrra) )

SCC = - sum(1e9 .* 𝝎_n .* parss.pp .* (parss.yy).^(1-parss.ηcrra) ) .* dWdec1w[1] .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

SCC_sum∆U = - sum(1e9 .* 𝝎_n .* parss.pp .* (parss.yy).^(1-parss.ηcrra) .* dWdec1[:,1] ) .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

# LCC  = −ci ∆E Ui
LCC_i = - 1e6 .* parss.pp .*  1e3 .* (parss.yy) .* dWdec1[:,1] .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

SCC_sumLCC = sum(𝝎_n .* (parss.yy).^(-parss.ηcrra) .* LCC_i)

println("Negishi weighted SCC: $SCC_sumLCC")

open(outfile * "SCC_ests.txt", "a") do f
    write(f, "Negishi weighted SCC: $SCC_sumLCC\n")
end

## Account for non linearity in Climate System (footnote exercise) ##

# change the point of linear approximation
parnl = deepcopy(parss) ; 
parnl.tt = parnl.tt .+= 1.7 
parnl.γ̄  = 2 .* parnl.γdam .* (parnl.tt - parnl.t̄) .* parnl.tt .* parnl.s_ES

# recompute suff stat experiment
dϵ̄1nl = 1.00 ; 
dp0nl = zeros(Ntot);
dtϵ0nl = zeros(Ntot);
#dtϵ0[1] = 1.0

dsϵ0nl = zeros(Ntot);
dtb0nl = 0.0;
dp1nl,dϵ̄2nl = solve_GE(parnl,dϵ̄1nl, dtϵ0nl, dsϵ0nl, dtb0nl, zeros(Ntot), false) ; 

m1nl = GE_solution(parnl,dp1nl, dϵ̄1nl, dtϵ0nl, dsϵ0nl, dtb0nl, zeros(Ntot));

dWdec1nl = m1nl.dWdec  
dWdec1wnl = m1nl.dX[:,6:10]

# first calculation: from world welfare cost:
SCC_nl = - sum(1e9 .* parss.pp .* (parss.yy).^(1-parss.ηcrra) ) .* dWdec1wnl[1] .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

# second calculation: from country welfare cost:
SCC_sum∆U_nl = - sum(1e9 * parss.pp .* (parss.yy).^(1-parss.ηcrra) .* dWdec1nl[:,1] ) .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

# third calculation: from country LCC_i 

# LCC  = −ci ∆E Ui
LCC_i_nl = - 1e6 .* parss.pp .*  1e3 .* (parss.yy) .* dWdec1nl[:,1] .* 12 /100 ./ (parss.EmEm .* 1e6 *12/100)

SCC_sumLCC_nl = sum((parss.yy).^(-parss.ηcrra) .* LCC_i_nl)

open(outfile * "SCC_ests.txt", "a") do f
    write(f, "SCC after 3 degrees of warming: $SCC_sumLCC_nl\n")
end

####################################
########### Carbon tax US ##########
####################################

print("USA: ", (1:193)'*(iso3_l.=="USA"))

dϵ̄1 = 0.00
dtϵ0 = zeros(Ntot)
dtϵ0[182] = 1.57
dsϵ0 = zeros(Ntot)
dtb0 = 0.0
dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot), true)

m2 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb0, zeros(Ntot));

dWdec2 = m2.dWdec  
dWdec2w = m2.dX[:,6:10]

df2_us = DataFrame([dWdec2 ;dWdec2w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df2_us.iso3 = [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_carbontaxUS_l_50.csv", df2_us)

df2_us_x = DataFrame(m2.dx,:auto)
rename!(df2_us_x, m2.name_varx)
df2_us_x.iso3 = iso3_l 

CSV.write("data/model_output/welfare_carbontaxUS_l_50_dx.csv", df2_us_x)

df2_us_X = DataFrame(m2.dX,:auto)
rename!(df2_us_X, m2.name_varX)

CSV.write("data/model_output/welfare_carbontaxUS_l_50_dXagg.csv", df2_us_X)



####################################
######## Carbon tax China ##########
####################################

print("CHN: ", (1:193)'*(iso3_l.=="CHN"))

# carbon tax of 1.57 as a ad valorem tax of 1.57*0.23*qef 
# with qef = 31.8 per MWh of oil-gas, implies a carbon tax of 50$/tCO2 
# since oil-gas has 

dϵ̄1 = 0.00
dtϵ0 = zeros(Ntot)
dtϵ0[35] = 1.57
dsϵ0 = zeros(Ntot)
dtb0 = 0.0
dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot), true)

m2 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb0, zeros(Ntot));

dWdec2 = m2.dWdec  
dWdec2w = m2.dX[:,6:10]

# save

df2_chn = DataFrame([dWdec2;dWdec2w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df2_chn.iso3 =  [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_carbontaxCHN_l_50.csv", df2_chn)

df2_chn_x = DataFrame(m2.dx,:auto)
rename!(df2_chn_x, m2.name_varx)
df2_chn_x.iso3 = iso3_l  # replace `rand(10)` with your data

CSV.write("data/model_output/welfare_carbontaxCHN_l_50_dx.csv", df2_chn_x)

df2_chn_X = DataFrame(m2.dX,:auto)
rename!(df2_chn_X, m2.name_varX)

CSV.write("data/model_output/welfare_carbontaxCHN_l_50_dXagg.csv", df2_chn_X)


####################################
##### Unilateral Carbon tax  #######
####################################

# preallocate output
dWdec4_i = zeros(193,5)
dWdec4_w = zeros(193,5)
dXagg4_w = zeros(193,15)

for i = 1:193

    dϵ̄1 = 0.00
    dtϵ0 = zeros(Ntot)
    dtϵ0[i] = 1.57
    dsϵ0 = zeros(Ntot)
    dtb0 = 0.0
    dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot), true)

    m4 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb0, zeros(Ntot));

    dWdec4_i[i,:] = m4.dWdec[i,:]  
    dWdec4_w[i,:] = m4.dX[:,6:10][:]
    dXagg4_w[i,:] = m4.dX[:]

end 

df4_i = DataFrame(dWdec4_i,[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df4_i.iso3 = iso3_l  

CSV.write("data/model_output/welfare_carbontax_unilateral_i_l_50.csv", df4_i)


df4_w = DataFrame(dWdec4_w,[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df4_w.iso3 = iso3_l  

CSV.write("data/model_output/welfare_carbontax_unilateral_w_l_50.csv", df4_w)

df4_X = DataFrame(dXagg4_w,:auto)
rename!(df4_X, m1.name_varX)

CSV.write("data/model_output/welfare_carbontax_unilateral_X_l_50_dXagg.csv", df4_X)

###########################
### Unilateral subsidy ####
###########################

dWdec5_i = zeros(193,5)
dWdec5_w = zeros(193,5)
dXagg5_w = zeros(193,15)

for i = 1:193

    dϵ̄1 = 0.00
    dtϵ0 = zeros(Ntot)
    dsϵ0 = zeros(Ntot)
    dsϵ0[i] = 0.426
    dtb0 = 0.0
    dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot), true)

    m5 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb0, zeros(Ntot));

    dWdec5_i[i,:] = m5.dWdec[i,:]  
    dWdec5_w[i,:] = m5.dX[:,6:10][:]
    dXagg5_w[i,:] = m5.dX[:]

end 

df5_i = DataFrame(dWdec5_i,[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df5_i.iso3 = iso3_l  

CSV.write("data/model_output/welfare_renewsubsidy_unilateral_i_l_50.csv", df5_i)


df5_w = DataFrame(dWdec5_w,[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df5_w.iso3 = iso3_l  
CSV.write("data/model_output/welfare_renewsubsidy_unilateral_w_l_50.csv", df5_w)

df5_X = DataFrame(dXagg5_w,:auto)
rename!(df5_X, m2.name_varX)

CSV.write("data/model_output/welfare_renewsubsidy_unilateral_X_l_50_dXagg.csv", df5_X)

##########################
## US renewable subsidy ##

dϵ̄1 = 0.00
dtϵ0 = zeros(Ntot)
dsϵ0 = zeros(Ntot)
dsϵ0[182] = 0.426
dtb0 = 0.0
dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot), true)

m5 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb0, zeros(Ntot));

dWdec5 = m5.dWdec
dWdec5w = m5.dX[:,6:10]

##################################
############# US CBAM ############

dϵ̄1 = 0.00
dtϵ0 = zeros(Ntot)
dsϵ0 = zeros(Ntot)
dtb0 = 50
J0 = zeros(Ntot)
J0[182] = 1 
dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0,J0 , true)

m6 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb0, J0);

dWdec6 = m6.dWdec
dWdec6w = m6.dX[:,6:10]

df6_us = DataFrame([dWdec6 ;dWdec6w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df6_us.iso3 = [iso3_l;"WLD"]  
CSV.write("data/model_output/welfare_carbontariff_US_l_50.csv", df6_us)

df6_us_x = DataFrame(m6.dx,:auto)
rename!(df6_us_x, m2.name_varx)
df6_us_x.iso3 = iso3_l 

CSV.write("data/model_output/welfare_carbontariff_l_50_dx.csv", df6_us_x)

df6_us_X = DataFrame(m6.dX,:auto)
rename!(df6_us_X, m2.name_varX)

CSV.write("data/model_output/welfare_carbontariff_l_50_dXagg.csv", df6_us_X)

############################################
###### Unilateral Carbon tax + CBAM ########
############################################

dWdec10_i = zeros(193,5)
dWdec10_w = zeros(193,5)
dXagg10_w = zeros(193,15)

for i = 1:193

    dϵ̄1 = 0.00
    dtϵ1 = zeros(Ntot)
    dtϵ1[i] = 1.57
    dsϵ0 = zeros(Ntot)
    dtb1 = 50
    J1 = zeros(Ntot)
    J1[i] = 1.0 

    dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ1, dsϵ0, dtb0, J1, true)

    m10 = GE_solution(parss,dp1, dϵ̄2, dtϵ1, dsϵ0, dtb0, J1);

    dWdec10_i[i,:] = m10.dWdec[i,:]  
    dWdec10_w[i,:] = m10.dX[:,6:10][:]
    dXagg10_w[i,:] = m10.dX[:]

end 

df10_i = DataFrame(dWdec10_i,[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df10_i.iso3 = iso3_l  

CSV.write("data/model_output/welfare_carbontax_carbontariffs_unilateral_i_l_50.csv", df10_i)

df10_w = DataFrame(dWdec10_i,[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df10_w.iso3 = iso3_l  

CSV.write("data/model_output/welfare_carbontax_carbontariffs_unilateral_w_l_50.csv", df10_w)

df10_X = DataFrame(dXagg10_w,:auto)
rename!(df10_X, m1.name_varX)

CSV.write("data/model_output/welfare_carbontax_carbontariffs_unilateral_X_l_50_dXagg.csv", df10_X)

########################################
############# EU carbon tax ############

dϵ̄1 = 0.00
dtϵ1 = zeros(Ntot)
dtϵ1[is_eu] .= 1.57
dsϵ0 = zeros(Ntot)
dtb0 = 0
Jeu_fl = Float64.(is_eu)
J0 = zeros(Ntot)

dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ1, dsϵ0, dtb0,J0 , true)

m7 = GE_solution(parss,dp1, dϵ̄2, dtϵ1, dsϵ0, dtb0, J0);

dWdec7 = m7.dWdec
dWdec7w = m7.dX[:,6:10]

df7_eu = DataFrame([dWdec7 ;dWdec7w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df7_eu.iso3 = [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_carbontax_EU_l_50.csv", df7_eu)

df7_eu_x = DataFrame(m7.dx,:auto)
rename!(df7_eu_x, m7.name_varx)
df7_eu_x.iso3 = iso3_l 

CSV.write("data/model_output/welfare_carbontax_EU_l_50_dx.csv", df7_eu_x)

df7_eu_X = DataFrame(m7.dX,:auto)
rename!(df7_eu_X, m7.name_varX)

CSV.write("data/model_output/welfare_carbontax_EU_l_50_dXagg.csv", df7_eu_X)

########################################
############# EU CBAM alone ############

dϵ̄1 = 0.00
dtϵ0 = zeros(Ntot)
dsϵ0 = zeros(Ntot)
dtb1 = 50
Jeu_fl = Float64.(is_eu)
J0 = zeros(Ntot)

dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb1,Jeu_fl , true)

m9 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb1, Jeu_fl);

dWdec9 = m9.dWdec
dWdec9w = m9.dX[:,6:10]

df9_eu = DataFrame([dWdec9 ;dWdec9w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df9_eu.iso3 = [iso3_l;"WLD"]  
CSV.write("data/model_output/welfare_carbontariff_EU_l_50.csv", df9_eu)

df9_eu_x = DataFrame(m9.dx,:auto)
rename!(df9_eu_x, m9.name_varx)
df9_eu_x.iso3 = iso3_l 
CSV.write("data/model_output/welfare_carbontariff_EU_l_50_dx.csv", df9_eu_x)

df9_eu_X = DataFrame(m9.dX,:auto)
rename!(df9_eu_X, m9.name_varX)
CSV.write("data/model_output/welfare_carbontariff_EU_l_50_dXagg.csv", df9_eu_X)

###############################################
############# EU carbon tax + CBAM ############

dϵ̄1 = 0.00
dtϵ1 = zeros(Ntot)
dtϵ1[is_eu] .= 1.57
dsϵ0 = zeros(Ntot)
dtb1 = 50
Jeu_fl = Float64.(is_eu)
J0 = zeros(Ntot)

dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ1, dsϵ0, dtb1,Jeu_fl , true)

m8 = GE_solution(parss,dp1, dϵ̄2, dtϵ1, dsϵ0, dtb1, Jeu_fl);

dWdec8 = m8.dWdec
dWdec8w = m8.dX[:,6:10]

df8_eu = DataFrame([dWdec8 ;dWdec8w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df8_eu.iso3 = [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_carbontax_carbontariff_EU_l_50.csv", df8_eu)

df8_eu_x = DataFrame(m8.dx,:auto)
rename!(df8_eu_x, m8.name_varx)
df8_eu_x.iso3 = iso3_l 

CSV.write("data/model_output/welfare_carbontax_carbontariff_EU_l_50_dx.csv", df8_eu_x)

df8_eu_X = DataFrame(m8.dX,:auto)
rename!(df8_eu_X, m8.name_varX)

CSV.write("data/model_output/welfare_carbontax_carbontariff_EU_l_50_dXagg.csv", df8_eu_X)

#####################################
############# Asean CBAM ############

dϵ̄1 = 0.00
dtϵ0 = zeros(Ntot)
dsϵ0 = zeros(Ntot)
dtb1 = 50
Jasean_fl = Float64.(is_asean)
J0 = zeros(Ntot)

dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb1,Jasean_fl , true)

m11 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb1, Jasean_fl);

dWdec11 = m11.dWdec
dWdec11w = m11.dX[:,6:10]

df11_as = DataFrame([dWdec11 ;dWdec11w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df11_as.iso3 = [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_carbontariff_As_l_50.csv", df11_as)

df11_as_x = DataFrame(m11.dx,:auto)
rename!(df11_as_x, m11.name_varx)
df11_as_x.iso3 = iso3_l 

CSV.write("data/model_output/welfare_carbontariff_As_l_50_dx.csv", df11_as_x)

df11_as_X = DataFrame(m11.dX,:auto)
rename!(df11_as_X, m11.name_varX)

CSV.write("data/model_output/welfare_carbontariff_As_l_50_dXagg.csv", df11_as_X)

########################################
########## Asean carbon tax ############

dϵ̄1 = 0.00
dtϵ1 = zeros(Ntot)
dtϵ1[is_asean] .= 1.57
dsϵ0 = zeros(Ntot)
dtb0 = 0
Jas_fl = Float64.(is_asean)
J0 = zeros(Ntot)

dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ1, dsϵ0, dtb0,J0 , true)

m12 = GE_solution(parss,dp1, dϵ̄2, dtϵ1, dsϵ0, dtb0, J0);

dWdec12 = m12.dWdec
dWdec12w = m12.dX[:,6:10]

df12_as = DataFrame([dWdec12 ;dWdec12w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df12_as.iso3 = [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_carbontax_As_l_50.csv", df12_as)

df12_as_x = DataFrame(m12.dx,:auto)
rename!(df12_as_x, m12.name_varx)
df12_as_x.iso3 = iso3_l 

CSV.write("data/model_output/welfare_carbontax_As_l_50_dx.csv", df12_as_x)

df12_as_X = DataFrame(m12.dX,:auto)
rename!(df12_as_X, m12.name_varX)

CSV.write("data/model_output/welfare_carbontax_As_l_50_dXagg.csv", df12_as_X)

###############################################
########## Asean carbon tax + CBAM ############

dϵ̄1 = 0.00
dtϵ1 = zeros(Ntot)
dtϵ1[is_asean] .= 1.57
dsϵ0 = zeros(Ntot)
dtb1 = 50
Jas_fl = Float64.(is_asean)
J0 = zeros(Ntot)

dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ1, dsϵ0, dtb1,Jas_fl , true)

m13 = GE_solution(parss,dp1, dϵ̄2, dtϵ1, dsϵ0, dtb1, Jas_fl);

dWdec13 = m13.dWdec
dWdec13w = m13.dX[:,6:10]

df13_as = DataFrame([dWdec13 ;dWdec13w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df13_as.iso3 = [iso3_l;"WLD"]  

CSV.write("data/model_output/welfare_carbontax_carbontariff_As_l_50.csv", df13_as)

df13_as_x = DataFrame(m13.dx,:auto)
rename!(df13_as_x, m13.name_varx)
df13_as_x.iso3 = iso3_l 

CSV.write("data/model_output/welfare_carbontax_carbontariff_As_l_50_dx.csv", df13_as_x)

df13_as_X = DataFrame(m13.dX,:auto)
rename!(df13_as_X, m13.name_varX)

CSV.write("data/model_output/welfare_carbontax_carbontariff_As_l_50_dXagg.csv", df13_as_X)

###########################
### Global carbon tax #####
###########################

dϵ̄1 = 0.00
dtϵ0 = ones(Ntot)
dsϵ0 = 1.57.*zeros(Ntot)
dtb0 = 0.0
dp1,dϵ̄2 = solve_GE(parss,dϵ̄1, dtϵ0, dsϵ0, dtb0, zeros(Ntot), true)

m3 = GE_solution(parss,dp1, dϵ̄2, dtϵ0, dsϵ0, dtb0, zeros(Ntot));

dWdec3 = m3.dWdec  
dWdec3w = m3.dX[:,6:10]

df3_wld = DataFrame([dWdec3;dWdec3w],[:dW, :dW_D, :dW_p, :dW_e,:dW_π])
df3_wld.iso3 = [iso3_l;"WLD"] 

CSV.write("data/model_output/welfare_carbontaxWLD_l_50.csv", df3_wld)

df3_wld_x = DataFrame(m3.dx,:auto)
rename!(df3_wld_x, m3.name_varx)
df3_wld_x.iso3 = iso3_l  

CSV.write("data/model_output/welfare_carbontaxWLD_l_50_dx.csv", df3_wld_x)

df3_wld_X = DataFrame(m3.dX,:auto)
rename!(df3_wld_X, m3.name_varX)

CSV.write("data/model_output/welfare_carbontaxWLD_l_50_dXagg.csv", df3_wld_X)











