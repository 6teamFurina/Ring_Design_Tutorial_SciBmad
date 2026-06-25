ring = include(joinpath(@__DIR__, "..", "common", "esr-main-18GeV-1IP.jl"))
include(joinpath(@__DIR__, "trombone_utils.jl"))
@kwdef mutable struct Controls
    # These are Union{TPS64,Float64} bc
    # we need to make all TPSs and see how they affect 
    # chromaticity in oder to compensate at end
    dksf1_1::Union{TPS64,Float64} = 0. 
    dksf2_1::Union{TPS64,Float64} = 0. 
    dksd1_1::Union{TPS64,Float64} = 0. 
    dksd2_1::Union{TPS64,Float64} = 0. 
    dksf1_3::Union{TPS64,Float64} = 0. 
    dksf2_3::Union{TPS64,Float64} = 0. 
    dksd1_3::Union{TPS64,Float64} = 0. 
    dksd2_3::Union{TPS64,Float64} = 0. 
    dksf1_5::Union{TPS64,Float64} = 0. 
    dksf2_5::Union{TPS64,Float64} = 0. 
    dksd1_5::Union{TPS64,Float64} = 0. 
    dksd2_5::Union{TPS64,Float64} = 0. 
    dksf1_7::Union{TPS64,Float64} = 0. 
    dksf2_7::Union{TPS64,Float64} = 0. 
    dksd1_7::Union{TPS64,Float64} = 0. 
    dksd2_7::Union{TPS64,Float64} = 0. 
    dksf1_9::Union{TPS64,Float64} = 0. 
    dksf2_9::Union{TPS64,Float64} = 0. 
    dksd1_9::Union{TPS64,Float64} = 0. 
    dksd2_9::Union{TPS64,Float64} = 0. 
    dksf1_11::Union{TPS64,Float64} = 0. 
    dksd1_11::Union{TPS64,Float64} = 0. 

    dksf2_11::Union{TPS64,Float64} = 0. # COMPENSATOR HORIZONTAL CHROM, defined near end
    dksd2_11::Union{TPS64,Float64} = 0. # COMPENSATOR VERTICAL CHROM, defined near end
            
    # Horizontal trombones:
    dnux_mlrf_6::Union{TPS64,Float64} = 0. 
    dnux_mlrr_6::Union{TPS64,Float64} = 0. 
    dnux_ip8::Union{TPS64,Float64}    = 0. 
    dnux_ip10::Union{TPS64,Float64}   = 0. 
    dnux_ip12::Union{TPS64,Float64}   = 0. 
    dnux_ip2::Union{TPS64,Float64}    = 0. 
            
    # Vertical trombones
    dnuy_mlrf_6::Union{TPS64,Float64} = 0. 
    dnuy_mlrr_6::Union{TPS64,Float64} = 0. 
    dnuy_ip8::Union{TPS64,Float64}    = 0. 
    dnuy_ip10::Union{TPS64,Float64}   = 0. 
    dnuy_ip12::Union{TPS64,Float64}   = 0. 
    dnuy_ip2::Union{TPS64,Float64}    = 0. 

    voltage::Union{TPS64,Float64}        = 0. # Cavity voltage which will be turned on/off
end

if !(@isdefined(CONTROLS))
    const CONTROLS = Controls()
    
    # Horizontal Compensator:
    const DNUX_IP4 = DefExpr(()->
        -(CONTROLS.dnux_mlrf_6 + CONTROLS.dnux_mlrr_6 + CONTROLS.dnux_ip8 + CONTROLS.dnux_ip10
            + CONTROLS.dnux_ip12 + CONTROLS.dnux_ip2)
    )
    
    # Vertical Compensator:
    const DNUY_IP4 = DefExpr(()->
        -(CONTROLS.dnuy_mlrf_6 + CONTROLS.dnuy_mlrr_6 + CONTROLS.dnuy_ip8 + CONTROLS.dnuy_ip10
            + CONTROLS.dnuy_ip12 + CONTROLS.dnuy_ip2)
    )
    
    # Global strengths for +1 chromaticity:
    const KSF = 2.3564355393201155
    const KSD = -3.1701524375068932

    # Cavity voltage for 0.05 synchrotron tune
    const VOLTAGE = 3.3478093176533515e6
end
cavities = filter(x->x.kind=="RFCavity", ring.line)
foreach(x->x.voltage = DefExpr(()->CONTROLS.voltage), cavities)

sextupoles = filter(x->x.kind=="Sextupole", ring.line)
arc1_sextupoles = filter(x->occursin(r"_1$", x.name), sextupoles)
arc3_sextupoles = filter(x->occursin(r"_3$", x.name), sextupoles)
arc5_sextupoles = filter(x->occursin(r"_5$", x.name), sextupoles)
arc7_sextupoles = filter(x->occursin(r"_7$", x.name), sextupoles)
arc9_sextupoles = filter(x->occursin(r"_9$", x.name), sextupoles)
arc11_sextupoles = filter(x->occursin(r"_11$", x.name), sextupoles)

# Arc sextupole layouts:
# 1:  D D F
# 3:  F D D
# 5:  D F F
# 7:  F F D
# 9:  D F F
# 11: F D D

# Easiest to separate F and D for each arc and apply according rules:
sfs_1 = filter(x->occursin(r"^sf", x.name), arc1_sextupoles)
sds_1 = filter(x->occursin(r"^sd", x.name), arc1_sextupoles)
sfs_3 = filter(x->occursin(r"^sf", x.name), arc3_sextupoles)
sds_3 = filter(x->occursin(r"^sd", x.name), arc3_sextupoles)
sfs_5 = filter(x->occursin(r"^sf", x.name), arc5_sextupoles)
sds_5 = filter(x->occursin(r"^sd", x.name), arc5_sextupoles)
sfs_7 = filter(x->occursin(r"^sf", x.name), arc7_sextupoles)
sds_7 = filter(x->occursin(r"^sd", x.name), arc7_sextupoles)
sfs_9 = filter(x->occursin(r"^sf", x.name), arc9_sextupoles)
sds_9 = filter(x->occursin(r"^sd", x.name), arc9_sextupoles)
sfs_11 = filter(x->occursin(r"^sf", x.name), arc11_sextupoles)
sds_11 = filter(x->occursin(r"^sd", x.name), arc11_sextupoles)

# 1:  D D F
for i in 1:4:length(sds_1)-3
    sds_1[i+0].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_1)
    sds_1[i+1].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_1)
    sds_1[i+2].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_1)
    sds_1[i+3].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_1)
end

for i in 1:2:length(sfs_1)-1
    sfs_1[i+0].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_1)
    sfs_1[i+1].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_1)
end

# 3:  F D D
for i in 1:4:length(sds_3)-3
    sds_3[i+0].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_3)
    sds_3[i+1].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_3)
    sds_3[i+2].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_3)
    sds_3[i+3].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_3)
end

for i in 1:2:length(sfs_3)-1
    sfs_3[i+0].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_3)
    sfs_3[i+1].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_3)
end

# 5:  D F F
for i in 1:2:length(sds_5)-1
    sds_5[i+0].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_5)
    sds_5[i+1].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_5)
end

for i in 1:4:length(sfs_5)-3
    sfs_5[i+0].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_5)
    sfs_5[i+1].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_5)
    sfs_5[i+2].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_5)
    sfs_5[i+3].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_5)
end

# 7:  F F D
for i in 1:2:length(sds_7)-1
    sds_7[i+0].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_7)
    sds_7[i+1].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_7)
end

for i in 1:4:length(sfs_7)-3
    sfs_7[i+0].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_7)
    sfs_7[i+1].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_7)
    sfs_7[i+2].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_7)
    sfs_7[i+3].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_7)
end

# 9:  D F F
for i in 1:2:length(sds_9)-1
    sds_9[i+0].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_9)
    sds_9[i+1].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_9)
end

for i in 1:4:length(sfs_9)-3
    sfs_9[i+0].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_9)
    sfs_9[i+1].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_9)
    sfs_9[i+2].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_9)
    sfs_9[i+3].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_9)
end

# 11: F D D
for i in 1:4:length(sds_11)-3
    sds_11[i+0].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_11)
    sds_11[i+1].Kn2 = DefExpr(()->KSD+CONTROLS.dksd1_11)
    sds_11[i+2].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_11)
    sds_11[i+3].Kn2 = DefExpr(()->KSD+CONTROLS.dksd2_11)
end

for i in 1:2:length(sfs_11)-1
    sfs_11[i+0].Kn2 = DefExpr(()->KSF+CONTROLS.dksf1_11)
    sfs_11[i+1].Kn2 = DefExpr(()->KSF+CONTROLS.dksf2_11)
end

# Attach phase trombones to marker locations. The map itself is defined in
# trombone_utils.jl so the chapter 12-specific mechanism is easy to inspect.
trombones = [mlrf_6, mlrr_6, ip8, ip10, ip12, ip2, ip4]
tw = twiss(ring)
t = tw.table

attach_trombone!(mlrf_6, ring, t, DefExpr(()->CONTROLS.dnux_mlrf_6), DefExpr(()->CONTROLS.dnuy_mlrf_6))
attach_trombone!(mlrr_6, ring, t, DefExpr(()->CONTROLS.dnux_mlrr_6), DefExpr(()->CONTROLS.dnuy_mlrr_6))
attach_trombone!(ip8, ring, t, DefExpr(()->CONTROLS.dnux_ip8), DefExpr(()->CONTROLS.dnuy_ip8))
attach_trombone!(ip10, ring, t, DefExpr(()->CONTROLS.dnux_ip10), DefExpr(()->CONTROLS.dnuy_ip10))
attach_trombone!(ip12, ring, t, DefExpr(()->CONTROLS.dnux_ip12), DefExpr(()->CONTROLS.dnuy_ip12))
attach_trombone!(ip2, ring, t, DefExpr(()->CONTROLS.dnux_ip2), DefExpr(()->CONTROLS.dnuy_ip2))
attach_trombone!(ip4, ring, t, DNUX_IP4, DNUY_IP4)

# Finally we need to define compensator families to keep chromaticity +1
# 6 arcs * 2 families * 2 planes = 24 families total 
# So 22 families free, 2 fix chromaticity
using GTPSA
# Finally we need to define compensator families to keep chromaticity +1
# 6 arcs * 2 families * 2 planes = 24 families total 
# So 22 families free, 2 fix chromaticity
using GTPSA
# First order all, first order parameter, 3rd order cross terms
dchrom = Descriptor([1,1,1,1,1,1], 3, ones(Int, 24), 1)
p = params(dchrom)
CONTROLS.dksf1_1  = p[1]
CONTROLS.dksf2_1  = p[2]
CONTROLS.dksd1_1  = p[3]
CONTROLS.dksd2_1  = p[4]
CONTROLS.dksf1_3  = p[5]
CONTROLS.dksf2_3  = p[6]
CONTROLS.dksd1_3  = p[7]
CONTROLS.dksd2_3  = p[8]
CONTROLS.dksf1_5  = p[9]
CONTROLS.dksf2_5  = p[10]
CONTROLS.dksd1_5  = p[11]
CONTROLS.dksd2_5  = p[12]
CONTROLS.dksf1_7  = p[13]
CONTROLS.dksf2_7  = p[14]
CONTROLS.dksd1_7  = p[15]
CONTROLS.dksd2_7  = p[16]
CONTROLS.dksf1_9  = p[17]
CONTROLS.dksf2_9  = p[18]
CONTROLS.dksd1_9  = p[19]
CONTROLS.dksd2_9  = p[20]
CONTROLS.dksf1_11 = p[21]
CONTROLS.dksd1_11 = p[22]

# Last two are compensators
CONTROLS.dksf2_11 = p[23]
CONTROLS.dksd2_11 = p[24]

tw = twiss(ring; at=[],GTPSA_descriptor=dchrom)
chromx_grad = GTPSA.gradient(par(tw.tunes[1], 6)[[0,0,0,0,0,0,:]], include_params=true)[7:end]
chromy_grad = GTPSA.gradient(par(tw.tunes[2], 6)[[0,0,0,0,0,0,:]], include_params=true)[7:end]
Mfamilychrom = vcat(chromx_grad[1:end-2]', chromy_grad[1:end-2]')
Mselfchrom = vcat(chromx_grad[end-1:end]', chromy_grad[end-1:end]')

# Reset controls back
CONTROLS.dksf1_1  = 0.
CONTROLS.dksf2_1  = 0.
CONTROLS.dksd1_1  = 0.
CONTROLS.dksd2_1  = 0.
CONTROLS.dksf1_3  = 0.
CONTROLS.dksf2_3  = 0.
CONTROLS.dksd1_3  = 0.
CONTROLS.dksd2_3  = 0.
CONTROLS.dksf1_5  = 0.
CONTROLS.dksf2_5  = 0.
CONTROLS.dksd1_5  = 0.
CONTROLS.dksd2_5  = 0.
CONTROLS.dksf1_7  = 0.
CONTROLS.dksf2_7  = 0.
CONTROLS.dksd1_7  = 0.
CONTROLS.dksd2_7  = 0.
CONTROLS.dksf1_9  = 0.
CONTROLS.dksf2_9  = 0.
CONTROLS.dksd1_9  = 0.
CONTROLS.dksd2_9  = 0.
CONTROLS.dksf1_11 = 0.
CONTROLS.dksd1_11 = 0.

# Now we want to set chromaticity compensator arcs
# We have
#      grad       other      compensator knob contribution
# 0 = [2 x 22] * [22 x 1] + [2 x 2] * X
# Equal to
# [2 x 2] X = -grad*other
# solve  X = [2 x 2] \ -grad*other
# We will need to solve
# 11: F D D

# use inv bc could be GTPSA

# Set chromaticity compensator arcs
let Mfamilychrom=Mfamilychrom, Mselfchrom=Mselfchrom
    global const DKSF2_11 = DefExpr(()->
        (inv(Mselfchrom) * -(Mfamilychrom * [
            CONTROLS.dksf1_1 
            CONTROLS.dksf2_1 
            CONTROLS.dksd1_1 
            CONTROLS.dksd2_1 
            CONTROLS.dksf1_3 
            CONTROLS.dksf2_3 
            CONTROLS.dksd1_3 
            CONTROLS.dksd2_3 
            CONTROLS.dksf1_5 
            CONTROLS.dksf2_5 
            CONTROLS.dksd1_5 
            CONTROLS.dksd2_5 
            CONTROLS.dksf1_7 
            CONTROLS.dksf2_7 
            CONTROLS.dksd1_7 
            CONTROLS.dksd2_7 
            CONTROLS.dksf1_9 
            CONTROLS.dksf2_9 
            CONTROLS.dksd1_9 
            CONTROLS.dksd2_9 
            CONTROLS.dksf1_11
            CONTROLS.dksd1_11
        ]))[1]
    )
    global const DKSD2_11 = DefExpr(()->
        (inv(Mselfchrom) * -(Mfamilychrom * [
            CONTROLS.dksf1_1 
            CONTROLS.dksf2_1 
            CONTROLS.dksd1_1 
            CONTROLS.dksd2_1 
            CONTROLS.dksf1_3 
            CONTROLS.dksf2_3 
            CONTROLS.dksd1_3 
            CONTROLS.dksd2_3 
            CONTROLS.dksf1_5 
            CONTROLS.dksf2_5 
            CONTROLS.dksd1_5 
            CONTROLS.dksd2_5 
            CONTROLS.dksf1_7 
            CONTROLS.dksf2_7 
            CONTROLS.dksd1_7 
            CONTROLS.dksd2_7 
            CONTROLS.dksf1_9 
            CONTROLS.dksf2_9 
            CONTROLS.dksd1_9 
            CONTROLS.dksd2_9 
            CONTROLS.dksf1_11
            CONTROLS.dksd1_11
        ]))[2]
    )

    # Set compensating arc
    for i in 1:4:length(sds_11)-3
        sds_11[i+2].Kn2 = KSD+DKSD2_11
        sds_11[i+3].Kn2 = KSD+DKSD2_11
    end

    for i in 1:2:length(sfs_11)-1
        sfs_11[i+1].Kn2 = KSF+DKSF2_11
    end
end


# Chapter 12 W-function optimization solution after iteration 7.
if !(@isdefined(CH12_W_OPTIMIZED_KNOBS_ITER7))
    const CH12_W_OPTIMIZED_KNOBS_ITER7 = (
        OF_5 = -0.698886562705769,
        OD_5 = -1.9941313963248843,
        OF_7 = -0.7056425776560357,
        OD_7 = 1.7937750928102283,
        TROMBONE_X1 = 0.006670065623374273,
        TROMBONE_X2 = 0.7723006173352076,
        TROMBONE_Y1 = -0.900678844996547,
        TROMBONE_Y2 = 0.5398947063177737,
    )
end

function apply_chapter12_w_solution!(knobs=CH12_W_OPTIMIZED_KNOBS_ITER7)
    CONTROLS.dksf1_5 = knobs.OF_5
    CONTROLS.dksf2_5 = -knobs.OF_5
    CONTROLS.dksd1_5 = knobs.OD_5
    CONTROLS.dksd2_5 = -knobs.OD_5
    CONTROLS.dksf1_7 = knobs.OF_7
    CONTROLS.dksf2_7 = -knobs.OF_7
    CONTROLS.dksd1_7 = knobs.OD_7
    CONTROLS.dksd2_7 = -knobs.OD_7

    CONTROLS.dnux_mlrf_6 = knobs.TROMBONE_X1
    CONTROLS.dnux_mlrr_6 = knobs.TROMBONE_X2
    CONTROLS.dnuy_mlrf_6 = knobs.TROMBONE_Y1
    CONTROLS.dnuy_mlrr_6 = knobs.TROMBONE_Y2

    return knobs
end

apply_chapter12_w_solution!()


ring
