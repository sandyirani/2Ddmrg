include("gridUtilities.jl")
include("tensorUtilities.jl")
include("dmrg.jl")

width = 6
len = 10
N = width*len

#Global variables
sz = Float64[0.5 0; 0 -0.5]
sp = Float64[0 1; 0 0]
sm = sp'
Htwosite = reshape(JK(sz,sz) + 0.5 * JK(sp,sm) + 0.5 * JK(sm,sp),2,2,2,2)
(hl, hr) = dosvdMid(JK(sz,sz) + 0.5 * JK(sp,sm) + 0.5 * JK(sm,sp))
# order for Htwosite is s1, s2, s1p, s2p

#  Make initial product state in up down up down up down pattern (Neel state)
# Make first tensor a 1 x 2 x m tensor; and last is m x 2 x 1  (rather than vectors)
A = [zeros(1,2,1) for i=1:N]
for i=1:N
    A[i][1,iseven(i) ? 2 : 1,1] = 1.0
end

HLR = [zeros(1,1) for i=1:N]	# Initialize to avoid errors on firs sweep
Aopen = [zeros(1,2,1,2) for i=1:N,  j=1:2*width]

function mainLoop()
  m = 3
  for swp = 0:10
    m = round(Int64,1.3*m)
    sweep(m)
  end
end
