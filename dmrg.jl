function dosvdtrunc(AA,m)		# AA a matrix;  keep at most m states
    (u,d,v) = svd(AA)
    prob = dot(d,d)		# total probability
    mm = min(m,length(d))	# number of states to keep
    d = d[1:mm]			# middle matrix in vector form
    trunc = prob - dot(d,d)
    U = u[:,1:mm]
    V = v[:,1:mm]'
    (U,d,V,trunc)		# AA == U * diagm(d) * V	with error trunc
end

function dosvdleftright(AA,m,toright)
    (U,d,V,trunc) = dosvdtrunc(AA,m)
    if toright
	V = diagm(d) * V
    else
	U = U * diagm(d)
    end
    (U,V,trunc)
end

function dosvd4(AA,m,toright)	# AA is ia * 2 * 2 * ib;  svd down the middle;  return two parts
    ia = size(AA,1)
    ib = size(AA,4)
    AA = reshape(AA,ia*2,2*ib)
    (U,V,trunc) = dosvdleftright(AA,m,toright)
    mm = size(U,2)
    U = reshape(U,ia,2,mm)
    V = reshape(V,mm,2,ib)
    (U,V,trunc)
end

using TensorOperations

function JK(a,b)	# Julia kron,  ordered for julia arrays; returns matrix
    (a1,a2) = size(a)
    (b1,b2) = size(b)
    reshape(Float64[a[i,ip] * b[j,jp] for i=1:a1, j=1:b1, ip=1:a2, jp=1:b2],a1*b1,a2*b2)
end

function JK4(a,b)	# Julia kron,  ordered for julia arrays, return expanded into 4 indices
    (a1,a2) = size(a)
    (b1,b2) = size(b)
    Float64[a[i,ip] * b[j,jp] for i=1:a1, j=1:b1, ip=1:a2, jp=1:b2]
end

sz = Float64[0.5 0; 0 -0.5]
sp = Float64[0 1; 0 0]
sm = sp'
Htwosite = reshape(JK(sz,sz) + 0.5 * JK(sp,sm) + 0.5 * JK(sm,sp),2,2,2,2)
# order for Htwosite is s1, s2, s1p, s2p

n = 28		# exact n=28 energy is -12.2254405486
#  Make initial product state in up down up down up down pattern (Neel state)
# Make first tensor a 1 x 2 x m tensor; and last is m x 2 x 1  (rather than vectors)
A = [zeros(1,2,1) for i=1:n]
for i=1:n
    A[i][1,iseven(i) ? 2 : 1,1] = 1.0
end

HLR = [zeros(1,1) for i=1:n]	# Initialize to avoid errors on firs sweep
m = 3
for swp = 0:10
    m = round(Int64,1.3*m)
    for ii=-n+1:n-1		# if negative, going right to left
	ii == 0 && continue
	i = abs(ii)
	toright = ii > 0

	println("\n sweep, i, dir, m = $swp, $i, ",toright ? "to right" : "to left"," $m")

	dleft = size(A[i],1)
	alpha = dleft * 2
	dright = size(A[i+1],3)
	beta = 2 * dright
	onesite = eye(2)

	HL = zeros(dleft,2,dleft,2)
	HR = zeros(2,dright,2,dright)
	if i > 1
	    Aim1 = A[i-1]
	    @tensor begin
		HL[a,si,ap,sip] := Htwosite[sim1,si,sim1p,sip] * Aim1[b,sim1,a] * Aim1[b,sim1p,ap]
	    end
	    i > 2 && ( HL += JK4(HLR[i-1],onesite) )
	end
	HL = reshape(HL,alpha,alpha)
	if i < n-1
	    Ai2 = A[i+2]
	    @tensor begin
		HR[si1,b,si1p,bp] := Htwosite[si1,si2,si1p,si2p] * Ai2[b,si2,a] * Ai2[bp,si2p,a]
	    end
	    i < n-2 && (HR += JK4(onesite,HLR[i+2]) )
	end
	HR = reshape(HR,beta,beta)

	Oleft =  Any[JK(eye(dleft),sz), 0.5*JK(eye(dleft),sp), 0.5*JK(eye(dleft),sm)]
	Oright = Any[JK(sz,eye(dright)),JK(sm,eye(dright)),JK(sp,eye(dright))]

	Ai = A[i]
	Ai1 = A[i+1]
	@tensor begin
	    AA[a,b,d,e] := Ai[a,b,c] * Ai1[c,d,e]
	end

#  Inefficient implementation:  m^4   Ham construction
	Ham = zeros(alpha*beta,alpha*beta)
	for j=1:length(Oleft)
	    Ham += JK(reshape(Oleft[j],alpha,alpha),reshape(Oright[j],beta,beta))
	end
	if i > 1
	    Ham += JK(HL,eye(beta))
	end
	if i < n-1
	    Ham += JK(eye(alpha),HR)
	end
	bigH = reshape(Ham,alpha*beta,alpha*beta)
	bigH = 0.5 * (bigH + bigH')
	evn = eigs(bigH;nev=1, which=:SR,ritzvec=true,v0=reshape(AA,alpha*beta))
	@show evn[1]
	@show size(evn[2])
	gr = evn[2][:,1]

	AA = reshape(gr,dleft,2,2,dright)

	(A[i],A[i+1],trunc) = dosvd4(AA,m,toright)
	@show trunc
	if toright
	    if 1 < i < n-1
		(i1,i2,i3) = size(A[i])
		Ai2 = reshape(A[i],i1*i2,i3)
		@tensor begin
		    hlri[b,bp] := HL[a,ap] * Ai2[a,b] * Ai2[ap,bp]
		end
		HLR[i] = hlri
	    end
	else
	    if 1 < i < n-1
		(i1,i2,i3) = size(A[i+1])
		Ai12 = reshape(A[i+1],i1,i2*i3)
		@tensor begin
		    hlri1[a,ap] := HR[b,bp] * Ai12[a,b] * Ai12[ap,bp]
		end
		HLR[i+1] = hlri1
	    end
	end
    end
end

function testFunction(n)
  Return(n*n)
end
