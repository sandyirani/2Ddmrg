

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
    Ham = zeros(alpha*beta,alpha*beta)

    HL = zeros(dleft,2,dleft,2)
    Hspan = zeros(dleft,2,2,dright,dleft,2,2,dright)
    Ileft = eye(dleft)
    Iright = eye(dright)

    neighs = getNeighbors(i,i+1)
    for j=1:length(neighs)
      if (neighs[j] < i)
        Aleft = Aopen[i-1,i-1-neighs[j]+1]
        @tensor begin
          HL[a,si,ap,sip] := Htwosite[sl,si,slp,sip] * Aleft[a,sl,ap,slp]
        end
        HL = reshape(HL,alpha,alpha)
        Ham += JK(HL,eye(beta))
      end

      if (neighs[j] > i)
        Aright = Aopen[i+2,neighs[j]-i-1]
        @tensor begin
          Hspan[sip1,b,a,si,sip1p,bp,ap,sip] := Htwosite[si,sr,sip,srp] * Aright[b,sr,bp,srp] * onesite(sip1, sip1p) * Iright[a,ap]
        end
        Ham += reshape(Hspan,alpha*beta,alpha*beta)
      end
    end

    i > 2 && ( Ham += JK(JK(onesite,HLR[i-1]),eye(beta)) )


    HR = zeros(2,dright,2,dright)
    if i < n-1
      Ai2 = A[i+2]
      @tensor begin
        HR[si1,b,si1p,bp] := Htwosite[si1,si2,si1p,si2p] * Ai2[b,si2,a] * Ai2[bp,si2p,a]
      end
      i < n-2 && (HR += JK4(onesite,HLR[i+2]) )
      HR = reshape(HR,beta,beta)
      Ham += JK(eye(alpha),HR)
    end


    Oleft =  Any[JK(eye(dleft),sz), 0.5*JK(eye(dleft),sp), 0.5*JK(eye(dleft),sm)]
    Oright = Any[JK(sz,eye(dright)),JK(sm,eye(dright)),JK(sp,eye(dright))]

    Ai = A[i]
    Ai1 = A[i+1]
    @tensor begin
      AA[a,b,d,e] := Ai[a,b,c] * Ai1[c,d,e]
    end

    #  Inefficient implementation:  m^4   Ham construction
    for j=1:length(Oleft)
      Ham += JK(reshape(Oleft[j],alpha,alpha),reshape(Oright[j],beta,beta))
    end



    bigH = 0.5 * (Ham + Ham')
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
