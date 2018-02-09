



function sweepFast(m)
  for ii=-N+1:N-1		# if negative, going right to left
    ii == 0 && continue
    i = abs(ii)
    toright = ii > 0

    println("\n i, dir, m = $i, ",toright ? "to right" : "to left"," $m")

    dleft = size(A[i],1)
    alpha = dleft * 2
    dright = size(A[i+1],3)
    beta = 2 * dright
    onesite = eye(2)
    Ham = zeros(alpha*beta,alpha*beta)

    HL = zeros(dleft,2,dleft,2)
    HR = zeros(2,dright,2,dright)
    #  Inefficient implementation:  m^4   Ham construction
    Hspan = zeros(dleft,2,2,dright,dleft,2,2,dright)
    Ileft = eye(dleft)
    Iright = eye(dright)
    if (toright)
      HLRupdate = zeros(alpha,alpha)
    else
      HLRupdate = zeros(beta,beta)
    end

    #Hamiltonian terms with one endpoint = i (except current edge)
    neighs = getNeighbors(i,i+1)
    for j=1:length(neighs)
      if (neighs[j] < i)
        Aleft = Aopen[i-1,i-1-neighs[j]+1]
        @tensor begin
          HL[a,si,ap,sip] := Htwosite[sl,si,slp,sip] * Aleft[a,sl,ap,slp]
        end
        HL = reshape(HL,alpha,alpha)
        if (toright)
          HLRupdate += HL
        end
        Ham += JK(HL,eye(beta))
      end

      if (neighs[j] > i+1)
        Aright = Aopen[i+2,neighs[j]-i-1]
        @tensor begin
          Hspan[a,si,sip1,b,ap,sip,sip1p,bp] := Htwosite[si,sr,sip,srp] * Aright[b,sr,bp,srp] * onesite[sip1, sip1p] * Ileft[a,ap]
        end
        Ham += reshape(Hspan,alpha*beta,alpha*beta)
      end
    end
    i > 2 && ( Ham += JK(JK(HLR[i-1],onesite),eye(beta)) )
    if (i > 2 && toright)
      HLRupdate += JK(HLR[i-1],onesite)
    end

    #Hamiltonian terms with one endpoint = i+1 (except current edge)
    neighs = getNeighbors(i+1,i)
    for j=1:length(neighs)
      if (neighs[j] > i+1)
        Aright = Aopen[i+2, neighs[j]-i-1]
        @tensor begin
          HR[sip1,b,sip1p,bp] := Htwosite[sip1,sr,sip1p,srp] * Aright[b,sr,bp,srp]
        end
        HR = reshape(HR,beta,beta)
        if (!toright)
          HLRupdate += HR
        end
        Ham += JK(eye(alpha),HR)
      end

      if (neighs[j] < i)
        Aleft = Aopen[i-1,i-1-neighs[j]+1]
        @tensor begin
          Hspan[a,si,sip1,b,ap,sip,sip1p,bp] := Htwosite[sl,sip1,slp,sip1p] * Aleft[a,sl,ap,slp] * onesite[si, sip] * Iright[b,bp]
        end
        Ham += reshape(Hspan,alpha*beta,alpha*beta)
      end
    end
    i < N-2 && ( Ham += JK(eye(alpha),JK(onesite,HLR[i+2])) )
    if (i < N-2 && !toright)
      HLRupdate += JK(onesite,HLR[i+2])
    end

    #Hamiltonian terms that span the current edge
    (left, right) = getSpanningPairs(i)
    for j=1:length(left)
      Aright = Aopen[i+2, right[j]-i-1]
      Aleft = Aopen[i-1,i-1-left[j]+1]
      @tensor begin
          Hspan[a,si,sip1,b,ap,sip,sip1p,bp] := Htwosite[sl,sr,slp,srp] * Aleft[a,sl,ap,slp] * Aright[sr,b,srp,bp] * onesite[si,sip] * onesite[sip1, sip1p]
      end
      Ham += reshape(Hspan,alpha*beta,alpha*beta)
    end


    #Hamiltonian term on current edge (i, i+1)
    Oleft =  Any[JK(eye(dleft),sz), 0.5*JK(eye(dleft),sp), 0.5*JK(eye(dleft),sm)]
    Oright = Any[JK(sz,eye(dright)),JK(sm,eye(dright)),JK(sp,eye(dright))]
    for j=1:length(Oleft)
      Ham += JK(reshape(Oleft[j],alpha,alpha),reshape(Oright[j],beta,beta))
    end

    #Current tensor is starting point for Lanczos eig algorithm
    Ai = A[i]
    Ai1 = A[i+1]
    @tensor begin
      AA[a,b,d,e] := Ai[a,b,c] * Ai1[c,d,e]
    end

    bigH = 0.5 * (Ham + Ham')
    evn = eigs(bigH;nev=1, which=:SR,ritzvec=true,v0=reshape(AA,alpha*beta))
    @show evn[1]
    @show size(evn[2])
    gr = evn[2][:,1]

    AA = reshape(gr,dleft,2,2,dright)

    #updates to A, HLR and Aopen
    (A[i],A[i+1],trunc) = dosvd4(AA,m,toright)
    @show trunc
    if toright && i < N-1
      updateToRight(i, HLRupdate)
    elseif !toright && i > 1
      updateToLeft(i, HLRupdate)
    end

  end #iteration = one edge update

end #end of function sweep

function updateToRight(i, HLRupdate)

  (i1,i2,i3) = size(A[i])
  Ai = A[i]
  @tensor begin
    Aopeni[b, si, bp, sip] := Ai[a,si,b]*Ai[a,sip,bp]
  end
  Aopen[i,1] = Aopeni
  for j = 2:min(i,2*width)
    Aopenim1 = Aopen[i-1,j-1]
    @tensor begin
      Aopeni[b,sl,bp,slp] := Aopenim1[a,sl,ap,slp]*Ai[a,si,b]*Ai[ap,si,bp]
    end
    Aopen[i,j] = Aopeni
  end

  Ai2 = reshape(A[i],i1*i2,i3)
  if 1 < i
    @tensor begin
      hlri[b,bp] := HLRupdate[a,ap] * Ai2[a,b] * Ai2[ap,bp]
    end
    HLR[i] = hlri
  end

end

function updateToLeft(i, HLRupdate)

  (i1,i2,i3) = size(A[i+1])
  Ai1 = A[i+1]
  @tensor begin
    Aopeni1[a, si, ap, sip] := Ai1[a,si,b]*Ai1[ap,sip,b]
  end
  Aopen[i+1,1] = Aopeni1
  for j = 2:min(N-i,2*width)
    Aopenip2 = Aopen[i+2,j-1]
    @tensor begin
      Aopenip1[a,sl,ap,slp] := Aopenip2[b,sl,bp,slp]*Ai1[a,si,b]*Ai1[ap,si,bp]
    end
    Aopen[i+1,j] = Aopenip1
  end

  Ai12 = reshape(A[i+1],i1,i2*i3)
  if  i < N-1
    @tensor begin
      hlri1[a,ap] := HLRupdate[b,bp] * Ai12[a,b] * Ai12[ap,bp]
    end
    HLR[i+1] = hlri1
  end

end
