



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
    numPairs = 0


    HL = zeros(dleft,2,dleft,2)
    HR = zeros(2,dright,2,dright)
    #  Inefficient implementation:  m^4   Ham construction

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
        numPairs+=1
        leftMats[numPairs] = HL
        rightMats[numPairs] = eye(beta)
      end

      if (neighs[j] > i+1)
        Aright = Aopen[i+2,neighs[j]-i-1]
        for k = 1:lrDim
            numPairs += 1
            leftMats[numPairs] = JK(Ileft,hl[:,k,:])
            @tensor begin
                HR[sip1,b,sip1p,bp] := Aright[b,sr,bp,srp] * onesite[sip1, sip1p] * hr[sr,k,srp]
            end
            rightMats[numPairs] = reshape(HR,beta,beta)
        end
      end
    end
    if (i > 2)
        numPairs += 1
        leftMats[numPairs] = JK(HLR[i-1],onesite)
        rightMats[numPairs] = eye(beta)
    end
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

    #Start updating here

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
    #This is specific to Hiesenberg will need to make generic
    Oleft =  Any[JK(eye(dleft),sz), 0.5*JK(eye(dleft),sp), 0.5*JK(eye(dleft),sm)]
    Oright = Any[JK(sz,eye(dright)),JK(sm,eye(dright)),JK(sp,eye(dright))]
    for j=1:length(Oleft)
      numPairs += 1
      leftMat[numPairs] = reshape(Oleft[j],alpha,alpha
      rightMat[numPairs] = reshape(Oright[j],beta,beta)
    end

    #Current tensor is starting point for Lanczos eig algorithm
    Ai = A[i]
    Ai1 = A[i+1]
    @tensor begin
      AA[a,b,d,e] := Ai[a,b,c] * Ai1[c,d,e]
    end

    params[1] = alpha
    params[2] = beta
    params[3] = numPairs
    bigH = LinearMap(applyH, alpha*beta; ismutating=false)
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

function applyH(v::AbstractVector)
  alpha = params[1]
  beta = params[2]
  numPairs = params[3]

end
