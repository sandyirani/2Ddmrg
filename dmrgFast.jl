



function sweepFast(m)
  energy = 0
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

    #testH = zeros(alpha*beta,alpha*beta)


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
        #testH += JK(leftMats[numPairs],rightMats[numPairs])
      end

      if (neighs[j] > i+1)
        Aright = Aopen[i+2,neighs[j]-i-1]
        for k = 1:lrDim
            numPairs += 1
            hrk = hr[k]
            @tensor begin
                HR[sip1,b,sip1p,bp] := Aright[b,sr,bp,srp] * onesite[sip1, sip1p] * hrk[sr,srp]
            end
            leftMats[numPairs] = JK(Ileft,hl[k])
            rightMats[numPairs] = reshape(HR,beta,beta)
            #testH += JK(leftMats[numPairs],rightMats[numPairs])
        end
      end
    end
    if (i > 2)
        numPairs += 1
        leftMats[numPairs] = JK(HLR[i-1],onesite)
        rightMats[numPairs] = eye(beta)
        #testH += JK(leftMats[numPairs],rightMats[numPairs])
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
        numPairs+=1
        leftMats[numPairs] = eye(alpha)
        rightMats[numPairs] = HR
        #testH += JK(leftMats[numPairs],rightMats[numPairs])
      end

      if (neighs[j] < i)
        Aleft = Aopen[i-1,i-1-neighs[j]+1]
        for k = 1:lrDim
            numPairs += 1
            hlk = hl[k]
            @tensor begin
                HL[a,si,ap,sip] := Aleft[a,sl,ap,slp] * onesite[si,sip] * hlk[sl,slp]
            end
            leftMats[numPairs] = reshape(HL,alpha,alpha)
            rightMats[numPairs] = JK(hr[k],Iright)
            #testH += JK(leftMats[numPairs],rightMats[numPairs])
        end
      end
    end
    if (i < N-2)
        numPairs += 1
        leftMats[numPairs] = eye(alpha)
        rightMats[numPairs] = JK(onesite,HLR[i+2])
        #testH += JK(leftMats[numPairs],rightMats[numPairs])
    end
    if (i < N-2 && !toright)
      HLRupdate += JK(onesite,HLR[i+2])
    end

    #Start updating here

    #Hamiltonian terms that span the current edge
    (left, right) = getSpanningPairs(i)
    for j=1:length(left)
      Aright = Aopen[i+2, right[j]-i-1]
      Aleft = Aopen[i-1,i-1-left[j]+1]
      for k = 1:lrDim
          hlk = hl[k]
          hrk = hr[k]
          numPairs += 1
          @tensor begin
              HL[a,si,ap,sip] := Aleft[a,sl,ap,slp] * onesite[si,sip] * hlk[sl,slp]
          end
          leftMats[numPairs] = reshape(HL,alpha,alpha)
          @tensor begin
              HR[sip1,b,sip1p,bp] := Aright[b,sr,bp,srp] * onesite[sip1, sip1p] * hrk[sr,srp]
          end
          rightMats[numPairs] = reshape(HR,beta,beta)
          #testH += JK(leftMats[numPairs],rightMats[numPairs])
      end
    end


    #Hamiltonian term on current edge (i, i+1)
    for k=1:lrDim
      numPairs += 1
      leftMats[numPairs] = JK(eye(dleft),hl[k])
      rightMats[numPairs] = JK(hr[k],eye(dright))
      #testH += JK(leftMats[numPairs],rightMats[numPairs])
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
    #evn = eigs(testH;nev=1, which=:SR,ritzvec=true,v0=reshape(AA,alpha*beta))
    @show evn[1][1]
    @show size(evn[2])
    gr = evn[2][:,1]
    energy = evn[1][1]

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
  energy

end #end of function sweepFast

function applyH(v::AbstractVector)
  alpha = Int64(params[1])
  beta = Int64(params[2])
  numPairs = Int64(params[3])

  vm = reshape(v, alpha, beta)
  sum = zeros(alpha, beta)
  for j = 1:numPairs
    sum += leftMats[j]' * vm * rightMats[j]
  end
  return(reshape(sum, alpha*beta))
end
