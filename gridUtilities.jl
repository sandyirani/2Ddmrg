function getCoords(n)
  if (n < 1 || n > N)
    return(0,0)
  end

  a = div(n-1, width)
  b = mod(n-1, width)
  j = a+1
  if (iseven(a))
    k = b+1
  else
    k = width-b
  end
  return(j, k)
end

function getIndex(j, k)
  if (iseven(j))
    b = width-k
  else
    b = k-1
  end
  return(((j-1)*width+b)+1)
end

function getNeighbors(n, mpsNeigh)
  # n is an index of the MPS. This testFunction
  # returns the neighbors of particle n besides
  # mpsNeigh.

  (j, k) = getCoords(n)
  count = 1 + Int8(j > 1) + Int8(j < length)
  neighbors = zeros(count)
  current = 1

  up = (k == width? 1: k+1)
  neigh = getIndex(j, up)
  if (neigh != mpsNeigh)
    neighbors[current] = neigh
    current+=1
  end

  down = (k == 1? width: k-1)
  neigh = getIndex(j, down)
  if (neigh != mpsNeigh)
    neighbors[current] = neigh
    current+=1
  end

  if (j > 1)
    neigh = getIndex(j-1,k)
    if (neigh != mpsNeigh)
      neighbors[current] = neigh
      current+=1
    end
  end

  if (j < length)
    neigh = getIndex(j+1,k)
    if (neigh != mpsNeigh)
      neighbors[current] = neigh
      current+=1
    end
  end
end

function getNumSpan(n)
  #returns the number of Hamiltonain terms (a, b)
  # where a < n and b > n+1
  # a and b are the indices of the particles, not the 2d grid coordinates

  (j, k) = getCoords(n)

  vert = ( k > 1 && k < width-1? 1: 0)
  if (j == length)
    return(vert)
  end

  if (iseven(j))
    return(vert+width-k)
  else
    return(vert+k-1)
  end
end

function getSpanningPairs(n)
  #returns the set of Hamiltonain terms (a, b)
  # where a < n and b > n+1
  # a and b are the indices of the particles, not the 2d grid coordinates
  # a's are stored in leftPoints and b's in rightPoints
  numSpan = getNumSpan(n)
  leftPoints = zeros(numSpan)
  rightPoints = zeros(numSpan)
  (j, k) = getCoords(n)

  if (j < length)
    (a,b) = (iseven(j)? (k+1,width): (1,k-1))
    for i = a:b
      leftPoints[i-a+1] = getIndex(j,a)
      rightPoints[i-a+1] = getIndex(j+1,a)
    end
  end

  if (k > 1 && k < width-1)
    low = getIndex(j,1)
    high = getIndex(j,width)
    if (high < low) (low, high) = (high,low) end
    leftPoints[numSpan] = low
    rightPoints[numSpan]  high
  end

end
