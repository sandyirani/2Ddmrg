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

  (j, k) = getCoords(n)
  if (j == length)
    return(0)
  end

  if (iseven(j))
    return(width-k)
  else
    return(k-1)
  end
end

function getSpanningPairs(n)
  numSpan = getNumSpan(n)

end
