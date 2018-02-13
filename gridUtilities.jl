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
  count = 1 + Int8(j > 1) + Int8(j < len)
  neighbors = zeros(Int64, count)
  current = 1

  up = (k == width? 1: k+1)
  neigh = getIndex(j, up)
  if (neigh != mpsNeigh)
    neighbors[current] = Int64(neigh)
    current+=1
  end

  down = (k == 1? width: k-1)
  neigh = getIndex(j, down)
  if (neigh != mpsNeigh)
    neighbors[current] = Int64(neigh)
    current+=1
  end

  if (j > 1)
    neigh = getIndex(j-1,k)
    if (neigh != mpsNeigh)
      neighbors[current] = Int64(neigh)
      current+=1
    end
  end

  if (j < len)
    neigh = getIndex(j+1,k)
    if (neigh != mpsNeigh)
      neighbors[current] = Int64(neigh)
      current+=1
    end
  end
  neighbors
end

function getNumSpan(n)
  #returns the number of Hamiltonain terms (a, b)
  # where a < n and b > n+1
  # a and b are the indices of the particles, not the 2d grid coordinates

  (j, k) = getCoords(n)

  if ((isodd(j) && (k == 1 || k == width-1)) || (iseven(j) && (k == 2 || k == width)))
      num = width-2
  else
      num = width-1
  end

  if (j == len)
    num -= (width-k)
  end

  if (j == 1)
      num -= max(0,width-k-1)
  end

  return(Int64(num))
end

function getSpanningPairs(n)
  #returns the set of Hamiltonain terms (a, b)
  # where a < n and b > n+1
  # a and b are the indices of the particles, not the 2d grid coordinates
  # a's are stored in leftPoints and b's in rightPoints
  numSpan = getNumSpan(n)
  leftPoints = zeros(Int64, numSpan)
  rightPoints = zeros(Int64, numSpan)
  (j, k) = getCoords(n)
  count = 0

  #horizontal edges to the right
  if (j < len)
    (a,b) = (iseven(j)? (k+1,width): (1,k-1))
    for i = a:b
      count += 1
      leftPoints[count] = Int64(getIndex(j,i))
      rightPoints[count] = Int64(getIndex(j+1,i))
    end
  end

  #horizonal edges to the left
  if (j > 1)
    (a,b) = (iseven(j)? (1,k-2): (k+2,width))
    for i = a:b
      count += 1
      leftPoints[count] = Int64(getIndex(j-1,i))
      rightPoints[count] = Int64(getIndex(j,i))
    end
  end

  #edge that wraps around the cylinder
  if (isodd(j) && k > 1 && k < width-1) || (iseven(j) && k > 2 && k < width)
    count += 1
    low = getIndex(j,1)
    high = getIndex(j,width)
    if (high < low)
        (low, high) = (high,low)
    end
    leftPoints[count] = Int64(low)
    rightPoints[count] = Int64(high)
  end

  if (count!=numSpan)
      @show("Error in caluclating number of spanning pairs")
  end

  (leftPoints, rightPoints)

end

function getLeftDepth(n)
    #returns the smallest index a such that there is a term
    #(a, b), where a <= n and b > n.
    #Note a <= n because (n, n+1) is a term
    #the number of left side terms to store is n-a+1

    a = n

    if (n >= N-1) return(n) end

    neighbors = getNeighbors(n+1, n+2)
    for j = 1:length(neighbors)
        a = min(a, neighbors[j])
    end
    neighbors = getNeighbors(n+2, n+1)
    for j = 1:length(neighbors)
        a = min(a, neighbors[j])
    end

    (left, right) = getSpanningPairs(n+1)
    for j = 1:length(left)
        a = min(a, left[j])
    end

    return(a)

end

function getRightDepth(n)
    #returns the largest index b such that there is a term
    #(a, b), where b >= n and a < n.
    #Note a <= n because (n, n+1) is a term
    #the number of left side terms to store is n-a+1

    b = n
    if (n <= 2) return(n) end

    neighbors = getNeighbors(n-1, n-2)
    for j = 1:length(neighbors)
        b = max(b, neighbors[j])
    end
    neighbors = getNeighbors(n-2, n-1)
    for j = 1:length(neighbors)
        b = max(b, neighbors[j])
    end

    (left, right) = getSpanningPairs(n-2)
    for j = 1:length(left)
        b = max(b, right[j])
    end

    return(b)

end
