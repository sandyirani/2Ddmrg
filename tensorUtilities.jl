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

function dosvdMid(AA)
    (U,d,V) = svd(AA)
    d = sqrt.(d)
    V = diagm(d) * V
	U = U * diagm(d)
    (U,V)
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
