module SSA

export decompSSA

function getphascor(U1,U2)
    phas = atan2(U1,U2)
    phas=phas[:]
    for i=1:length(phas)-1
        if phas[i+1]-phas[i] < -pi
            phas[i+1:end]=phas[i+1:end].+2*pi
        elseif phas[i+1]-phas[i] > pi
            phas[i+1:end]=phas[i+1:end].-2*pi
        end
    end
    return(cor(phas,[1:length(phas)]))
end

function decompSSA(v::Array{Float64,1},ncomp::Integer=10,M::Integer=ifloor(length(v)/2))
    
    N=length(v)
    hank=Float64[v[i+j-1] for i=1:M, j=1:(N-M+1)];
    S=A_mul_Bt(hank,hank);
    (U,l,V)=svd(S);
    
    xout = Array(Float64,N,ncomp)
    #Allocate some temporary arrays
    V1=zeros(Float64,N-M+1,1)
    V2=zeros(Float64,N-M+1,1)
    X1=zeros(Float64,M,N-M+1)
    X2=zeros(Float64,M,N-M+1)
    ipair=1
    icomp=1
    while icomp <= ncomp
        U1=sub(U,1:M,(ipair):(ipair))
        U2=sub(U,1:M,(ipair+1):(ipair+1))
        At_mul_B!(V1,hank,U1)
        A_mul_Bt!(X1,U1,V1)
        dopair = abs(getphascor(U1,U2)) > 0.8
        if dopair
            dopair && At_mul_B!(V2,hank,U2)
            A_mul_Bt!(X2,U2,V2)
            for f = 1:length(X1) X1[f]=X1[f]+X2[f] end
        end 
        #Diagonal averaging
        for k=1:N
            j1 = k>M ? k-M+1 : 1
            jn = k>(N-M+1) ? N-M+1 : k
            # Calculate mean
            s  = 0.0
            nc = 0
            for j=j1:jn
                s=s+X1[k-j+1,j]
                nc = nc+1
            end
            xout[k,icomp]=s/nc
        end
        println("$icomp $dopair $(getphascor(U1,U2)) $(var(xout[:,icomp]))")
        ipair = dopair ? ipair+2 : ipair+1
        icomp = icomp+1
    end
    xout
end

end # module
