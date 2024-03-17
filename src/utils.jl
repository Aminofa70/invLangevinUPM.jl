
"""
Ref: # H.R. Warner , Kinetic theory and rheology of dilute suspensions of finitely ex- tendible dumbbells, Ind. Eng. Chem. Fund 11 (1972) 379–387 .
"""
function invWarner()
    function eval_inv(x)
        out = (3*x)/(1-x^2)
        return out
    end
    return eval_inv::Function
end
#################################
"""
Ref:  M. Puso , Mechanistic Constitutive Models for Rubber Elasticity and Viscoelas- ticity, Phd thesis, University of California, Davis, 2003 .

"""
function invPuso()
    function eval_inv(x)
        out = (3*x)/(1-x^3)
        return out
    end
        return eval_inv::Function
end
 
#################################
"""
Ref: Cohen, A. (1991). "A Padé approximant to the inverse Langevin function". Rheologica Acta. 30 (3): 270–273. doi:10.1007/BF00366640
"""
function invCohen()
    function eval_inv(x)
        out = x * ((3 - x^2) / (1 - x^2))
        return out
    end
    return eval_inv::Function
end 
#################################
"""
Ref: Jedynak, R. (2015). "Approximation of the inverse Langevin function revisited". Rheologica Acta. 54 (1): 29–39. doi:10.1007/s00397-014-0802-2.

"""
function invJedynak()
    function eval_inv(x)
        out =(x * (3.0 - 2.6x + 0.7x^2)) / ((1 - x) * (1 + 0.1x))
        return out
    end
    return eval_inv::Function
end 

#################################
"""
Ref: Petrosyan, R. (2016). "Improved approximations for some polymer extension models". Rheologica Acta. 56: 21–26. 
"""
function invPetrosyan()
    function eval_inv(x)
        out = 3x + (x^2/5) * sin(7x/2) + x^3/(1 - x)
        return out
    end
    return eval_inv::Function
end

#################################
"""
Ref: Jedynak, R. (2017). "New facts concerning the approximation of the inverse Langevin function". Journal of Non-Newtonian Fluid Mechanics. 249: 8–25. doi:10.1016/j.jnnfm.2017.09.003
"""

function invJedynak_modified()
    function eval_inv(x)
        out = x * (3 - 1.00651x^2 - 0.962251x^4 + 1.47353x^6 - 0.48953x^8) / ((1 - x) * (1 + 1.01524x))
        return out
    end
    return eval_inv::Function
end

#################################

function invMarchiArruda()
    function eval_inv(y)
        out = (y*((-3 + 2.96295*y - 0.96292*y^2) / (y - 1))) + (0.28701*y)^11.33414 - (1.40114*y^3.42076) * (y - 0.78833)

        return out
    end
    return eval_inv::Function
end
#################################
"""
Ref: Benítez, J.M.; Montáns, F.J. (2018). "A simple and efficient numerical procedure to compute the inverse Langevin function with high accuracy". Journal of Non-Newtonian Fluid Mechanics. 261: 153–163.

"""
function invBenitezMontans()
    np  = 100001.; xir = 0.980; yir = 50.      
    dyir= yir^2/(1+yir^2-yir^2*(coth(yir))^2)  

    y1 = collect(1e-14:(2yir-1e-14)/(2np-2):2yir)  
    x1 = coth.(y1)-1.0./y1         
    x1[1] = 0.0; y1[1] = 0.0       

   
    sp1 = csplinef(x1,y1) 

   
    x2 = collect(0:xir/(np+1):xir)
    dx = x2[2] - x2[1]       
    y2 = sp1.(x2)             
    sp2= csplinef(x2,y2,1,3.0,1,dyir) 

    
    btemp = yir*(1.0 -xir^2)
    aILrf = -2.0*btemp*xir/(1.0-xir^2) + dyir*(1.0-xir^2)
    bILrf = btemp - aILrf*xir

    
    function InvLangevin(xval)
        xval = convert(Float64, xval)
            if xval <= 0.
                yval = 0.
            elseif xval >=1.
                yval = Inf64
            elseif xval < xir
                yval = sp2(xval)
            else
                yval = (aILrf*xval + bILrf) / (1.0 -xval^2)
            end
        return yval::Float64
    end
    return InvLangevin::Function
end
##############################################
##############################################
"""
Ref: Average-chain behavior of isotropic incompressible polymers obtained from macroscopic experimental data. A simple structure-based WYPiWYG model in Julia language
Advances in Engineering Software 130, 41-57

"""
function csplinef(x::Array{Float64,1}, y::Array{Float64,1},
                  EndCond1::Int64 = 2, EC1::Float64 = 0.,
                  EndCond2::Int64 = 2,EC2::Float64 = 0.)

    N = length(x)
    Ny= length(y)
    (N != Ny) && error("Leng of y (=$Ny) no equal to length of x (=$N)")
    (N < 3)   && error("Length of data (=$N) must be > 3")
    for i=2:N
        (x[i] <= x[i-1]) && error("x values must be monotonically increasing")
    end

    
    D=zeros(Float64,N)            
    h=zeros(Float64,N-1)          
    l=zeros(Float64,N)           
    u=zeros(Float64,N)            
    di=zeros(Float64,N)           
    Y =zeros(Float64,N)           
                                     
    a = zeros(Float64,N-1)       
    b = zeros(Float64,N-1)        
    c = zeros(Float64,N-1)       
    d = zeros(Float64,N-1)        

 
    for i=N-1:-1:1
        h[i] = x[i+1] - x[i]       
    end
    for i=N-1:-1:1
        d[i] = (y[i+1]-y[i])/(x[i+1]-x[i])  
    end

    
    if EndCond1==1        
        D[1] = 6.0 * ( d[1] - EC1 )
        u[1] = h[1]
        di[1]= 2.0*h[1]
    elseif EndCond1==2    
        D[1] = EC1
        u[1] = 0.0
        di[1]= 1.0
    else
        error("EndCond1 cannot be $EndCond1")
    end
    if EndCond2==1        
        D[N] = 6.0 * ( EC2 - d[N-1] )
        l[N] = h[N-1]
        di[N]= 2.0*h[N-1]
    elseif EndCond2==2    
        D[N] = EC2
        l[N] = 0.0
        di[N]= 1.0
    else
        error("EndCond2 cannot be $EndCond2")
    end

    
    for i=N-1:-1:2
        D[i] = 6.0 * ( d[i] - d[i-1] )
        l[i] = h[i-1]
        u[i] = h[i]
        di[i]= 2.0 * ( h[i-1] + h[i] )
    end

    
    u[1] /= di[1]
    D[1] /= di[1]
    for i=2:N-1
        z    = di[i] - l[i] * u[i-1]
        u[i] /= z
        D[i] = ( D[i] - l[i] * D[i-1] ) / z
    end
    D[N] = ( D[N] - l[N] * D[N-1] ) / ( di[N] - l[N] * u[N-1] )
    Y[N] = D[N]
    for i=N-1:-1:1
        Y[i] = D[i] - u[i] * Y[i+1]
    end

   
    for i=1:N-1
        a[i]  = y[i]
        b[i]  = d[i] - h[i] / 6.0 * (2.0 * Y[i] + Y[i+1])
        c[i]  = Y[i] / 2.0
        d[i]  = (Y[i+1] - Y[i]) / (6.0 * h[i])
    end

    
    aN = 0.0; if all(h[1] .== h); aN = 1.0; end

   
    C = [[a;aN] [b;y[1]] [c;[y[N]]] [d;0.0] x]

    function csplinefval(xx)
 
        N = size(C)[1]       
        xx = convert(Float64, xx)         
    
        if C[N,1] == 1.
       
            x0 = C[1,5]
            h = C[2,5] - x0
                
                if xx <= x0
                    s = 1
                elseif xx >= x0 + h * (N-1)
                    s = N-1
                else
                    s = floor(Int,(xx - x0)/h) + 1
                end
                
                z = xx - C[s,5]; z2 = z*z
                yy = C[s,1] + C[s,2] * z + C[s,3] * z2 + C[s,4] * z2 * z
        else

                s0 = 1
                s1 = N-1
                if xx <= C[s0,5]
                    s = s0
                elseif xx >= C[s1,5]
                    s = s1
                else
                    s0 = 1; s1 = N-1
                    while s1-s0 > 1
                        si = s0 + floor(Int,(s1-s0)/2)
                        (xx == C[si,5]) && (s0 = si; s1 = si)
                        (xx > C[si,5])  && (s0 = si)
                        (xx < C[si,5])  && (s1 = si)
                    end
                    s = s0
                end
                
                z = xx - C[s,5]; z2 = z*z
                yy = C[s,1] + C[s,2] * z + C[s,3] * z2 + C[s,4] * z2 * z
        end
        return yy
    end # of csplineval
  
    return csplinefval::Function 
end # of csplinef
