# Integrate (AX+B)/C between a and b
function _integrate_partial10(A::Real, B::Real, C::Real, a::Real, b::Real)
    @assert !iszero(C)
    (b-a)/C * (A*(b+a)/2 + B)
end

# Integrate A/(BX + C) between a and b
function _integrate_partial01(A::Real, B::Real, C::Real, a::Real, b::Real)
    @assert !iszero(B) || !iszero(C)
    NewT = promote_type(typeof(A), typeof(B), typeof(C), typeof(a), typeof(b)) 
    if iszero(A)
        return zero(NewT)
    elseif iszero(B)
        return A/C *(b-a)
    elseif iszero(C)
        @assert !(a≤0≤b) "You want to integrate $A / $B X  over ($a, $b)" 
        return A/B * log(abs(b/a))
    else
        @assert !(a≤-C/B≤b) "You want to integrate $A / ($B X + $C) over ($a, $b)" 
        return A/B* log(abs((B*b + C)/(B*a + C)))
    end   
end

# Integrate (AX + B)/(CX + D) between a and b
function _integrate_partial11(A::Real, B::Real, C::Real, D::Real, a::Real, b::Real)
    @assert !iszero(C) || !iszero(D)
    if iszero(A)
        return _integrate_partial01(B, C, D, a, b)
    elseif iszero(C)
        return _integrate_partial10(A, B, D, a, b)
    else
        return A/C * (b-a)  + (B/C - (A*D)/C^2)*log(abs((C*b+D)/(C*a+D)))
    end
end

# Integrate A/(BX²+CX + D) between a and b
function _integrate_partial02(A::Real, B::Real, C::Real, D::Real, a::Real, b::Real)
    @assert !iszero(C) || !iszero(D)
    NewT = promote_type(typeof(A), typeof(B), typeof(C), typeof(D), typeof(a), typeof(b)) 
    if iszero(A)
        return zero(NewT)
    elseif iszero(B)
        return _integrate_partial01(B, C, D, a, b)
    else

    end
end

# Integrate (AX+B)/(CX² + DX + E) between a and b
function _integrate_partial12(A::Real, B::Real, C::Real, D::Real, E::Real, a::Real, b::Real; enforceNullDelta::Bool = false )
    @assert !iszero(C) || !iszero(D) || !iszero(E)
    NewT = promote_type(typeof(A), typeof(B), typeof(C), typeof(D), typeof(E), typeof(a), typeof(b)) 
    if iszero(A) && iszero(B)
        return zero(NewT)
    elseif iszero(B) && iszero(E)
        return  _integrate_partial01(A, C, D, a, b)
    elseif iszero(C)
        return _integrate_partial11(A, B, D, E, a, b)
    else
        Δ = enforceNullDelta ? zero(NewT) : D^2 - 4*E*C
        if Δ > 0
            r₁ = (-D - sqrt(Δ))/(2*C)
            r₂ = (-D + sqrt(Δ))/(2*C)
            @assert (r₁ > b) || (r₂ < a) || (r₂ > b && r₁ < a) "You want to integrate $A X + $B/ ($C X^2 + $D X + $E) over ($a, $b)"
        elseif Δ == 0
            r₀ = -D/(2*C)
            @assert (r₀ > b) || (r₀ < a) "You want to integrate $A X + $B/ ($C X^2 + $D X + $E) over ($a, $b)"
        end
        C1 = D/(2*C)
        C2 = -Δ/(4*C^2)
        C3 = B - A*C1
        if enforceNullDelta || C2 == 0
            return A/C * log( abs((b+C1) / (a+C1) ) ) + C3/C * (1/(a+C1) - 1/(b+C1))
        elseif C2 > 0
            sqrtC2 = sqrt(C2)
            return A/(2*C) * log( abs( ((b+C1)^2 + C2) / ((a+C1)^2 + C2) ) ) + C3/(C*sqrtC2) * (atan((b + C1)/sqrtC2) - atan((a + C1)/sqrtC2))
        elseif C2 < 0
            sqrtC2 = sqrt(-C2)
            return A/(2*C) * log( abs( ((b+C1)^2 - C2) / ((a+C1)^2 - C2) ) ) + C3/(C*2*sqrtC2) * (log(abs((b+C1-sqrtC2)/(b+C1+sqrtC2))) - log(abs((a+C1-sqrtC2)/(a+C1+sqrtC2))))      
        end
    end
end

# Integrate X^k/(AX + B) between a and b for k≥0 
@inline function _integration_monome_over_deg1(k::Int, A::Real, B::Real, a::Real, b::Real)
    @assert !iszero(A) || !iszero(B)
    @assert k ≥ 0
    if iszero(B)
        if k == 0
            return 1/A * log(abs(b/a))
        else
            return 1/(A * k) * (b^k - a^k)
        end
    elseif iszero(A)
        return 1/(B * (k+1)) * (b^(k+1) - a^(k+1))
    else
        C = B/A
        right =  b^(k+1)*pFq((1,k+1),(k+2,),-b/C)
        left = a^(k+1)*pFq((1,k+1),(k+2,),-a/C)
        return (right-left)/(A*(C*k+C))
    end
end

# Integrate X^k/(AX² + BX + C) between a and b for k≥0 
@inline function _integration_monome_over_deg2(k::Int, A::Real, B::Real, C::Real, a::Real, b::Real; enforceNullDelta::Bool = false)
    @assert !iszero(A) || !iszero(B) || !iszero(C)
    @assert k ≥ 0
    if iszero(A)
        return _integration_monome_over_deg1(k, B, C, a, b)
    elseif k == 0
        return _integrate_partial12(0, 1, A, B, C, a, b; enforceNullDelta = enforceNullDelta)
    elseif iszero(C)
        return _integration_monome_over_deg1(k-1, A, B, a, b)
    else
        Δ = B^2 - 4*A*C
        if enforceNullDelta || Δ == 0
            rac = -B/(2*A)
            right = b^(k+1)*pFq((2,k+1),(k+2,), b/rac) 
            left  = a^(k+1)*pFq((2,k+1),(k+2,), a/rac) 
            return 1/(C * (k+1)) * (right - left)
        elseif Δ > 0
            sqrtΔ = sqrt(Δ)
            var1 = sqrtΔ + B
            var2 = sqrtΔ - B
            right = b^(k+1) * (var1 * pFq((1,k+1),(k+2,), 2*A*b/var2) + var2 * pFq((1,k+1),(k+2,), -2*A*b/var1))
            left  = a^(k+1) * (var1 * pFq((1,k+1),(k+2,), 2*A*a/var2) + var2 * pFq((1,k+1),(k+2,), -2*A*a/var1))
            return (right - left)/(2*C*(k+1)*sqrtΔ)
         
        elseif Δ <0
            sqrtΔ = sqrt(-Δ)
            var1 = sqrtΔ + B*im
            var2 = sqrtΔ - B*im
            right = b^(k+1) * (var1 * pFq((1,k+1),(k+2,), 2*A*b*im/var2) + var2 * pFq((1,k+1),(k+2,), -2*A*b*im/var1))
            left  = a^(k+1) * (var1 * pFq((1,k+1),(k+2,), 2*A*a*im/var2) + var2 * pFq((1,k+1),(k+2,), -2*A*a*im/var1))
            return (right - left)/(2*C*(k+1)*sqrtΔ)
        end
    end
end