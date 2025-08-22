using Statistics

# petite régression linéaire (renvoie slope, intercept, SSE)
function linfit(x::AbstractVector, y::AbstractVector)
    μx, μy = mean(x), mean(y)
    vx = sum((x .- μx).^2)
    vx == 0 && return (0.0, mean(y), sum((y .- mean(y)).^2))
    slope = sum((x .- μx) .* (y .- μy)) / vx
    intercept = μy - slope*μx
    ŷ = intercept .+ slope .* x
    sse = sum((y .- ŷ).^2)
    return (slope, intercept, sse)
end

# pente locale par régression glissante sur z = log(y)
function local_slopes(r::AbstractVector, z::AbstractVector; w::Int=15)
    n = length(r)
    s = similar(r, Float64)
    half = w ÷ 2
    for i in 1:n
        a = max(1, i - half)
        b = min(n, i + half)
        s[i], _, _ = linfit(view(r, a:b), view(z, a:b))
    end
    return s
end

# lissage 1D simple (moyenne glissante)
function smooth1d(x; w::Int=11)
    w = isodd(w) ? w : w+1
    n = length(x)
    y = similar(x, Float64)
    half = w ÷ 2
    for i in 1:n
        a = max(1, i - half)
        b = min(n, i + half)
        y[i] = mean(@view x[a:b])
    end
    return y
end

# détection de ruptures sur la pente locale par test de contraste gauche/droite
function detect_changepoints(sl::AbstractVector; L::Int=25, τ::Float64=0.1, minsep::Int=30)
    n = length(sl)
    score = zeros(Float64, n)
    for i in (L+1):(n-L)
        sL = mean(@view sl[i-L:i-1])
        sR = mean(@view sl[i+1:i+L])
        score[i] = abs(sR - sL)
    end
    # non-max suppression + distance minimale
    cps = Int[]
    i = L+1
    while i <= n-L
        if score[i] > τ
            # cherche le maximum local dans une fenêtre ±minsep/2
            a = i
            b = min(n-L, i + minsep)
            j = argmax(@view score[a:b]) + a - 1
            push!(cps, j)
            i = j + minsep
        else
            i += 1
        end
    end
    return cps
end

"""
Segmente une décroissance exponentielle en régimes ~ C_i * exp(-α_i r).

Arguments:
- r: vecteur croissant des abscisses
- y: données (strictement > 0)
Paramètres:
- w: fenêtre pour pente locale
- L: fenêtre pour contraste G/D
- τ: seuil de détection de rupture (sur différence de pente)
- minseg: longueur minimale d’un segment (en points)

Retourne: (alphas, breaks, fitinfo)
- alphas: vecteur des α_i
- breaks: indices de début de chaque segment (ex: [1, b2, b3, ...])
- fitinfo: vecteur de NamedTuples par segment (slope, intercept, sse, r², range)
"""
function segment_exponential(r::AbstractVector, y::AbstractVector;
        w::Int=15, L::Int=25, τ::Float64=0.1, minseg::Int=40, smoothw::Int=11)

    @assert length(r)==length(y)
    @assert all(y .> 0) "Les y doivent être > 0 pour le log."

    z = log.(y)
    # lissage doux (sur z puis sur la pente)
    zsm = smooth1d(z; w=smoothw)
    s  = local_slopes(r, zsm; w=w)
    ssm = smooth1d(s; w=smoothw)

    # changepoints sur la pente
    cps = detect_changepoints(ssm; L=L, τ=τ, minsep=minseg)

    # bornes de segments (inclure extrémités)
    bounds = sort(unique(vcat(1, cps, length(r))))
    # assure minseg
    pruned = [bounds[1]]
    for k in 2:length(bounds)-1
        if (bounds[k] - pruned[end] >= minseg) &&
           (bounds[k+1] - bounds[k] >= minseg)
            push!(pruned, bounds[k])
        end
    end
    push!(pruned, bounds[end])

    # fit par segment et collecte des α_i
    alphas = Float64[]
    info = NamedTuple[]
    starts = pruned[1:end-1]
    stops  = pruned[2:end]
    for (a,b) in zip(starts, stops)
        slope, intercept, sse = linfit(view(r,a:b), view(z,a:b))
        ŷ = intercept .+ slope .* r[a:b]
        sst = sum((z[a:b] .- mean(z[a:b])).^2)
        r2 = sst ≈ 0 ? 1.0 : 1 - sse/sst
        push!(alphas, -slope) # alpha = -pente
        push!(info, (slope=slope, intercept=intercept, sse=sse, r2=r2, range=(a,b)))
    end

    return (alphas, starts, info)
end
