using CairoMakie
using Random

CairoMakie.activate!()

# Symboles des 118 éléments
symbols = [
    "H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe",
    "Cs","Ba",
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra",
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",
    "Nh","Fl","Mc","Lv","Ts","Og"
]

# Numéros atomiques
Z = 1:118

# Nombres aléatoires associés
ioni =[
    0.70, 0.38, 0.67, 0.37, 0.55, 0.67, 0.63, 0.78, 0.87, 0.35,
    0.68, 0.35, 0.57, 0.70, 0.70, 0.82, 0.98, 0.37, 0.70, 0.42,
    0.46, 0.46, 0.72, 0.59, 0.55, 0.74, 0.82, 0.61, 0.80, 0.31,
    0.52, 0.70, 0.67, 0.80, 1.00, 0.37, 0.70, 0.48, 0.46, 0.65,
    0.55, 0.82, 0.80, 0.72, 0.61, 0.63, 0.78, 0.35, 0.57, 0.76,
    0.72, 0.85, 1.00, 0.38, 0.70, 0.52, 0.61, 0.53, 0.46, 0.40,
    0.61, 0.44, 0.46, 0.50, 0.35, 0.35, 0.33, 1.06, 0.10, 0.44,
    0.48, 0.53, 0.57, 0.67, 0.83, 0.74, 0.59, 0.63, 0.82, 0.37,
    0.57, 0.67, 0.74, 0.80, 1.00, 0.38, 0.72, 0.55, 0.59, 0.57,
    0.61, 0.61, 0.74, 0.74, 0.70, 0.65, 1.19, 1.02, 0.59, 0.23,
    0.70, 0.40, 0.57, 0.67, 0.72, 0.76, 0.59, 0.50, 0.65, 0.61,
    0.80, 0.38, 0.59, 0.74, 0.87, 0.97, 1.00, 0.38
]
using Printf
randvals = [@sprintf("%.2f", x) for x in ioni]

# Positionnement du tableau (col, row) de chaque élément
# (simplifié : lanthanides et actinides placés sous forme de 2 rangées séparées)
positions = Dict{Int,Tuple{Int,Int}}()

# --- Période 1
positions[1] = (1,1)   # H
positions[2] = (18,1)  # He

# --- Période 2
positions[3] = (1,2);  positions[4] = (2,2)   # Li, Be
for (i,z) in enumerate(5:10)                  # B → Ne
    positions[z] = (i+12,2)
end

# --- Période 3
positions[11] = (1,3); positions[12] = (2,3)  # Na, Mg
for (i,z) in enumerate(13:18)                 # Al → Ar
    positions[z] = (i+12,3)
end

# --- Période 4
for (i,z) in enumerate(19:36)
    positions[z] = (i,4)
end

# --- Période 5
for (i,z) in enumerate(37:54)
    positions[z] = (i,5)
end

# --- Période 6 (Cs, Ba, Hf → Rn)
positions[55] = (1,6)   # Cs
positions[56] = (2,6)   # Ba
for (i,z) in enumerate(72:86)   # Hf → Rn
    positions[z] = (i+3,6)      # décalé de 3 colonnes (col=4 → Hf)
end
# Lanthanides (57–71), affichés en ligne 9
for (i,z) in enumerate(57:71)
    positions[z] = (i+3,9)      # placés sous le trou (col=4 → La)
end

# --- Période 7 (Fr, Ra, Rf → Og)
positions[87] = (1,7)   # Fr
positions[88] = (2,7)   # Ra
for (i,z) in enumerate(104:118) # Rf → Og
    positions[z] = (i+3,7)      # décalé de 3 colonnes (col=4 → Rf)
end
# Actinides (89–103), affichés en ligne 10
for (i,z) in enumerate(89:103)
    positions[z] = (i+3,10)     # col=4 → Ac
end


### COULEURS

# Familles par numéro atomique
families = Dict(
    # Groupe 1 : alcalins
    1=>:nonmetal, 3=>:alkali, 11=>:alkali, 19=>:alkali, 37=>:alkali, 55=>:alkali, 87=>:alkali,

    # Groupe 2 : alcalino-terreux
    4=>:alkaline, 12=>:alkaline, 20=>:alkaline, 38=>:alkaline, 56=>:alkaline, 88=>:alkaline,

    # Lanthanides
    57=>:lanthanide, 58=>:lanthanide, 59=>:lanthanide, 60=>:lanthanide, 61=>:lanthanide,
    62=>:lanthanide, 63=>:lanthanide, 64=>:lanthanide, 65=>:lanthanide, 66=>:lanthanide,
    67=>:lanthanide, 68=>:lanthanide, 69=>:lanthanide, 70=>:lanthanide, 71=>:lanthanide,

    # Actinides
    89=>:actinide, 90=>:actinide, 91=>:actinide, 92=>:actinide, 93=>:actinide,
    94=>:actinide, 95=>:actinide, 96=>:actinide, 97=>:actinide, 98=>:actinide,
    99=>:actinide, 100=>:actinide, 101=>:actinide, 102=>:actinide, 103=>:actinide,

    # Groupe 17 : halogènes
    9=>:halogen, 17=>:halogen, 35=>:halogen, 53=>:halogen, 85=>:halogen, 117=>:halogen,

    # Groupe 18 : gaz nobles
    2=>:noble, 10=>:noble, 18=>:noble, 36=>:noble, 54=>:noble, 86=>:noble, 118=>:noble,

    # Métalloïdes
    5=>:metalloid, 14=>:metalloid, 32=>:metalloid, 33=>:metalloid, 51=>:metalloid, 52=>:metalloid, 84=>:metalloid,

    # Métaux pauvres (post-transition)
    13=>:post, 31=>:post, 49=>:post, 50=>:post, 81=>:post, 82=>:post, 83=>:post, 113=>:post, 114=>:post, 115=>:post, 116=>:post,

    # Non-métaux
    6=>:nonmetal, 7=>:nonmetal, 8=>:nonmetal, 15=>:nonmetal, 16=>:nonmetal,

    # Hydrogène déjà mis en nonmetal plus haut

    # Métaux de transition (les autres qui ne sont pas explicitement classés)
)
# Attribution générique des métaux de transition
for z in [21:30; 39:48; 72:80; 104:112]
    families[z] = :transition
end

# --- Couleurs par famille ---
family_colors = Dict(
    :nonmetal    => :lightgreen,
    :alkali      => :orange,
    :alkaline    => :yellow,
    :transition  => :lightskyblue,
    :post        => :lightgray,
    :metalloid   => :khaki,
    :halogen     => :lightseagreen,
    :noble       => :violet,
    :lanthanide  => :pink,
    :actinide    => :red
)

# Définir familles aux couleurs "claires"
light_families = Set([:alkaline, :metalloid, :post, :alkali, :lanthanide])

function text_color_for(fam)
    if fam in light_families
        return :black
    else
        return :white
    end
end




# --- Dessin ---
f = Figure(resolution=(2800,1600))
ax = Axis(f[1,1], xticksvisible=false, yticksvisible=false, xgridvisible=false, ygridvisible=false,
          xlabel="Tableau périodique", ylabel="", aspect=DataAspect())

# Désactiver axes
hidedecorations!(ax)

for z in 1:118
    (col,row) = positions[z]
    # Dessiner la case
    #poly!(ax, Rect(col, -row-1, 1, 1), color=:white, strokecolor=:black)
    fam = get(families, z, :post)              # famille de l’élément z
    color = family_colors[fam]                    # couleur associée
    #poly!(ax, Rect(col_, -row-1, 1, 1), color=col, strokecolor=:black)

    poly!(ax, Rect(col, -row-1, 1, 1), color=color, strokecolor=:black, strokewidth=1)


    # Texte : symbole en haut gauche
    text!(ax, symbols[z], position=(col+0.05, -row-0.15), align=(:left,:top), fontsize=30, color=:blue)

    # Texte : Z en haut droit
    text!(ax, "z=$(z)", position=(col+0.95, -row-0.15), align=(:right,:top), fontsize=30, color=:black)

    # Nombre aléatoire au centre
    #text!(ax, string(randvals[z]), position=(col+0.5, -row-0.55), align=(:center,:center), fontsize=20, color=:white)

    # Texte du centre (nombre aléatoire) avec couleur adaptée
    txtcol = text_color_for(fam)
    text!(ax, string(randvals[z]),
        position=(col+0.5, -row-0.7),
        align=(:center,:center), fontsize=50, color=txtcol)
end

# Case de légende (taille 2×2, placée en bas à gauche)
poly!(ax, Rect(1, -11, 2, 2), color=:white, strokecolor=:black, strokewidth=2)

# Exemple : symbole
text!(ax, "Symbol", position=(1.2, -9.1), align=(:left,:top), fontsize=34, color=:blue)

# Exemple : numéro atomique
text!(ax, "z", position=(2.9, -9.1), align=(:right,:top), fontsize=34, color=:black)

# Exemple : nombre au centre
text!(ax, " Maximal \n Ionization", position=(2, -10.0), align=(:center,:center), fontsize=45, color=:black)



display(f)
save("tableau_periodique.png", f)
