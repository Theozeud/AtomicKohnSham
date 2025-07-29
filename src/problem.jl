# --------------------------------------------------------------------
#                         PROBLEM STRUCTURE
# --------------------------------------------------------------------
"""
    AtomProblem

Defines an atomic simulation problem for SCF (self-consistent field) resolution.

This structure gathers all the physical and numerical parameters needed to solve
the ground-state electronic structure of an atom or ion.

# Fields
- `lh::Int`: Maximum orbital angular momentum quantum number `l`.
- `nh::Int`: Maximum principal quantum number `n`.
- `alg::SCFAlgorithm`: SCF algorithm (e.g. `ODA`, `Quadratic`).
- `z::TZ`: Nuclear charge.
- `N::TN`: Number of electrons.
- `hartree`: Prefactor in front of the Hartree term.
- `integration_method::FM`: Integration method for radial quadrature.
- `Rmax::T`: Radial domain cut-off.
- `Nmesh::Int`: Number of discretization points.
- `typemesh`: Mesh constructor (e.g. `ExponentialMesh`).
- `optsmesh::NamedTuple`: Mesh options.
- `typebasis`: Basis constructor (e.g. `IntLegendre`).
- `optsbasis::NamedTuple`: Basis options.
- `name::AbstractString`: Optional problem name.
- `solveropts`: Additional solver options.
"""
struct AtomProblem{T <: Real,
    A <: SCFAlgorithm,
    M <: KSEModel,
    S <: AbstractString,
    TZ, TN, FM,
    OM <: NamedTuple,
    OB <: NamedTuple,
    OS,
    TM,
    TB}
    lh::Int
    nh::Int
    alg::A
    z::TZ
    N::TN
    hartree::Any
    integration_method::FM
    Rmax::T
    Nmesh::Int
    typemesh::TM
    optsmesh::OM
    typebasis::TB
    optsbasis::OB
    name::S
    solveropts::OS
end

"""
    AtomProblem(; T, lh, nh, alg, z, N, hartree, integration_method,
                  Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis,
                  name = "", kwargs...)

Constructor for `AtomProblem` from keyword arguments.

# Keyword

- `T`: Number type (e.g. `Float64`, `Double64`).
- `lh`, `nh`: Cutoffs for angular and principal quantum numbers.
- `alg`: SCF algorithm.
- `z`: Nuclear charge.
- `N`: Number of electrons.
- `hartree`: Coefficient in front of Hartree term.
- `integration_method`: Quadrature integration method.
- `Rmax`, `Nmesh`: Domain size and discretization.
- `typemesh`, `optsmesh`: Mesh constructor and options.
- `typebasis`, `optsbasis`: Basis constructor and options.
- `name`: Optional string name.
- `kwargs`: Additional solver options.
"""
function AtomProblem(;
        T,
        lh,
        nh,
        alg,
        z,
        N,
        hartree,
        integration_method,
        Rmax,
        Nmesh,
        typemesh,
        optsmesh,
        typebasis, optsbasis,
        name = "", kwargs...)
    return AtomProblem{
        T,
        typeof(alg),
        typeof(z),
        typeof(N),
        typeof(hartree),
        typeof(integration_method),
        typeof(optsmesh),
        typeof(optsbasis),
        typeof(kwargs),
        typeof(typemesh),
        typeof(typebasis),
        typeof(name)}(
        lh,
        nh,
        alg,
        z,
        N,
        hartree,
        integration_method,
        Rmax,
        Nmesh,
        typemesh,
        optsmesh,
        typebasis,
        optsbasis,
        name,
        kwargs
    )
end

"""
    AtomProblem(prob::AtomProblem; kwargs...)

Copy constructor for `AtomProblem` with optional overrides.

Useful for creating a new problem from an existing one while updating fields.
"""
function AtomProblem(prob::AtomProblem;
        T = typeof(prob.Rmax),
        lh = prob.lh, nh = prob.nh,
        alg = prob.alg, z = prob.z, N = prob.N,
        hartree = prob.hartree,
        integration_method = prob.integration_method,
        Rmax = prob.Rmax, Nmesh = prob.Nmesh,
        typemesh = prob.typemesh, optsmesh = prob.optsmesh,
        typebasis = prob.typebasis, optsbasis = prob.optsbasis,
        name = prob.name, kwargs = prob.solveropts)
    return AtomProblem{T, typeof(alg), typeof(z), typeof(N), typeof(hartree),
        typeof(integration_method), typeof(optsmesh),
        typeof(optsbasis), typeof(kwargs),
        typeof(typemesh), typeof(typebasis), typeof(name)}(
        lh, nh, alg, z, N, hartree, integration_method,
        Rmax, Nmesh, typemesh, optsmesh,
        typebasis, optsbasis, name, kwargs
    )
end

datatype(::AtomProblem{T}) where {T} = T

const ATOMIC_NUMBER_TO_NAME = Dict(
    1 => "Hydrogen",
    2 => "Helium",
    3 => "Lithium",
    4 => "Beryllium",
    5 => "Boron",
    6 => "Carbon",
    7 => "Nitrogen",
    8 => "Oxygen",
    9 => "Fluorine",
    10 => "Neon",
    11 => "Sodium",
    12 => "Magnesium",
    13 => "Aluminium",
    14 => "Silicon",
    15 => "Phosphorus",
    16 => "Sulfur",
    17 => "Chlorine",
    18 => "Argon",
    19 => "Potassium",
    20 => "Calcium",
    21 => "Scandium",
    22 => "Titanium",
    23 => "Vanadium",
    24 => "Chromium",
    25 => "Manganese",
    26 => "Iron",
    27 => "Cobalt",
    28 => "Nickel",
    29 => "Copper",
    30 => "Zinc",
    31 => "Gallium",
    32 => "Germanium",
    33 => "Arsenic",
    34 => "Selenium",
    35 => "Bromine",
    36 => "Krypton",
    37 => "Rubidium",
    38 => "Strontium",
    39 => "Yttrium",
    40 => "Zirconium",
    41 => "Niobium",
    42 => "Molybdenum",
    43 => "Technetium",
    44 => "Ruthenium",
    45 => "Rhodium",
    46 => "Palladium",
    47 => "Silver",
    48 => "Cadmium",
    49 => "Indium",
    50 => "Tin",
    51 => "Antimony",
    52 => "Tellurium",
    53 => "Iodine",
    54 => "Xenon",
    55 => "Cesium",
    56 => "Barium",
    57 => "Lanthanum",
    58 => "Cerium",
    59 => "Praseodymium",
    60 => "Neodymium",
    61 => "Promethium",
    62 => "Samarium",
    63 => "Europium",
    64 => "Gadolinium",
    65 => "Terbium",
    66 => "Dysprosium",
    67 => "Holmium",
    68 => "Erbium",
    69 => "Thulium",
    70 => "Ytterbium",
    71 => "Lutetium",
    72 => "Hafnium",
    73 => "Tantalum",
    74 => "Tungsten",
    75 => "Rhenium",
    76 => "Osmium",
    77 => "Iridium",
    78 => "Platinum",
    79 => "Gold",
    80 => "Mercury",
    81 => "Thallium",
    82 => "Lead",
    83 => "Bismuth",
    84 => "Polonium",
    85 => "Astatine",
    86 => "Radon",
    87 => "Francium",
    88 => "Radium",
    89 => "Actinium",
    90 => "Thorium",
    91 => "Protactinium",
    92 => "Uranium",
    93 => "Neptunium",
    94 => "Plutonium",
    95 => "Americium",
    96 => "Curium",
    97 => "Berkelium",
    98 => "Californium",
    99 => "Einsteinium",
    100 => "Fermium",
    101 => "Mendelevium",
    102 => "Nobelium",
    103 => "Lawrencium",
    104 => "Rutherfordium",
    105 => "Dubnium",
    106 => "Seaborgium",
    107 => "Bohrium",
    108 => "Hassium",
    109 => "Meitnerium",
    110 => "Darmstadtium",
    111 => "Roentgenium",
    112 => "Copernicium",
    113 => "Nihonium",
    114 => "Flerovium",
    115 => "Moscovium",
    116 => "Livermorium",
    117 => "Tennessine",
    118 => "Oganesson"
)

const ATOMIC_NUMBER_TO_SYMBOL = Dict(
    1 => "H",
    2 => "He",
    3 => "Li",
    4 => "Be",
    5 => "B",
    6 => "C",
    7 => "N",
    8 => "O",
    9 => "F",
    10 => "Ne",
    11 => "Na",
    12 => "Mg",
    13 => "Al",
    14 => "Si",
    15 => "P",
    16 => "S",
    17 => "Cl",
    18 => "Ar",
    19 => "K",
    20 => "Ca",
    21 => "Sc",
    22 => "Ti",
    23 => "V",
    24 => "Cr",
    25 => "Mn",
    26 => "Fe",
    27 => "Co",
    28 => "Ni",
    29 => "Cu",
    30 => "Zn",
    31 => "Ga",
    32 => "Ge",
    33 => "As",
    34 => "Se",
    35 => "Br",
    36 => "Kr",
    37 => "Rb",
    38 => "Sr",
    39 => "Y",
    40 => "Zr",
    41 => "Nb",
    42 => "Mo",
    43 => "Tc",
    44 => "Ru",
    45 => "Rh",
    46 => "Pd",
    47 => "Ag",
    48 => "Cd",
    49 => "In",
    50 => "Sn",
    51 => "Sb",
    52 => "Te",
    53 => "I",
    54 => "Xe",
    55 => "Cs",
    56 => "Ba",
    57 => "La",
    58 => "Ce",
    59 => "Pr",
    60 => "Nd",
    61 => "Pm",
    62 => "Sm",
    63 => "Eu",
    64 => "Gd",
    65 => "Tb",
    66 => "Dy",
    67 => "Ho",
    68 => "Er",
    69 => "Tm",
    70 => "Yb",
    71 => "Lu",
    72 => "Hf",
    73 => "Ta",
    74 => "W",
    75 => "Re",
    76 => "Os",
    77 => "Ir",
    78 => "Pt",
    79 => "Au",
    80 => "Hg",
    81 => "Tl",
    82 => "Pb",
    83 => "Bi",
    84 => "Po",
    85 => "At",
    86 => "Rn",
    87 => "Fr",
    88 => "Ra",
    89 => "Ac",
    90 => "Th",
    91 => "Pa",
    92 => "U",
    93 => "Np",
    94 => "Pu",
    95 => "Am",
    96 => "Cm",
    97 => "Bk",
    98 => "Cf",
    99 => "Es",
    100 => "Fm",
    101 => "Md",
    102 => "No",
    103 => "Lr",
    104 => "Rf",
    105 => "Db",
    106 => "Sg",
    107 => "Bh",
    108 => "Hs",
    109 => "Mt",
    110 => "Ds",
    111 => "Rg",
    112 => "Cn",
    113 => "Nh",
    114 => "Fl",
    115 => "Mc",
    116 => "Lv",
    117 => "Ts",
    118 => "Og"
)
