# --------------------------------------------------------------------
#                         PROBLEM STRUCTURE
# --------------------------------------------------------------------
"""
    AtomProblem

Defines an atomic simulation problem for SCF (self-consistent field) resolution.

This structure gathers all the physical and numerical parameters needed to solve
the ground-state electronic structure of an atom or ion.

# Fields
- `z::Real`: Nuclear charge.
- `N::Real`: Number of electrons.
- `hartree::Real` : Prefactor in front of the Hartree term.
- `ex` : Exchange functional.
- `ec` : Correlation functional.
- `lh::Int`: Maximum orbital angular momentum quantum number `l`.
- `nh::Int`: Maximum principal quantum number `n`.
- `Nmesh::Int`: Number of discretization points.
- `Rmax::Real`: Radial domain cut-off.
- `typemesh`: Mesh constructor (e.g. `ExponentialMesh`).
- `optsmesh::NamedTuple`: Mesh options.
- `typebasis`: Basis constructor (e.g. `P1IntLegendreBasis`).
- `optsbasis::NamedTuple`: Basis options.
- `integration_method::FEMIntegrationMethod`: Integration method for radial quadrature.
- `optsintegration::NamedTuple`: FEL Integration method options.
- `alg::SCFAlgorithm`: SCF algorithm (e.g. `ODA`, `Quadratic`).
- `name::AbstractString`: Optional problem name.
- `solveropts`: Additional solver options.
"""
struct AtomProblem{T <: Real,
                   EX,
                   EC,
                   TM,
                   OM <: NamedTuple,
                   TB,
                   OB <: NamedTuple,
                   FM,
                   OFM <: NamedTuple,
                   A <: SCFAlgorithm,
                   S <: AbstractString,
                   OS }
    # Model parameters
    z::T
    N::T
    hartree::T
    ex::EX
    ec::EC
    # Discretization parameters
    lh::Int
    nh::Int
    Nmesh::Int
    Rmax::T
    typemesh::TM
    optsmesh::OM
    typebasis::TB
    optsbasis::OB
    integration_method::FM
    optsintegration::OFM
    # Algorithm
    alg::A
    # Metadata
    name::S
    solveropts::OS
end
# --------------------------------------------------------------------
#                        CONSTRUCTORS
# --------------------------------------------------------------------

"""
    AtomProblem(; T, z, N, hartree, ex, ec,
                  lh, nh, Nmesh, Rmax,
                  typemesh, optsmesh,
                  typebasis, optsbasis,
                  integration_method, optsintegration,
                  alg;
                  name = "", kwargs...)

Constructs an `AtomProblem` from keyword arguments.

# Keyword Arguments

## Physical model
- `z::Real`: Nuclear charge.
- `N::Real`: Number of electrons.
- `hartree::Real`: Hartree prefactor.
- `ex`: Exchange functional.
- `ec`: Correlation functional.

## Discretization
- `T::DataType`: Number type used for numerical computations (e.g., `Float64`, `Double64`).
- `lh::Int`: Angular quantum number truncations.
- `nh::Int`: Principal quantum number truncations.
- `Nmesh::Int`: Number of discretization points.
- `Rmax::Real`: Radial cut-off.
- `typemesh`: Mesh constructor.
- `optsmesh::NamedTuple`: Mesh options.
- `typebasis`: Basis constructor.
- `optsbasis::NamedTuple`: Basis options.
- `integration_method`: Quadrature method.
- `optsintegration::NamedTuple`: FEL Integration method options.

## Algorithm
- `alg::SCFAlgorithm`: SCF algorithm.

## Optional
- `name::AbstractString`: Name of the problem.
- `kwargs`: Additional solver options.
"""
function AtomProblem(; T::DataType, z::Real, N::Real, hartree::Real, ex::EX, ec::EC,
                      lh::Int, nh::Int, Nmesh::Int, Rmax::Real,
                      typemesh, optsmesh::NamedTuple,
                      typebasis, optsbasis::NamedTuple,
                      integration_method, optsintegration::NamedTuple,
                      alg::SCFAlgorithm,
                      name::AbstractString = "", kwargs...) where {EX, EC}
    _name = if name == "" && floor(z) == z
        charge_diff = z - N
        if iszero(charge_diff)
            ATOMIC_NUMBER_TO_NAME[Int(z)]
        else
            symbol_ion = charge_diff > 0 ? "+" : "-"
            ATOMIC_NUMBER_TO_NAME[Int(z)]*string(abs(charge_diff))*symbol_ion
        end
    else
        name
    end
    return AtomProblem{T, typeof(ex), typeof(ec),
                       typeof(typemesh), typeof(optsmesh),
                       typeof(typebasis), typeof(optsbasis),
                       typeof(integration_method), typeof(optsintegration),
                       typeof(alg), typeof(_name), typeof(kwargs)}(
        T(z), T(N), T(hartree), ex, ec,
        lh, nh, Nmesh, T(Rmax),
        typemesh, optsmesh,
        typebasis, optsbasis,
        integration_method, optsintegration,
        alg, _name, kwargs
    )
end

"""
    AtomProblem(prob::AtomProblem; overrides...)

Copy constructor for `AtomProblem`, allowing field overrides.

Useful for creating a modified version of an existing problem.

All fields default to `prob.<field>`, and can be selectively overridden.
"""
function AtomProblem(prob::AtomProblem;
                     z::Real = prob.z,
                     N::Real = prob.N,
                     hartree::Real = prob.hartree,
                     ex::Ex = prob.ex,
                     ec::Ec = prob.ec,
                     lh::Int = prob.lh,
                     nh::Int = prob.nh,
                     Nmesh::Int = prob.Nmesh,
                     Rmax::Real = prob.Rmax,
                     typemesh = prob.typemesh,
                     optsmesh::NamedTuple = prob.optsmesh,
                     typebasis = prob.typebasis,
                     optsbasis::NamedTuple = prob.optsbasis,
                     integration_method = prob.integration_method,
                     optsintegration::NamedTuple = prob.optsintegration,
                     alg::SCFAlgorithm = prob.alg,
                     name::AbstractString = prob.name,
                     kwargs = prob.solveropts) where {Ex, Ec}
    T = _datatype(prob)
    return AtomProblem{T, typeof(ex), typeof(ec),
                       typeof(typemesh), typeof(optsmesh),
                       typeof(typebasis), typeof(optsbasis),
                       typeof(integration_method), typeof(optsintegration),
                       typeof(alg), typeof(name), typeof(kwargs)}(
        T, z, N, hartree, ex, ec,
        lh, nh, Nmesh, Rmax,
        typemesh, optsmesh,
        typebasis, optsbasis,
        integration_method, optsintegration,
        alg, name, kwargs
    )
end


_datatype(::AtomProblem{T}) where {T} = T


"""
    ATOMIC_NUMBER_TO_NAME

Dictionary mapping atomic numbers (Z) to element names in English.

Useful for labeling atoms in periodic table computations or visualizations.
"""
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


"""
    ATOMIC_NUMBER_TO_SYMBOL

Dictionary mapping atomic numbers (Z) to atomic symbols.

Useful for displaying or labeling atoms in scientific applications.
"""
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
