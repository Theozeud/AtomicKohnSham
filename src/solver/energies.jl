"""
    Energies{T<:Real}

Container for the energy components computed during the SCF iterations.

Fields (all of type `T`):
- `Etot`    : total energy
- `Ekin`    : kinetic energy
- `Ecou`    : Coulomb energy
- `Ehar`    : Hartree energy
- `Eexc`    : exchange–correlation energy
- `Ekincor` : kinetic correlation energy (when applicable)
"""
mutable struct Energies{T<:Real}
    Etot::T
    Ekin::T
    Ecou::T
    Ehar::T
    Eexc::T
    Ekincor::T
    function Energies(::Type{T}) where {T<:Real}
        z = zero(T)
        new{T}(z, z, z, z, z, z)
    end
end
