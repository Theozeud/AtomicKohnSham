#--------------------------------------------------------------------
#                               LogConfig
#--------------------------------------------------------------------
struct LogConfig
    occupation_number::Bool
    orbitals_energy::Bool
    stopping_criteria::Bool
    energy::Bool
    density::Bool

    function LogConfig(;occupation_number = false, orbitals_energy = false, stopping_criteria = true, energy = false, density = false)
        new(occupation_number, orbitals_energy, stopping_criteria, energy, density)
    end
end


function Base.show(io::IO, lc::LogConfig)
    println(io, "LogConfig:")
    println(io, "  occupation_number  = ", lc.occupation_number)
    println(io, "  orbitals_energy    = ", lc.orbitals_energy)
    println(io, "  stopping_criteria  = ", lc.stopping_criteria)
    println(io, "  energy             = ", lc.energy)
    println(io, "  density            = ", lc.density)
end


#--------------------------------------------------------------------
#                               LogBook
#--------------------------------------------------------------------
struct LogBook
    config::LogConfig
    occupation_number_log
    orbitals_energy_log
    stopping_criteria_log
    energy_log
    density_log

    function LogBook(config, T)
        occupation_number_log = []                  # To see
        orbitals_energy_log   = []                  # To see
        stopping_criteria_log = T[]                 # To see
        energy_log            = T[]                 # To see
        density_log           = []                  # To see
        new(config, occupation_number_log, orbitals_energy_log, stopping_criteria_log, energy_log, density_log)
    end
end