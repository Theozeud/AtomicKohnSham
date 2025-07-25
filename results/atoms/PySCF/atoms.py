from pyscf import gto, scf, ci, tools
import numpy as np

# Définir la molécule LiCN
mol = gto.Mole()
mol.build(
    atom = '''                      
    Li 0.000000 0.000000 -3.6830000
    C  0.000000 0.000000  0.000000 
    N  0.000000 0.000000  2.168
    ''',
    basis = 'ccpvdz',
)

mol.build(symmetry=True, symmetry_subgroup = 'C2v')

# Calcul Hartree-Fock
mf = scf.RHF(mol)
mf.run(max_cycle=100)

# CISD
myci = ci.CISD(mf)
myci.conv_tol = 1e-8
myci.run()

# Énergie totale CISD
print(f"Énergie totale CISD : {myci.e_tot:.6f} Hartree")

# Moment dipolaire (approximé avec les orbitales HF)
dipole_moment = mf.dip_moment()
print(f"Moment dipolaire (approx. HF) : {dipole_moment} Debye")

# Obtenir la matrice densité à une particule (RHF : 1 matrice)
rdm1 = myci.make_rdm1()

# Diagonalisation pour obtenir les orbitales naturelles
# Cela donne les orbitales naturelles (colonnes de `natorbs`) et leurs occupations
occ, natorbs = np.linalg.eigh(rdm1)

# Ordonner par occupation décroissante
idx = np.argsort(-occ)
occ = occ[idx]
natorbs = natorbs[:, idx]

# Afficher les occupations naturelles
print("Occupations des orbitales naturelles (CISD) :")
for i, o in enumerate(occ):
    print(f"  Orb {i+1:2d}: {o:.6f}")

# Pour les orbitales naturelles (après diagonalisation de la densité) :
with open("natorb.molden", "w") as f:
    tools.molden.header(mol, f)
    tools.molden.orbital_coeff(mol, f, natorbs)


print(f"Énergie de répulsion nucléaire : {mol.energy_nuc():.6f} Hartree")


print("\nRésumé des propriétés de la molécule :")
print(f"  Nombre d'atomes       : {mol.natm}")
print(f"  Charge totale         : {mol.charge}")
print(f"  Multiplicité de spin  : {mol.spin + 1}")
print(f"  Nombre d'électrons    : {mol.nelectron}")
print(f"  Nombre d'orbitales MO : {mf.mo_coeff.shape[1]}")
print(f"  Énergie HF totale     : {mf.e_tot:.6f} Hartree")

print(f"Coordinates : {mol.atom_coord()}")