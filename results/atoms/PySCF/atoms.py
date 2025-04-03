from pyscf import gto, scf, ci

# Définir la molécule LiCN
mol = gto.Mole()
mol.build(
    # Position des atomes (cartésienne)
    atom = '''                      
    Li 0.000000 0.000000 0.000000
    C  1.100000 0.000000 0.000000
    N  2.200000 0.000000 0.000000
    ''',
    # Base utilisée pour les calculs                        
    basis = '6-31G',
    # Pas de symétrie  
    symmetry = False  
)

# Effectuer un calcul Hartree-Fock
mf = scf.RHF(mol)
mf.kernel()

# Créer un objet CI avec le calcul Hartree-Fock
myci = ci.CISD(mf)

# Effectuer le calcul CI
myci.kernel()

# Afficher les résultats
print(f"Énergie totale CI (CISD) : {myci.e_tot:.6f} Hartree")

# Calculer le moment dipolaire
dipole_moment = mf.dip_moment()

# Afficher le moment dipolaire
print(f"Le moment dipolaire de la molécule LiCN est : {dipole_moment}")