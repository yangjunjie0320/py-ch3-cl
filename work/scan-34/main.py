import sys
import pyscf
from pyscf import tdscf

def main(xyz_file, basis="631g*", xc="b3lyp", output="ene.out"):
    xyz = None

    with open(xyz_file, "r") as xyz:
        xyz = xyz.readlines()[2:]
        xyz = "".join(xyz)

    mol = pyscf.gto.Mole()
    mol.atom = xyz
    mol.basis = basis
    mol.verbose = 5
    mol.build()

    if xc is not None and xc != "hf":
        mf = pyscf.scf.RKS(mol)
        mf.xc = xc
        mf.grids.level = 9
        mf.kernel()
    else:
        mf = pyscf.scf.RHF(mol)
        mf.kernel()

    ene_mf = mf.e_tot

    with open(output, "w") as f:
        tda = pyscf.tdscf.TDA(mf)
        tda.nstates = 5
        tda.max_space = 8
        tda.max_cycle = 2000
        tda.kernel()

        con_tda = tda.converged
        ene_tda = tda.e[:5]
        dip_tda = tda.transition_dipole()[:5]

        f.write(f"RKS energy: {ene_mf: 12.8f}\n")
        f.write("TDA results:\n")
        for i, (con, ene, dip) in enumerate(zip(con_tda, ene_tda, dip_tda)):
            f.write(f" {i:2d} {con} {ene: 12.8f} {dip[0]: 12.8f} {dip[1]: 12.8f} {dip[2]: 12.8f}\n")

        # rpa = pyscf.tdscf.RPA(mf)
        # rpa.nstates = 5
        # rpa.max_space = 8
        # rpa.max_cycle = 2000
        # rpa.kernel()

        # con_rpa = rpa.converged
        # ene_rpa = rpa.e[:5]
        # dip_rpa = rpa.transition_dipole()[:5]
        
        # f.write("\nRPA results:\n")
        # for i, (con, ene, dip) in enumerate(zip(con_rpa, ene_rpa, dip_rpa)):
        #     f.write(f" {i:2d} {con} {ene: 12.8f} {dip[0]: 12.8f} {dip[1]: 12.8f} {dip[2]: 12.8f}\n")

if __name__ == "__main__":
    xyz = sys.argv[1]
    main(xyz, basis="6311g**", xc="b3lyp", output="out.log")
