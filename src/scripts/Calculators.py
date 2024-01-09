import paths
import numpy as np
from ase.calculators.eam import EAM

# TODO: For M3GNet and CHGNet is stress_weight=1.0*GPa or stress_weight=1.0?
from ase.units import GPa


# No longer needed using conda LAMMPS, which provides the needed stress.
class CustomEAM(EAM):
    """
    Override EAM to zero out stress for ASE Database support
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_stress(self, atoms=None, **kwargs):
        return np.zeros(6)


def get_ase_calculator(model="MutterASE"):
    """
    Create and return an ASE calculator object.

    Args:
        model (string): The name of the interatomic potential model

    Returns:
        ASE.calculators

    Notes:
        For the EAM potentials there is the option to select the
        internal ASE calculator for EAM. This works but does not
        provide the stress there for we return a customEAM object
        which returns a zero stress array.

        EAM/MEAM potentials are referenced via first author last
        name, this means no specification of chemical system is
        required.

    References:
        - Mutter, D.; Nielaba, Phys. Rev. B 2010, 82 (22), 224201. https://doi.org/10.1103/PhysRevB.82.224201.
        - Zhong, Y.; Gall, K.; Zhu, T. Atomistic Study of Nanotwins in NiTi Shape Memory Alloys. Journal of Applied Physics 2011, 110 (3), 033532. https://doi.org/10.1063/1.3621429.
        - Kim, J.-S.; et al.  Calphad 2017, 59, 131–141. https://doi.org/10.1016/j.calphad.2017.09.005.
        - Ko, W.-S.; Grabowski, B.; Neugebauer, J. Phys. Rev. B 2015, 92 (13), 134107. https://doi.org/10.1103/PhysRevB.92.134107.
        - Kavousi, S.; et al. Modelling Simul. Mater. Sci. Eng. 2019, 28 (1), 015006. https://doi.org/10.1088/1361-651X/ab580c.
        - Chen, C.; Ong, S. P. Nat Comput Sci 2022, 2 (11), 718–728. https://doi.org/10.1038/s43588-022-00349-3.
        - Deng, B.; et al. Nat Mach Intell 2023, 5 (9), 1031–1041. https://doi.org/10.1038/s42256-023-00716-3.
        - Batatia, I.; et al. MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields; 2022.
        - Tang, H.; et al. Acta Materialia 2022, 238, 118217. https://doi.org/10.1016/j.actamat.2022.118217.

    """
    if model == "Mutter":
        from ase.calculators.lammpslib import LAMMPSlib

        potential_file = str(paths.static / "NiTi_Mutter.eam.fs")
        cmds = ["pair_style eam/fs", f"pair_coeff * *  {potential_file} Ni Ti"]
        amds = ["thermo_style custom etotal lx ly lz vol pxx pyy pzz pxy pxz pyz press"]
        asecalc = LAMMPSlib(lmpcmds=cmds, amendments=amds, keep_alive=True)

    elif model == "MutterASE":
        potential_file = str(paths.static / "NiTi_Mutter.eam.fs")
        asecalc = CustomEAM(potential=potential_file)

    elif model == "Zhong":
        from ase.calculators.lammpslib import LAMMPSlib

        potential_file = str(paths.static / "NiTi_Zhong.eam.fs")
        cmds = ["pair_style eam/fs", f"pair_coeff * *  {potential_file} Ni Ti"]
        amds = ["thermo_style custom etotal lx ly lz vol pxx pyy pzz pxy pxz pyz press"]
        asecalc = LAMMPSlib(lmpcmds=cmds, amendments=amds, keep_alive=True)

    elif model == "ZhongASE":
        potential_file = str(paths.static / "NiTi_Zhong.eam.fs")
        asecalc = CustomEAM(potential=potential_file)

    elif model == "Ko":
        from ase.calculators.lammpslib import LAMMPSlib

        library_file = str(paths.static / "NiTi_Ko.meam.library")
        potential_file = str(paths.static / "NiTi_Ko.meam.potential")
        cmds = [
            "pair_style meam",
            f"pair_coeff * *  {library_file} Ni Ti {potential_file} Ni Ti",
        ]
        amds = ["thermo_style custom etotal lx ly lz vol pxx pyy pzz pxy pxz pyz press"]
        asecalc = LAMMPSlib(lmpcmds=cmds, amendments=amds, keep_alive=True)

    elif model == "Kavousi":
        from ase.calculators.lammpslib import LAMMPSlib

        library_file = str(paths.static / "NiTi_Kavousi.meam.library")
        potential_file = str(paths.static / "NiTi_Kavousi.meam.potential")
        cmds = [
            "pair_style meam",
            f"pair_coeff * *  {library_file} Ni Ti {potential_file} Ni Ti",
        ]
        amds = ["thermo_style custom etotal lx ly lz vol pxx pyy pzz pxy pxz pyz press"]
        asecalc = LAMMPSlib(lmpcmds=cmds, amendments=amds, keep_alive=True)

    elif model == "Kim":
        from ase.calculators.lammpslib import LAMMPSlib

        library_file = str(paths.static / "PtTi_Kim.meam.library")
        potential_file = str(paths.static / "PtTi_Kim.meam.potential")
        cmds = [
            "pair_style meam",
            f"pair_coeff * *  {library_file} Ni Ti {potential_file} Ni Ti",
        ]
        amds = ["thermo_style custom etotal lx ly lz vol pxx pyy pzz pxy pxz pyz press"]
        asecalc = LAMMPSlib(lmpcmds=cmds, amendments=amds, keep_alive=True)

    elif model == "M3GNet":
        import matgl
        from matgl.ext.ase import M3GNetCalculator

        m3gnet_params = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        asecalc = M3GNetCalculator(m3gnet_params, stress_weight=1.0 * GPa)

    elif model == "CHGNet":
        from chgnet.model.dynamics import CHGNetCalculator

        asecalc = CHGNetCalculator(stress_weight=1.0 * GPa)

    elif model == "MACE":
        from mace.calculators import MACECalculator

        #chkpoint = "2023-08-14-mace-universal.model"
        chkpoint = "2023-12-03-mace-128-L1_epoch-199.model"
        asecalc = MACECalculator(
            model_paths=str(paths.static / chkpoint),
            device="cpu",
            default_dtype="float32",
        )

        #Alternative API for newer MACE
        #from mace.calculators import mace_mp
        #asecalc = mace_mp() # return the default medium ASE calculator equivalent to mace_mp(model="medium")
        #asecalc = mace_mp(model="large") # return a larger model

    elif model == "ALIGNN":
        from alignn.ff.ff import AlignnAtomwiseCalculator, default_path

        alignn_params = default_path()
        asecalc = AlignnAtomwiseCalculator(path=alignn_params, device="cpu",stress_wt = 1.0 * GPa)

    elif model == "DeepMD":
        from deepmd.calculator import DP

        model_params = str(paths.static / "NiTi_DeepMD.pb")
        asecalc = DP(model=model_params)
        
    return asecalc


if __name__ == "__main__":
    calculator = get_ase_calculator()
    assert isinstance(calculator, CustomEAM)
