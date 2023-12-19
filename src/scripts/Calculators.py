import paths
import numpy as np
from ase.calculators.eam import EAM


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

    elif model == "M3GNet":
        import matgl
        from matgl.ext.ase import M3GNetCalculator

        m3gnet_params = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        asecalc = M3GNetCalculator(m3gnet_params, stress_weight=1.0)

    elif model == "CHGNet":
        from chgnet.model.dynamics import CHGNetCalculator

        asecalc = CHGNetCalculator(stress_weight=1.0)

    elif model == "MACE":
        from mace.calculators import MACECalculator

        asecalc = MACECalculator(
            model_path=str(paths.static / "2023-08-14-mace-universal.model"),
            device="cpu",
            default_dtype="float32",
        )

    elif model == "ALIGNN":
        from alignn.ff.ff import AlignnAtomwiseCalculator, default_path

        alignn_params = default_path()
        asecalc = AlignnAtomwiseCalculator(path=alignn_params, device="cpu")

    elif model == "Ko":
        from ase.calculators.lammpslib import LAMMPSlib

        library_file = str(paths.static / "NiTi_Ko.meam.library")
        potential_file = str(paths.static / "NiTi_Ko.meam.potential")
        cmds = ["pair_style meam", f"pair_coeff * *  {library_file} Ni Ti {potential_file} Ni Ti"]
        amds = ["thermo_style custom etotal lx ly lz vol pxx pyy pzz pxy pxz pyz press"]
        asecalc = LAMMPSlib(lmpcmds=cmds, amendments=amds,  log_file='test.log',keep_alive=True)
        
    return asecalc


if __name__ == "__main__":
    calculator = get_ase_calculator()
    assert isinstance(calculator, CustomEAM)
