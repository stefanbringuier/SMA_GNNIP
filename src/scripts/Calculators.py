import paths
import numpy as np
from ase.calculators.eam import EAM


class CustomEAM(EAM):
    '''
    Override EAM to zero out stress for ASE Database support
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def get_stress(self, atoms=None, **kwargs):
        return np.zeros(6)

def GetCalculator(args):
    if args.model == 'Mutter':
        asecalc = CustomEAM(potential= str(paths.static / "NiTi_mutter.eam.fs"))
    elif args.model == 'Zhong':
        asecalc = CustomEAM(potential= str(paths.static / "NiTiZhong.eam.fs"))
    elif args.model == "M3GNet":
        import matgl
        from matgl.ext.ase import M3GNetCalculator
        m3gnet_params = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        asecalc = M3GNetCalculator(m3gnet_params, stress_weight=1.0)
    elif args.model == "CHGNet":
        from chgnet.model.dynamics import CHGNetCalculator
        asecalc = CHGNetCalculator(stress_weight=1.0)
    elif args.model == "MACE":
        from mace.calculators import MACECalculator
        asecalc = MACECalculator(model_path= str(paths.static / '2023-08-14-mace-universal.model'),
                      device="cpu",
                      default_dtype="float32")
    elif args.model == "ALIGNN":
        from alignn.ff.ff import AlignnAtomwiseCalculator, default_path
        alignn_params = default_path()
        asecalc = AlignnAtomwiseCalculator(path=alignn_params, device='cpu')

    return asecalc
