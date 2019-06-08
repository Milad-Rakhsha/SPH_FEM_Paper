from proteus.default_so import *
import fsi_validation
from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

if fsi_validation.sedimentDynamics:
    pnList = [("vos_p",               "vos_n"),#0
              ("vof_p",               "vof_n"),#1
              ("ls_p",                "ls_n"),#2
              ("redist_p",            "redist_n"),#3
              ("ls_consrv_p",         "ls_consrv_n"),#4
              ("threep_navier_stokes_sed_p", "threep_navier_stokes_sed_n"),#5
              ("twp_navier_stokes_p", "twp_navier_stokes_n"),#6
              ("pressureincrement_p", "pressureincrement_n"),#7
              ("pressure_p", "pressure_n"),#8
              ("pressureInitial_p", "pressureInitial_n")]#9
    fsi_validation.VOS_model=0
    fsi_validation.VOF_model=1
    fsi_validation.LS_model=2
    fsi_validation.RD_model=3
    fsi_validation.MCORR_model=4
    fsi_validation.SED_model=5
    fsi_validation.V_model=6
    fsi_validation.PINC_model=7
    fsi_validation.PRESSURE_model=8
    fsi_validation.PINIT_model=9
else:
    pnList = [("vof_p",               "vof_n"),#0
              ("ls_p",                "ls_n"),#1
              ("redist_p",            "redist_n"),#2
              ("ls_consrv_p",         "ls_consrv_n"),#3
              ("twp_navier_stokes_p", "twp_navier_stokes_n"),#4
              ("pressureincrement_p", "pressureincrement_n"),#5
              ("pressure_p", "pressure_n"),#6
              ("pressureInitial_p", "pressureInitial_n")]#7
    fsi_validation.VOS_model=None
    fsi_validation.SED_model=None
    fsi_validation.VOF_model=0
    fsi_validation.LS_model=1
    fsi_validation.RD_model=2
    fsi_validation.MCORR_model=3
    fsi_validation.V_model=4
    fsi_validation.PINC_model=5
    fsi_validation.PRESSURE_model=6
    fsi_validation.PINIT_model=7

if fsi_validation.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "fsi_validation"

#modelSpinUpList = [fsi_validation.VOF_model, fsi_validation.LS_model, fsi_validation.V_model, fsi_validation.PINIT_model]
modelSpinUpList = [fsi_validation.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]


systemStepControllerType = Sequential_MinAdaptiveModelStepPS
# systemStepControllerType = Sequential_FixedStep#Sequential_FixedStep #Sequential_FixedStep_Simple # uses time steps in so.tnList
# dt_system_fixed = fsi_validation.dT_Chrono;
# systemStepExact=False;


# Time stepping/


tnList = fsi_validation.tnList
info = open("TimeList.txt","w")



needEBQ_GLOBAL = False
needEBQ = False
