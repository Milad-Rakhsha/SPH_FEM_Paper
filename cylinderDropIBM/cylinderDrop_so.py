from proteus.default_so import *
import cylinderDrop

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

if cylinderDrop.sedimentDynamics:
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
    cylinderDrop.VOS_model=0
    cylinderDrop.VOF_model=1
    cylinderDrop.LS_model=2
    cylinderDrop.RD_model=3
    cylinderDrop.MCORR_model=4
    cylinderDrop.SED_model=5
    cylinderDrop.V_model=6
    cylinderDrop.PINC_model=7
    cylinderDrop.PRESSURE_model=8
    cylinderDrop.PINIT_model=9
else:
    pnList = [("vof_p",               "vof_n"),#0
              ("ls_p",                "ls_n"),#1
              ("redist_p",            "redist_n"),#2
              ("ls_consrv_p",         "ls_consrv_n"),#3
              ("twp_navier_stokes_p", "twp_navier_stokes_n"),#4
              ("pressureincrement_p", "pressureincrement_n"),#5
              ("pressure_p", "pressure_n"),#6
              ("pressureInitial_p", "pressureInitial_n")]#7
    cylinderDrop.VOS_model=None
    cylinderDrop.SED_model=None
    cylinderDrop.VOF_model=0
    cylinderDrop.LS_model=1
    cylinderDrop.RD_model=2
    cylinderDrop.MCORR_model=3
    cylinderDrop.V_model=4
    cylinderDrop.PINC_model=5
    cylinderDrop.PRESSURE_model=6
    cylinderDrop.PINIT_model=7

if cylinderDrop.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "cylinderDrop"

#modelSpinUpList = [cylinder.VOF_model, cylinder.LS_model, cylinder.V_model, cylinder.PINIT_model]
modelSpinUpList = [cylinderDrop.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]


systemStepControllerType = Sequential_MinAdaptiveModelStepPS
# systemStepControllerType = Sequential_FixedStep #Sequential_FixedStep_Simple # uses time steps in so.tnList
# dt_system_fixed = cylinderDrop.dt_fixed ;
# systemStepExact=False;




needEBQ_GLOBAL = False
needEBQ = False
