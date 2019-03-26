#this is python, so you can loop over the particles
# for i in range(nParticles): ...["q:('phis',{i})".format(i)]
#import pdb
#pdb.set_trace()
#note: this has to be attached to correct model
#navier-stokes model has index 4 in this model
#from so file:
    # pnList = [("vof_p",               "vof_n"),#0
    #           ("ls_p",                "ls_n"),#1
    #           ("redist_p",            "redist_n"),#2
    #           ("ls_consrv_p",         "ls_consrv_n"),#3
    #           ("twp_navier_stokes_p", "twp_navier_stokes_n"),#4
    #           ("pressureincrement_p", "pressureincrement_n"),#5
    #           ("pressure_p", "pressure_n"),#6
    #           ("pressureInitial_p", "pressureInitial_n")]#7

#simFlagsList[4]['storeQuantities'] = ["q:('phis')"]
simFlagsList[4]['storeQuantities'] = ["q:('cfl', 0)"]

start
end
