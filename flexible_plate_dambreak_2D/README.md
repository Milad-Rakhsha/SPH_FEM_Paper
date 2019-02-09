github.com/erdc-cm/proteus
branch: add_fem_dem
commit: 3d8077bb27729b2933cf1d439b8cd2ed1b94e908
run with
parun -l 5 -v dambreak_so.py


Important note:
This test case has not been fully tested for parallel processing.
#mpirun -n 8 parun -l 6 -v dambreak_so.py -D results -O ~/Repos/mrInternship/inputs/petsc.options.asm
parun -l 6 -v dambreak_so.py -D results -O ~/Repos/mrInternship/inputs/petsc.options.asm


