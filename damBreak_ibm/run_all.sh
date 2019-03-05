mpirun -n 10 parun -v -l 5 dambreak_so.py -O petsc.options.asm -D res3 -C "use_r_ls_consrv=0 T=6.0"
mpirun -n 10 parun -v -l 5 dambreak_so.py -O petsc.options.asm -D res4 -C "use_r_ls_consrv=1 T=6.0"
mpirun -n 10 parun -v -l 5 dambreak_so.py -O petsc.options.asm -D res5 -C "use_r_ls_consrv=0 T=3.0"
mpirun -n 10 parun -v -l 5 dambreak_so.py -O petsc.options.asm -D res6 -C "use_r_ls_consrv=1 T=3.0"
