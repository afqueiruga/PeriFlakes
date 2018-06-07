from subprocess import call
import os
p = os.path.dirname(os.path.abspath(__file__))
call(['make'], cwd=p)
try:
    from . import bonds_lib
except ImportError:
    call(['rm','-r','CMakeFiles'], cwd=p)
    call(['rm','CMakeCache.txt'], cwd=p)
    call(['cmake',p], cwd=p)
    call(['make'],  cwd=p)
    from . import bonds_lib
kernel_line_intersection = bonds_lib.cvar.kernel_line_intersection
kernel_bond_pressure = bonds_lib.cvar.kernel_bond_pressure
