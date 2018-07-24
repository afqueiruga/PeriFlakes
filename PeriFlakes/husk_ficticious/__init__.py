from subprocess import call
import os
p = os.path.dirname(os.path.abspath(__file__))
call(['make'], cwd=p)
try:
    from . import ficticious_lib
except ImportError:
    call(['rm','-r','CMakeFiles'], cwd=p)
    call(['rm','CMakeCache.txt'], cwd=p)
    call(['cmake',p], cwd=p)
    call(['make'],  cwd=p)
    from . import ficticious_lib
kernel_bobaru_n = ficticious_lib.cvar.kernel_bobaru_n
kernel_bobaru_F3 = ficticious_lib.cvar.kernel_bobaru_F3
kernel_bobaru_F = ficticious_lib.cvar.kernel_bobaru_F
kernel_bobaru_n3 = ficticious_lib.cvar.kernel_bobaru_n3
kernel_bobaru_y = ficticious_lib.cvar.kernel_bobaru_y
