from subprocess import call
import os
p = os.path.dirname(os.path.abspath(__file__))
call(['make'], cwd=p)
try:
    from . import peridynamics_lib
except ImportError:
    call(['rm','-r','CMakeFiles'], cwd=p)
    call(['rm','CMakeCache.txt'], cwd=p)
    call(['cmake',p], cwd=p)
    call(['make'],  cwd=p)
    from . import peridynamics_lib
kernel_Fbased_inv = peridynamics_lib.cvar.kernel_Fbased_inv
kernel_Fbased_linear = peridynamics_lib.cvar.kernel_Fbased_linear
kernel_Oterkus2_quadr = peridynamics_lib.cvar.kernel_Oterkus2_quadr
kernel_Fbased_cubic = peridynamics_lib.cvar.kernel_Fbased_cubic
kernel_Fbased_quadr = peridynamics_lib.cvar.kernel_Fbased_quadr
kernel_Fstab_Silling_quarticA = peridynamics_lib.cvar.kernel_Fstab_Silling_quarticA
kernel_Silling_quarticA = peridynamics_lib.cvar.kernel_Silling_quarticA
kernel_Oterkus2_linear = peridynamics_lib.cvar.kernel_Oterkus2_linear
kernel_Silling_cubic = peridynamics_lib.cvar.kernel_Silling_cubic
kernel_Fstab_Silling_cubic = peridynamics_lib.cvar.kernel_Fstab_Silling_cubic
kernel_smooth_inv = peridynamics_lib.cvar.kernel_smooth_inv
kernel_smooth_linear = peridynamics_lib.cvar.kernel_smooth_linear
kernel_Fstab_Littlewood_cubic = peridynamics_lib.cvar.kernel_Fstab_Littlewood_cubic
kernel_Oterkus2_inv = peridynamics_lib.cvar.kernel_Oterkus2_inv
kernel_Fbased_quarticA = peridynamics_lib.cvar.kernel_Fbased_quarticA
kernel_Fstab_Silling_linear = peridynamics_lib.cvar.kernel_Fstab_Silling_linear
kernel_smooth_const = peridynamics_lib.cvar.kernel_smooth_const
kernel_Silling_linear = peridynamics_lib.cvar.kernel_Silling_linear
kernel_Fstab_Littlewood_quarticA = peridynamics_lib.cvar.kernel_Fstab_Littlewood_quarticA
kernel_Oterkus2_quarticA = peridynamics_lib.cvar.kernel_Oterkus2_quarticA
kernel_smooth_quarticA = peridynamics_lib.cvar.kernel_smooth_quarticA
kernel_Fstab_Littlewood_const = peridynamics_lib.cvar.kernel_Fstab_Littlewood_const
kernel_Fstab_Silling_const = peridynamics_lib.cvar.kernel_Fstab_Silling_const
kernel_Oterkus2_cubic = peridynamics_lib.cvar.kernel_Oterkus2_cubic
kernel_Silling_inv = peridynamics_lib.cvar.kernel_Silling_inv
kernel_Fstab_Littlewood_quadr = peridynamics_lib.cvar.kernel_Fstab_Littlewood_quadr
kernel_Fbased_const = peridynamics_lib.cvar.kernel_Fbased_const
kernel_Fstab_Littlewood_linear = peridynamics_lib.cvar.kernel_Fstab_Littlewood_linear
kernel_Fstab_Silling_quadr = peridynamics_lib.cvar.kernel_Fstab_Silling_quadr
kernel_Silling_quadr = peridynamics_lib.cvar.kernel_Silling_quadr
kernel_Fstab_Silling_inv = peridynamics_lib.cvar.kernel_Fstab_Silling_inv
kernel_Oterkus2_const = peridynamics_lib.cvar.kernel_Oterkus2_const
kernel_smooth_cubic = peridynamics_lib.cvar.kernel_smooth_cubic
kernel_Fstab_Littlewood_inv = peridynamics_lib.cvar.kernel_Fstab_Littlewood_inv
kernel_Silling_const = peridynamics_lib.cvar.kernel_Silling_const
kernel_smooth_quadr = peridynamics_lib.cvar.kernel_smooth_quadr
