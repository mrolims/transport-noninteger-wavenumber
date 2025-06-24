import os
import sys

compiler = "ifx"

if len(sys.argv) > 1:
    compiler = sys.argv[1]

if compiler not in ["gfortran", "ifx"]:
    raise ValueError("Invalid compiler. Use 'gfortran' or 'ifx'.")
if compiler == "ifx":
    compiler = "ifx -fp-model precise -ipo"
else:
    compiler = "gfortran -flto"

file_prefix = "manifolds"
f90_file = file_prefix + ".f90"
exe_file = file_prefix + ".x"

os.system("rm -rf %s" % exe_file)
print("Compiling the program...")
comm = f"{compiler} params_qp.f90 functions.f90 {f90_file} -o {exe_file}"

print("$", comm)
os.system(comm)
if os.path.isfile(exe_file):
    print("Compilation successful.")
else:
    print("Problems in the compilation. Stopping execution...")
    sys.exit()


num_ic = int(1e5)
ms = [1, -1, 2, -2, 3, -3, 4, -4, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 3.4, -3.6, 3.5, -9]

for m in ms:

    comm = "./%s %.5f %i &" % (exe_file, m, num_ic)
    print("$", comm)
    os.system(comm)
