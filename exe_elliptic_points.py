import os
import sys

compiler = "gfortran"

if len(sys.argv) > 1:
    compiler = sys.argv[1]

if compiler not in ["gfortran", "ifx"]:
    raise ValueError("Invalid compiler. Use 'gfortran' or 'ifx'.")
if compiler == "ifx":
    compiler = "ifx -fp-model precise -ipo"
else:
    compiler = "gfortran -flto"

file_prefix = "elliptic_points"
f90_file = file_prefix + ".f90"
exe_file = file_prefix + ".x"

m = [1, -1, 2, -2, 3, -3, 4, -4, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 3.4, -3.6, 3.5, -9]

x_range = {"upper": [0.47, 0.53], "lower": [0.51, 0.55]}
y_range = {"upper": [0.3, 0.4], "lower": [0.05, 0.2]}

os.system("rm -rf %s" % exe_file)
print("Compiling the program...")
comm = f"{compiler} params_qp.f90 functions.f90 {f90_file} -o {exe_file}" % (f90_file, exe_file)

comm = f"{compiler} params_dp.f90 functions.f90 %s -o %s" % (f90_file, exe_file)
print("$", comm)
os.system(comm)
if os.path.isfile(exe_file):
    print("Compilation successful.")
else:
    print("Problems in the compilation. Stopping execution...")
    sys.exit()

print("sbatch script generated successfully.")
print("Executing script...")
for region in ["upper", "lower"]:
    for i in range(len(m)):
        comm = "./%s %.5f %.5f %.5f %.5f %.5f 0 %s" % (exe_file, m[i], x_range[region][0], x_range[region][1], y_range[region][0], y_range[region][1], region)
        print("$", comm)
        os.system(comm)
print("-" * len(comm))
