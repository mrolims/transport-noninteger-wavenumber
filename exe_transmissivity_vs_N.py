import os
import sys
import numpy as np

compiler = "ifx"

if len(sys.argv) > 1:
    compiler = sys.argv[1]

if compiler not in ["gfortran", "ifx"]:
    raise ValueError("Invalid compiler. Use 'gfortran' or 'ifx'.")
if compiler == "ifx":
    compiler = "ifx -ipo"
else:
    compiler = "gfortran -flto"

file_prefix = "transmissivity_vs_N"
f90_file = file_prefix + ".f90"
exe_file = file_prefix + ".x"

os.system("rm -rf %s" % exe_file)
print("Compiling the program...")
comm = f"{compiler} params_dp.f90 functions.f90 {f90_file} -o {exe_file}"
print("$", comm)
os.system(comm)
if os.path.isfile(exe_file):
    print("Compilation successful.")
else:
    print("Problems in the compilation. Stopping execution...")
    sys.exit()


num_ic = int(1e4)
y0s = [1, -1]
esc_ys = [-1, 1]
N = int(1e6)

golden_ratio = (1 - np.sqrt(5)) / 2  # example definition if needed

ms = [golden_ratio, -golden_ratio]

# list of numbers whose square roots we want
nums = [2.0, 3.0, 5.0, 6.0, 7.0, 8.0, 10.0, 11.0, 12.0, 13.0, 14.0]

# extend values with sqrt and -sqrt of each number
for n in nums:
    root = np.sqrt(n)
    ms.append(root)
    ms.append(-root)

for i in range(len(y0s)):
    y0 = y0s[i]
    esc_y = esc_ys[i]
    for m in ms:
        comm = f"./{exe_file} {m} {y0} {esc_y} {N}"
        print("$", comm)
        os.system(comm)
