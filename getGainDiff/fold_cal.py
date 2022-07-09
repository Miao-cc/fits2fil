import sys
import glob
import math
import numpy as np
import matplotlib.pyplot as plt

def getNum(filename):
    filenum = filename.split("M01_")[-1].split(".fits")[0]
    return int(filenum)

tsamp = float(sys.argv[1])*1E-6
fold_period = float(sys.argv[2])

fold_step = int(fold_period/tsamp)

print(fold_period, tsamp, fold_period/tsamp)

pol_AA_list = sorted(glob.glob("*AA.txt"), key=getNum)
pol_BB_list = sorted(glob.glob("*BB.txt"), key=getNum)
pol_AA = np.array([])
pol_BB = np.array([])
for pol_file in pol_AA_list:
    pol_AA = np.append(pol_AA, np.loadtxt(pol_file))

for pol_file in pol_BB_list:
    pol_BB = np.append(pol_BB, np.loadtxt(pol_file))

fold_round = math.floor(len(pol_AA)/fold_step)
print(pol_AA.shape, len(pol_AA))
print(pol_BB.shape, len(pol_BB))

plt.plot(pol_AA[:fold_step], 'r')
plt.plot(pol_BB[:fold_step], 'b')
plt.show()

pol_AA = pol_AA[:fold_round*fold_step].reshape(-1, fold_step).mean(axis=0)
pol_BB = pol_BB[:fold_round*fold_step].reshape(-1, fold_step).mean(axis=0)

plt.plot(pol_AA, 'r')
plt.plot(pol_BB, 'b')
plt.show()

onCal_idx = np.logical_and(pol_AA > np.mean(pol_AA), pol_BB > np.mean(pol_BB))
offCal_idx = np.logical_and(pol_AA < np.mean(pol_AA), pol_BB < np.mean(pol_BB))

onCal_mean_AA = np.mean(pol_AA[onCal_idx])
offCal_mean_AA = np.mean(pol_AA[offCal_idx])

onCal_mean_BB = np.mean(pol_BB[onCal_idx])
offCal_mean_BB = np.mean(pol_BB[offCal_idx])
print(onCal_mean_AA, offCal_mean_AA, onCal_mean_AA-offCal_mean_AA)
print(onCal_mean_BB, offCal_mean_BB, onCal_mean_BB-offCal_mean_BB)

gainDiff = (onCal_mean_AA-offCal_mean_AA) / (onCal_mean_BB-offCal_mean_BB)
baseLineDiff = offCal_mean_AA - offCal_mean_BB*gainDiff
#print(onCal_mean_AA, offCal_mean_AA, onCal_mean_AA-offCal_mean_AA)
#print(onCal_mean_BB, offCal_mean_BB, onCal_mean_BB-offCal_mean_BB)
plt.plot(pol_AA, 'r')
plt.plot(pol_BB*gainDiff+baseLineDiff, 'b')
plt.show()
print("gainAA/gainBB: ", gainDiff, baseLineDiff)
