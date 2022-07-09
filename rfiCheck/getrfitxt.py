import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

rfilist = ''
count = 0

if len(sys.argv[:]) == 1:
    filelist = glob.glob("*rfi-channel-list.txt")
elif len(sys.argv[:]) == 2:
    filelist = [sys.argv[1]]
else:
    sys.exit(0)

data = np.array([])
for filename in filelist:
    #print("loading file: ", filename)
    data = np.append(data, np.loadtxt(filename,dtype=int))

def runs(iterable):
    iterator = iter(iterable)
    start = end = next(iterator)

    for item in iterator:
        if item != end + 1:
            yield (start, end)
            start = item

        end = item

    yield (start, end)

data = np.array(list(set(data)), dtype=int)
lows, highs = zip(*runs(data))
#print(lows, highs)
for i,j in zip(lows, highs):
    #print(i,j)
    rfilist = rfilist+str(i)+':'+str(j)+','

print(rfilist[:-1])
