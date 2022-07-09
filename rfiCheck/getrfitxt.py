import sys
import numpy as np
import matplotlib.pyplot as plt

rfilist = ''
count = 0


name = sys.argv[1]
data = np.loadtxt(name+'-rfi-channel-list.txt',dtype=int)


#plt.scatter(data, np.zeros(len(data)))
#
#for i in rfilist:
#    plt.axvline(x=i)
#plt.show()


def runs(iterable):
    iterator = iter(iterable)
    start = end = next(iterator)

    for item in iterator:
        if item != end + 1:
            yield (start, end)
            start = item

        end = item

    yield (start, end)


lows, highs = zip(*runs(data))
#print(lows, highs)
for i,j in zip(lows, highs):
    #print(i,j)
    rfilist = rfilist+str(i)+':'+str(j)+','

print(rfilist[:-1])
