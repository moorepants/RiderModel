from RiderLeanAB import RiderLeanAB
import pickle
from scipy.io.matlab.mio import savemat

ppdir = '/media/Data/Documents/School/TU Delft/PhysicalParameters/'

# load in the base data file
f = open(ppdir + 'data/data.p', 'r')
data = pickle.load(f)
f.close()

v = [2.5, 5., 7.5]
# write the parameter files
for i, name in enumerate(data['bikes']):
    fname = ''.join(name.split())
    directory = '../data/BikeLeanAB/'
    f = open(directory + fname + 'LeanAB.txt', 'w')
    f.write('States (x) are: [Roll,Steer,Lean,Roll Rate,Steer Rate,Lean Rate]\n')
    f.write('Inputs (u) are: [Roll Torque,Steer Torque,Lean Torque]\n')
    for speed in v:
        A, B, nothing = RiderLeanAB(i, speed)
        mdict = {'A':A, 'B':B}
        savemat(directory + 'mat/' + fname + 'LeanAB_' + str(speed) + '.mat', mdict)
        f.write('A' + ' (v = ' + str(speed) + ')\n')
        f.write(str(A) + '\n')
        f.write('B' + ' (v = ' + str(speed) + ')\n')
        f.write(str(B) + '\n')
    f.close()
