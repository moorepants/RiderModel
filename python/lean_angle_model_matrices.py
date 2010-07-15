from LeanAngleAB import LeanAngleAB
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
    directory = '../data/LeanAngleAB/'
    f = open(directory + fname + 'LeanAngleAB.txt', 'w')
    f.write('States (x) are: [Roll,Steer,Roll Rate,Steer Rate]\n')
    f.write('Inputs (u) are: [Roll Torque, Steer Torque, Lean Angle, Lean Rate, Lean Accel]\n')
    for speed in v:
        A, B, nothing = LeanAngleAB(i, speed)
        mdict = {'A':A, 'B':B}
        savemat(directory + 'mat/' + fname + 'LeanAngleAB_' + str(speed) + '.mat', mdict)
        f.write('A' + ' (v = ' + str(speed) + ')\n')
        f.write(str(A) + '\n')
        f.write('B' + ' (v = ' + str(speed) + ')\n')
        f.write(str(B) + '\n')
    f.close()
