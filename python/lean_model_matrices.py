from RiderLeanAB import RiderLeanAB
import pickle

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
    f.write('States (x) are: [Roll,Steer,Lean,Roll Rate,Steer Rate,Lean Rate]\n\r')
    f.write('Inputs (u) are: [Roll Torque,Steer Torque,Lean Torque]\n\r')
    for speed in v:
        A, B, nothing = RiderLeanAB(i, speed)
        f.write('A' + ' (v = ' + str(speed) + ')\n\r')
        f.write(str(A) + '\n\r')
        f.write('B' + ' (v = ' + str(speed) + ')\n\r')
        f.write(str(B) + '\n\r')
    f.close()
