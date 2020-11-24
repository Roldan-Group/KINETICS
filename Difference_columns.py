


import sys
import numpy as np

file_name = [ str(i) for i in sys.argv[1:]]
output = file_name[0][:-4]
for i in file_name[1:]:
    output = str(output + "_" + i[:-4])

if len(file_name) > 2:
    print("   Only accept two files")
    sys.exit()


#for name in file_name[1:]:
#    if len(np.loadtxt(name, comments="#")) != len(np.loadtxt(file_name[0], comments="#")):
#        print("   Files of different lenght")
#        break



X = np.loadtxt(file_name[0])[:,0]
Y = np.loadtxt(file_name[0])[:,1]

DataSet = np.loadtxt(file_name[1], comments="#")                            # import and reads data
for i in range(len(DataSet)):
    if DataSet[i,0] not in X:
        X.append(DataSet[i,0])
        Y.append(DataSet[i,1]*(-1))
    else:
        if DataSet[i,1] != 0:
            Y[i] = (Y[i] - DataSet[i,1])


XY = [ (X[i], Y[i]) for i in range(len(X)) ]
XY = sorted(XY, key=lambda x: x[0])


ofile = open("DIFF_" + output + file_name[0][-4:],'w+')
ofile.write ("# Second column difference of ")
for name in file_name:
    ofile.write (str(name))
ofile.write ("\n")
for i in range(len(DataSet)):
    ofile.write ("%.5f %.5f\n" %(XY[i]))
ofile.close()

