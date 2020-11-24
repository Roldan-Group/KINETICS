


import sys
import numpy as np

file_name = [ str(i) for i in sys.argv[1:]]
output = file_name[0][:-4]
for i in file_name[1:]:
    output = str(output + "_" + i[:-4])

#for name in file_name[1:]:
#    if len(np.loadtxt(name, comments="#")) != len(np.loadtxt(file_name[0], comments="#")):
#        print("   Files of different lenght")
#        break

X = np.loadtxt(file_name[0])[:,0]
Y = [0]*len(np.loadtxt(file_name[0])[:,0])
for n, name in enumerate(file_name):
    DataSet = np.loadtxt(name, comments="#")                            # import and reads data
    for i in range(len(DataSet)):
        if DataSet[i,0] not in X:
            X.append(DataSet[i,0])
            Y.append(DataSet[i,1])
        else:
            if DataSet[i,1] != 0:
                Y[i] = (Y[i] + DataSet[i,1])/2


XY = [ (X[i], Y[i]) for i in range(len(X)) ]
XY = sorted(XY, key=lambda x: x[0])


ofile = open("SUM_" + output + file_name[0][-4:],'w+')
ofile.write ("# Second column sum of ")
for name in file_name:
    ofile.write (str(name))
ofile.write ("\n")
for i in range(len(DataSet)):
    ofile.write ("%.5f %.5f\n" %(XY[i]))
ofile.close()

