"""

In order to create the global stiffness matrix (denoted KG) , first the connectivity matrix must be created.
If the global displacement vector of an element ELE_i is U_ELE_i = [u1 , u2 ... , un] , where n the DOFs of said
element , the row of the connectivity matrix corresponding to ELE_i MUST be IX(ELE_i) = [ 1 , 2 , ... , n].

"""



import numpy as np

NELE = 10
NNOD = 6
NDOFS = 2*NNOD

E1 = 200e9
E2 = 200e9
E3 = 200e9
A1 = 10e-2
A2 = 10e-2
A3 = 10e-2
L1 = 1
L2 = 1
L3 = 1
th5 = np.arctan(360/280)
th4 = 0
th1 = 0
th2 = 0
th3 = 0
th6 = np.pi/2
th9 = np.pi/2
th10 = np.pi - th5
th7 = np.arctan(3600/500)
th8 = np.pi - th7
th = [th1,th2,th3,th4,th5,th6,th7,th8,th9,th10]

Kloc = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(NELE)]

for thi,i in zip(th,range(NELE)):
    c = np.cos(thi)
    s = np.sin(thi)
    Kloc[i] = [[c**2,c*s,-c**2,-c*s],[c*s,s**2,-c*s,-s**2],[-c**2,-c*s,c**2,c*s],[-c*s,-s**2,c*s,s**2]]

for ele in range(NELE):
    for i in range(4):
        for j in range(4):
            Kloc[ele][i][j] *= 2E9

KG = [[0 for _ in range(NDOFS)] for _ in range(NDOFS)]

IX = [[1,2,5,6],
      [5,6,9,10],
      [9,10,11,12],
      [3,4,7,8],
      [1,2,3,4],
      [5,6,3,4],
      [5,6,7,8],
      [9,10,3,4],
      [9,10,7,8],
      [11,12,7,8]]

print(len(IX))
for i in range(len(IX)):
    for j in range(len(IX[0])):
        IX[i][j] -=1


for ele in range(NELE):
    for i in range(4):
        for j in range(4):
            KG[IX[ele][i]][IX[ele][j]] +=  Kloc[ele][i][j]

F = [0 for _ in range(NDOFS)]

#####ENTER CONSTRAINTS/KNOWNS##########
########################################
F[2] = 40e3
F[5] = 20e3
F[8] = 20e3

FREE_NODES = [3-1,4-1,7-1,8-1,5-1,6-1,9-1,10-1,11-1]
KNOWN_FORCES = FREE_NODES
########################################

UDOFs = len(FREE_NODES) #UNCONSTRAINED NODES
KU = [[0 for _ in range(UDOFs)] for _ in range(UDOFs)]
FK = [F[i] for i in KNOWN_FORCES]

count_r = 0
count_c = 0

for i in KNOWN_FORCES:
    for j in FREE_NODES:
        KU[count_r][count_c] = KG[i][j]
        count_c += 1
    count_c = 0
    count_r +=1

UU = np.linalg.solve(KU,FK)

U = [0 for _ in range(NDOFS)]
count = 0
for i in FREE_NODES:
    U[i] = UU[count]
    count +=1

Forces = np.dot(np.array(KG),np.array(U))
Disps  = np.array(U)
print(np.round(Forces))
print(Disps)