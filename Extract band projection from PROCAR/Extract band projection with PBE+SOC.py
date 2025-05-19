#Contributor:  Fulong Hu
#First version : 2023/5/17
# E-mail : 1500148414@qq.com (Fulong Hu)

import linecache
import numpy as np
import sys
import pandas as pd
from matplotlib import pyplot as plt

def Cal_inverse_space_vector():
    datas = np.loadtxt('POSCAR',dtype=np.float64,skiprows=2,max_rows=3,usecols=(0,1,2))
    unit_bohr = 1
    datas_bohr = datas/unit_bohr
    a1 = datas_bohr[0,:]
    a2 = datas_bohr[1,:]
    a3 = datas_bohr[2,:]
    V = np.dot(np.cross(a1,a2),a3)
    b1 = 2*np.pi*np.cross(a2,a3)/V
    b2 = 2*np.pi*np.cross(a3,a1)/V
    b3 = 2*np.pi*np.cross(a1,a2)/V
    b4 = [b1,b2,b3]
    b = np.matrix(b4)
    return b, V
b, V = Cal_inverse_space_vector()

def FERMI_ENERGY():
    Ef = linecache.getline('FERMI_ENERGY',2).split()
    E_f = float(Ef[0])
    return E_f
#E_f = float(input('请输入费米能级(单位eV):'))
    
def obtain_kpoints():
    info = linecache.getline('PROCAR',2).split()
    nk = int(info[3])  # the number of k-points
    nb=int(info[7])  # the number of bands
    nion=int(info[11]) # the number of ions
    print("The computing system contains", nion, 'atoms')    
    m = 4 #for including soc is equal to 4,else is equal to 1  
    line_spin_up=(((nion+1)*m+4)*nb+3)*nk+1
    if nion==1:
        nion=0
    x = 0.0
    X = []
    k0=linecache.getline('PROCAR',4).replace('-',' -').split()
    k=np.array([[float(k0[4]),float(k0[5]),float(k0[6])]])
    for kp in range(0,nk):
        k_line=4+kp*(((nion+1)*m+4)*nb+3)
        kpoint=linecache.getline('PROCAR',k_line)
        kpoint=kpoint.replace('-',' -').split()
        k = np.append(k,[[float(kpoint[4]),float(kpoint[5]),float(kpoint[6])]],axis=0) #The information about all k-points.
        dk=(k[1]-k[0])*b
        dk=np.sqrt(np.dot(dk,np.transpose(dk)))
        if kp in same_high_kpoints_number:   #Number of repeated points
            x = x
        else:
            x = x+dk[0,0]
        k=np.delete(k,0,axis=0)
        X.append(float('%.5f'%x))
    new_X = X * nb
    np.savetxt(r'kpoint.txt',X,fmt='%.5f')
    return nion, nk, m, nb

def KLABELS_and_same_high_kpoints_number(b):
    with open('KPATH.in','r') as reader:
        lines=reader.readlines()[0:]
    kpoint_name = []
    for i in range(4,len(lines)):
        if i%3 == 0:
            pass
        else:
            kpoint_name.append(lines[i].split()[3])
    same_high_kpoints_number = []
    for i in range((len(kpoint_name)//2)-1):
        if kpoint_name[2*i+1] == kpoint_name[2*(i+1)]:
            pass
        else:
            same_high_kpoints_number.append(int(2*(i+1)))
    datas_KPATH = np.loadtxt('KPATH.in',dtype=np.float64,skiprows=(4),usecols=(0,1,2))
    x = 0.0
    X, new_X = [], []
    k=np.array([[float(datas_KPATH[0,0]),float(datas_KPATH[0,1]),float(datas_KPATH[0,2])]])
    for kp in range(0,len(datas_KPATH[:,0])):
        kpoint=datas_KPATH[kp,:]
        k = np.append(k,[[float(kpoint[0]),float(kpoint[1]),float(kpoint[2])]],axis=0) #The information about all k-points.
        dk=(k[1]-k[0])*b
        dk=np.sqrt(np.dot(dk,np.transpose(dk)))
        if kp in same_high_kpoints_number:
            x = x
        else:
            x = x + dk[0,0]
        k = np.delete(k,0,axis=0)
        X.append(float('%.5f'%x))
    np.savetxt(r'kpoint.txt',X,fmt='%.5f')
    Hight_pointname, Hight_pointname1 = [], []
    for i in range(len(datas_KPATH[:,0])):
        if i%2 == 0:
            j = 5 + (i//2)*3
        else:
            j = 6 + ((i-1)//2)*3
        info = linecache.getline('KPATH.in',j).split()
        Hight_pointname.append(info[3])
    for i in range(1+len(X)//2):
        if i == 0:
            Hight_pointname1.append(Hight_pointname[i])
            new_X.append(X[i])
        elif i == len(X)//2:
            Hight_pointname1.append(Hight_pointname[len(X)-1])
            new_X.append(X[len(X)-1])
        else:
            if Hight_pointname[2*i] == Hight_pointname[2*i-1]:
                Hight_pointname1.append(Hight_pointname[2*i])
            else:
                new_name = Hight_pointname[2*i-1] + '|' + Hight_pointname[2*i]
                Hight_pointname1.append(new_name)
            new_X.append(X[2*i])

    
    X = new_X
    data1 = pd.DataFrame({'X': Hight_pointname1,'Y': X})
    X = np.transpose(np.array(X))

    filename = "KLABELS"
    with open(filename, 'w') as file:
        pass
    file.close()
    space_string = "                  "
    with open('KLABELS', mode = 'a+', newline='\n') as file:
        file.write('K-Label    K-Coordinate in band-structure plots \n')
        for i in range(len(X)):
            file.write(data1['X'][i])
            file.write(space_string)
            M = str(X[i])
            file.write(M)
            file.write('\n')
        file.write('\n')
        file.write('\n')
        file.write('* Give the label for each high symmetry point in KPOINTS (KPATH.in) file. Otherwise, they will be identified as ''Undefined'' in KLABELS file')
    return(same_high_kpoints_number)

def form_energy_and_projection_file():
    nion, nk, m, nb = obtain_kpoints()
    E_f = FERMI_ENERGY()
    En_number, En, En_total, En_real, En_rubsh, projection_of_atom_number, Nubern =[], [], [], [], [], [], []    
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is between ",'1',' to ',nion,", it indicates the total projection contribution of the energy band on the ",'1',' to ',nion,' atoms')
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is ",nion+1," it indicates the total projection contribution of the energy band as a whole.")
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is between ",(nion+1)+1,'to',2*(nion+1)-1,"it means the x-direction projection contribution of the energy band on the",'1','to',nion,' atoms')
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is ",2*(nion+1)," it indicates the total projection contribution of the energy band as a whole in the x direction")
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is between ",2*(nion+1)+1,'to',3*(nion+1)-1,"it means the y-direction projection contribution of the energy band on the",'1','to',nion,' atoms')
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is ",3*(nion+1)," it indicates the total projection contribution of the energy band as a whole in the y direction")
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is between ",3*(nion+1)+1,'to',4*(nion+1)-1,"it means the z-direction projection contribution of the energy band on the",'1','to',nion,' atoms')
    print("The atomic proportion indicates the projection of the specific atom in the output. When the input value is ",4*(nion+1)," it indicates the total projection contribution of the energy band as a whole in the z direction")
    atom_number = int(input('Please enter the atomic proportion to be calculated:'))    
    if atom_number in range((nion+1)+1,2*(nion+1)+1):
        last_name = 'X'
    elif atom_number in range(2*(nion+1)+1,3*(nion+1)+1):
        last_name = 'Y'
    else:
        last_name = 'Z'
    for kp in range(0,nk):
        k_line=4+kp*(((nion+1)*m+4)*nb+3)
        for kb in range(0,nb):
            line_band=2+k_line+(4+m*(nion+1))*kb
            line_projection_of_atom=4+atom_number+k_line+(4+m*(nion+1))*kb
            projection_of_atom_number.append(line_projection_of_atom)
            En_number.append(line_band)

    for i in range(0,len(En_number)):
        Band=linecache.getline('PROCAR',En_number[i]).split()
        en=float(Band[4])-E_f
        En.append(float('%.6f'%en))
    np.savetxt(r'coefficient.txt',En,fmt='%.6f')
    atom_contribution_set = int(input('Please enter the contribution of the part to be calculated, where 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 correspond to the contribution of s, py, pz, px, dxy, dyz, dz2, dxz, x2-y2, tot respectively.'))
    for i in range(0,len(projection_of_atom_number)):
        projection_of_atom=linecache.getline('PROCAR',projection_of_atom_number[i]).split()
        atoms_number=float(projection_of_atom[atom_contribution_set]) 
        Nubern.append(float('%.6f'%atoms_number))
    np.savetxt(r'projection.txt',Nubern,fmt='%.3f')

    a=[]
    a_projection =[]
    for i in range(0,nb):
        for j in range(0,nk):
            number = i + j*nb
            a.append(En[number])
            a_projection.append(Nubern[number])
    m = nb
    n = len(a)//m
    chunks = [a[i:i+n] for i in range(0,len(a),n)]
    chunks1 = [a_projection[i:i+n] for i in range(0,len(a),n)]
    chunks = np.matrix(chunks)
    chunks1 = np.matrix(chunks1)
    chunks = np.mat(chunks.T)
    chunks1 = np.mat(chunks1.T)
    np.savetxt(r'banddat1.txt',chunks,fmt='%.6f')
    np.savetxt(r'projection'+last_name+'.txt',chunks1,fmt='%.3f')
    return


if __name__ == "__main__":
    b,V = Cal_inverse_space_vector()
    same_high_kpoints_number = KLABELS_and_same_high_kpoints_number(b)
    form_energy_and_projection_file()
    print(same_high_kpoints_number)
