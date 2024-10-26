#   Fulong Hu
#First version : 2024/10/2
#Second version : 2024/10/26
#E-mail : 1500148414@qq.com


import linecache
import numpy as np
import cmath
import scipy.constants as C
import os

def select_input_filenames():
    current_folder = os.getcwd()
    filenames = os.listdir(current_folder)

    for i in range(len(filenames)):
        print('The', i, 'file name is','"',filenames[i],'"')
    print('    ')
    number = int(input("Please enter the input file number:",))
    select_filenames = filenames[number]
    return select_filenames

def solve_complex_equations(coefficients, constants):
    constants = [10**-12 if x == 0 else x for  x in constants]
    # Converting factors and constants to complex form
    coefficients = np.array(coefficients)
    constants = np.array(constants)
    sum_S2 = constants[0]**2 + constants[1]**2 + constants[2]**2
    if constants[0]>0:
        a_new = np.sqrt(constants[0]**2/sum_S2)
    else:
        a_new = -np.sqrt(constants[0]**2/sum_S2)

    if constants[1]>0:
        b_new = np.sqrt(constants[1]**2/sum_S2)
    else:
        b_new = -np.sqrt(constants[1]**2/sum_S2)

    if constants[2]>0:
        c_new = np.sqrt(constants[2]**2/sum_S2)
    else:
        c_new = -np.sqrt(constants[2]**2/sum_S2)

    constants = [a_new,b_new,c_new,1]
    
    # Solving a system of equations
    solutions = np.linalg.solve(coefficients, constants)
    β = cmath.sqrt(list(solutions)[0]*list(solutions)[3]/list(solutions)[1])
    α = cmath.sqrt(list(solutions)[2])   
    return α, β

#The matrix of the equation <Sx> = α∗β + β∗α,<Sy> = i(β∗α − α∗β),<Sz> = α^2 − β^2, α^2 + β^2=1, which are listed according to the values ​​of Sx, Sy and Sz and the normalization conditions
coefficients = [[1, 1, 0, 0],
                [-1j, 1j, 0, 0],
                [0, 0, 1, -1],
                [0, 0, 1, 1]]

def calculate_spin_overlap(file_path):
    f = open('Results.txt','w')
    print('The source code in Python is free to use and to improve, where only request is that this work is cited.',file=f)
    info=linecache.getline(file_path,2).split()
    iscalculated_with_both_cbm_and_vbm = input("Please select whether to calculate spin overlap—— Y(yes) or N(no)? \n")
    if iscalculated_with_both_cbm_and_vbm not in ['Y','y','N','n']:
        print("Error, please input E, WL or WN.")
        exit()
    elif iscalculated_with_both_cbm_and_vbm in ['Y','y']:
        for i in range(1,len(info)):
            data_line0 = 3*i
            data_line1 = 3*i + 1
            data_line2 = 3*i + 2
            point = linecache.getline(file_path,data_line0).split()
            S_i = linecache.getline(file_path,data_line1).split()
            S_j = linecache.getline(file_path,data_line2).split()
            constants = [float(S_i[0]),float(S_i[1]),float(S_i[2])]
            constants1 = [float(S_j[0]),float(S_j[1]),float(S_j[2])]
            solutions = solve_complex_equations(coefficients, constants)
            solutions1 = solve_complex_equations(coefficients, constants1)
            α = solutions[0]
            β = solutions[1]
            
            α1 = solutions1[0]
            β1 = solutions1[1]
            
            β1_conj = β1.conjugate()
            α1_conj = α1.conjugate()
            results = α1_conj * α + β1_conj * β
            spin_overlap = (results.real**2 + results.imag**2)*100
            
            print(point,file=f)
            print('The spinor wave function (α) of this material at CBM is','%.3f'%α.real,'+','%.3f'%α.imag,'i','(','%.1f'%(100*(α.real**2+α.imag**2)),'%',')',file=f)
            print('The spinor wave function (β) of this material at CBM  is','%.3f'%β.real,'+','%.3f'%β.imag,'i','(','%.1f'%(100*(β.real**2+β.imag**2)),'%',')',file=f)
            print('The spinor wave function (α) of this material at VBM is','%.3f'%α1.real,'+','%.3f'%α1.imag,'i','(','%.1f'%(100*(α1.real**2+α1.imag**2)),'%',')',file=f)
            print('The spinor wave function (β) of this material at VBM is','%.3f'%β1.real,'+','%.3f'%β1.imag,'i','(','%.1f'%(100*(β1.real**2+β1.imag**2)),'%',')',file=f)
            print('The spin overlap of the material is equal to：', '%.2f'%spin_overlap,'%',file=f)
            print('',file=f)
    else:
        for i in range(1,len(info)):
            data_line1 = 2*i + 2
            data_line0 = 2*i + 1
            point = linecache.getline(file_path,data_line0).split()
            S_i = linecache.getline(file_path,data_line1).split()
            constants = [float(S_i[0]),float(S_i[1]),float(S_i[2])]
            solutions = solve_complex_equations(coefficients, constants)
            α = solutions[0]
            β = solutions[1]
            print(point,file=f)
            print('The spinor wave function (α) of this material at CBM or VBM is','%.3f'%α.real,'+','%.3f'%α.imag,'i','(','%.1f'%(100*(α.real**2+α.imag**2)),'%',')',file=f)
            print('The spinor wave function (β) of this material at CBM or VBM is','%.3f'%β.real,'+','%.3f'%β.imag,'i','(','%.1f'%(100*(β.real**2+β.imag**2)),'%',')',file=f)
    f.close()
    return f


if __name__ == "__main__":
    file_path = select_input_filenames()  # Replace with your actual file path
    result = calculate_spin_overlap(file_path)



