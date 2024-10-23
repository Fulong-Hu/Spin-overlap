#   Fulong Hu
#First version : 2024/10/2
#e-mail: 1500148414@qq.com

import numpy as np
import cmath

def solve_complex_equations(coefficients, constants):
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

###########################################Solve for the wave function value at CBM##############################################
constants = [0.008, 0.238, -0.245, 1]   #需要改变的数值 CBM K

solutions = solve_complex_equations(coefficients, constants)
β =solutions[1]
α = solutions[0]
print('The spinor wave function (α) of this material at CBM is','%.3f'%α.real,'+','%.3f'%α.imag,'i','(','%.1f'%(100*(α.real**2+α.imag**2)),'%',')')
print('The spinor wave function (β) of this material at CBM  is','%.3f'%β.real,'+','%.3f'%β.imag,'i','(','%.1f'%(100*(β.real**2+β.imag**2)),'%',')')

###########################################Solve for the wave function value at VBM##############################################
constants = [0.063, 0.322, -0.524, 1]   #需要改变的数值 VBM K
solutions = solve_complex_equations(coefficients, constants)
β1 =solutions[1]
α1 = solutions[0]
print()
print('The spinor wave function (α) of this material at VBM is','%.3f'%α1.real,'+','%.3f'%α1.imag,'i','(','%.1f'%(100*(α1.real**2++α1.imag**2)),'%',')')
print('The spinor wave function (β) of this material at VBM is','%.3f'%β1.real,'+','%.3f'%β1.imag,'i','(','%.1f'%(100*(β1.real**2+β1.imag**2)),'%',')')

###########################################Calculating the value of spin overlap##############################################
β1_conj = β1.conjugate()
α1_conj = α1.conjugate()

results = α1_conj * α + β1_conj * β
spin_overlap = (results.real**2 + results.imag**2)*100
print('The spin overlap of the material is equal to：', '%.2f'%spin_overlap,'%')

