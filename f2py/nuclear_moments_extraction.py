#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 15:52:30 2020

@author: asimina
"""
import numpy as np
from nuclear_moments import *

#print(matrices.__doc__)
#print(solver.__doc__)

#A = 237.0
A1 = int(input('Give the mass number of the reference isotope: '))  #238
A2 = int(input('Give the mass number of the target isotope: '))   #236
A = (A1+A2)/2

#lithium-like example
#F = np.array([[-1.849792025703267*(10**5),2.425446341436609*(10**2),-0.635925306138426,\
#               0.001037472142435],[-7.2740620581839*(10**3),9.0085521902942,-0.023664207999006,\
#                                   0.000038259925557]], 'd', order='F')
                                
#dnu = np.array([27422.148184519512, 1084.9898508226213], 'd', order='F')   

#beryllium-like example
F = np.array([[-1.573512465813949*(10**5),2.075124708925421*(10**2),-0.545067800809402,\
               0.000896047878353],[2.260087431698599*(10**5),-2.939843311017920*(10**2),\
                                   0.772387440181695,-0.001266437615842]], 'd', order='F')
                                   
dnu = np.array([-23407.79512057857, 33676.28639191137], 'd', order='F')
#dnu = np.array([200370.1866062507,-288194.0873848016], 'd', order='F')
                                   
dr2 = 0.1638
dr4 = 13.7693
#dr2 = -1.4005
#dr4 = -115.6626

#er_order = 0.001     
er_order = float(input('Give the order of magnitude of the errors in the pseudo - field shifts: '))              

concat_F, C, D = matrices.coefficient_matrices(F,A)

r = solver.matrix_equation(concat_F,dnu)

sigma_f = solver.error_estimation(dnu,concat_F,er_order)

y = solver.matrix_equation(C, dnu)

r_new = solver.matrix_equation(D,y)

sigma_f_new = solver.error_estimation2(dnu,C,D,er_order)

moder_dr2 = r[0] - dr2
moder_dr4 = r[1] - dr4

moder_dr2_new = r_new[0] - dr2
moder_dr4_new = r_new[1] - dr4


print('\n')
print('-'*77)
print(' original  summation: <dr^2> = {:7.4f} (±{:6.4f}), <dr^4> = {:8.4f} (±{:6.4f})'\
      .format(r[0],sigma_f[0,0],r[1],sigma_f[1,1]))
print('-'*77)
print('     model  error   :          {:7.4f},                   {:9.4f}'.format(moder_dr2,moder_dr4))
print('\n')
print('-'*77)
print('rearranged summation: <dr^2> = {:7.4f} (±{:6.4f}), <dr^4> = {:8.4f} (±{:6.4f})'\
      .format(r_new[0],sigma_f_new[0,0],r_new[1],sigma_f_new[1,1]))
print('-'*77)
print('     model  error   :          {:7.4f},                   {:9.4f}'.format(moder_dr2_new,moder_dr4_new))
print('-'*77)




