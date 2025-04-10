#%%
import numpy as np
from scipy.integrate import odeint


"""
//////////////////////////////////////////////////
/// Create kinetic traces from a kinetic model ////
//////////////////////////////////////////////////

"""

def Create_kinetc_signal(matrix_with_k, time, t_0, time_before_zero):

    def Resolve_ODE(z, time, matrix_with_k, t_0):

        dzdt = np.zeros(matrix_with_k.shape[0]) 

        irf = np.exp(-2*((time-t_0)/0.1)**2)

        for i in range (matrix_with_k.shape[0]):

            if i == 0:              
                dzdt[i] = np.matmul(matrix_with_k[i,:],z.T) # + irf
            else:
                dzdt[i] = np.matmul(matrix_with_k[i,:],z.T)
        return dzdt

    if time_before_zero == True:

    # initial condiction
        Initial = np.zeros(matrix_with_k.shape[0])


        # resolve differential equation
        z = np.zeros([len(time),matrix_with_k.shape[0]])

        for i in range(1,len(time)): #this for loop is usefull to stablish max(z) = z(t_0)
        
            t_span = [time[i-1],time[i]]
            z_step = odeint(Resolve_ODE, Initial, t_span, args=(matrix_with_k,t_0))

            z[i,:] = z_step[1]
            Initial = z_step[1]

    else:

        Initial = np.zeros(matrix_with_k.shape[0])
        Initial[0] = 1
        
        z = odeint(Resolve_ODE, Initial, time, args=(matrix_with_k,t_0))
        
    return z

"""
////////////////////////////////////////////////////////
/// associate to each pathway the kinetic rate constant, ////
/// generating the kinetic matrix with rate constants    ////
/////////////////////////////////////////////////////////

"""

def Generate_matrix_withrate(N_species, Model_pred, k_values):

    matrix_with_k = np.zeros([N_species,N_species])

    k_use = k_values.copy()

    for i in range (N_species):

        for j in range (N_species):
            
            probe = Model_pred[j,i]

            if probe <= -1:

                for n in range(int(abs(probe))):

                    matrix_with_k[j,i]  = matrix_with_k[j,i]  - k_use[n]

            if probe >= 1:

                for n in range(int(abs(probe))):

                    matrix_with_k[j,i]  = matrix_with_k[j,i]  + k_use[n]
                    k_use = np.delete(k_use, 0) 

    return matrix_with_k
