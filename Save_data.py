import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

"""every data is saved always for A,B,C,B and the time constants for tau1,tau2,tau3..."""

def Save_analysis(ind,path, pre_amp, pre_tau, Kinetic_signal, residuals):

    inner_folder = os.path.join(path, f'Solution_{ind}')
    os.makedirs(inner_folder, exist_ok=True) 

    amp_path = os.path.join(inner_folder, 'Amplitude_solution.csv')
    df_amplitudes = pd.DataFrame(pre_amp, columns=['amp. A', 'amp. B', 'amp. C', 'amp. D'])
    df_amplitudes.to_csv(amp_path, sep='\t', index=False)

    tau_path = os.path.join(inner_folder, 'Tau_solution.csv')
    df_taus = pd.DataFrame([pre_tau], columns=[f'tau{i+1}' for i in range(7)])
    df_taus.to_csv(tau_path, sep='\t', index=False)

    kinetic_path = os.path.join(inner_folder, 'Kinetic_solution.csv')
    df_kinetic = pd.DataFrame(Kinetic_signal[:,:-1], columns=['trace A', 'trace B', 'trace C', 'trace D'])
    df_kinetic.to_csv(kinetic_path, sep='\t', index=False)

    kinetic_path = os.path.join(inner_folder, 'Residuals.txt')
    np.savetxt(kinetic_path, residuals)
    
    plot_file_path = os.path.join(inner_folder, 'Final_analysis_graph.png')
    plt.savefig(plot_file_path)
