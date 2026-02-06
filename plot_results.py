import matplotlib.pyplot as plt
from draw_model import model_drawing
import numpy as np

def plotting_results(model_binary, pre_model, ind, top_3_indices, pre_amp, Time, Kinetic_signal, data, Fit, score):

    wl = np.arange(0,256,1)

    fig, axs = plt.subplots(1, 4,figsize=(10, 3), dpi=150, num="Solution {}".format(ind+1))

    #plot the model prediction
    model_drawing(axs[0], model_binary, pre_model[0,top_3_indices[ind]])  # Replace prob with your actual probability

    # plot the amplitude predictions

    colors = ["red","purple","green","orange"]
    names1 = ["A","B","C","D"]

    for j in range (4):

        axs[1].plot(pre_amp[0,:,j], linewidth=2, color= colors[j], label = names1[j]) #

    axs[1].set_title("amplitude prediction", fontsize=12)
    axs[1].legend(ncol=2, loc='upper right', fontsize=8)            
    axs[1].set_xlabel('wl (unit)', fontsize=10)  # Set xlabel for the subplot
    axs[1].set_ylabel('norm. I', fontsize=10)  # Set ylabel for the subplot
    axs[1].set_xlim(0, 256)
    axs[1].set_ylim(0, 1.1)

    # plot the population profile

    names3 = ["A","B","C","D"]

    for j in range (4):

        axs[2].plot(Time, Kinetic_signal[:,j], linewidth=2, color= colors[j], label = names3[j]) #

    axs[2].set_title("population profile", fontsize=12)
    axs[2].legend(ncol=2, loc='upper right', fontsize=8)            
    axs[2].set_xlabel('time (unit)', fontsize=10)  # Set xlabel for the subplot
    axs[2].set_ylabel('norm. pop.', fontsize=10)
    axs[2].set_xlim(0, 1000)
    axs[2].set_ylim(0, 1.1)

    # plot the residuals

    levels = np.arange(-0.3, 0.31, 0.02)
    res_2D = axs[3].contourf(Time,wl,(data[0,:,:,0] - Fit[:,:]), levels = levels, cmap ="rainbow")
    axs[3].set_xlabel('time (unit)', fontsize=10)
    axs[3].set_ylabel('wl (unit)', fontsize=10)
    axs[3].set_title("residual = {}".format(np.trunc(score)), fontsize=10)
    fig.colorbar(res_2D, ax=axs[3], ticks=[-0.3, 0, 0.3])

    #plt.subplots_adjust(wspace=0.3, hspace=0.4)  # You can adjust these values as needed
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.show(block=False)

    plt.show()

def plotting_results_test(model_binary, pre_model, ind, top_3_indices, pre_amp, Time, Kinetic_signal, data, Fit, score, true_amp):

    wl = np.arange(0,256,1)

    fig, axs = plt.subplots(1, 4,figsize=(10, 3), dpi=150, num="Solution {}".format(ind+1))

    #plot the model prediction
    model_drawing(axs[0], model_binary, pre_model[0,top_3_indices[ind]])  # Replace prob with your actual probability

    # plot the amplitude predictions

    colors = ["red","purple","green","orange"]
    names1 = ["Pr. A","Pr. B","Pr. C","Pr. D"]
    names2 = ["Ex. A","Ex. B","Ex. C","Ex. D"]

    for j in range (4):

        axs[1].plot(pre_amp[0,:,j], linewidth=2, color= colors[j], label = names1[j]) #
        axs[1].plot(true_amp[:,j],"--", linewidth=2, color= colors[j], label = names2[j]) # 

    axs[1].set_title("amplitude prediction", fontsize=12)
    axs[1].legend(ncol=2, loc='upper right', fontsize=8)            
    axs[1].set_xlabel('wl (index)', fontsize=10)  # Set xlabel for the subplot
    axs[1].set_ylabel('norm. I', fontsize=10)  # Set ylabel for the subplot
    axs[1].set_xlim(0, 256)
    axs[1].set_ylim(0, 1.1)

    # plot the population profile

    names3 = ["A","B","C","D"]

    for j in range (4):

        axs[2].plot(Time, Kinetic_signal[:,j], linewidth=2, color= colors[j], label = names3[j]) #

    axs[2].set_title("population profile", fontsize=12)
    axs[2].legend(ncol=2, loc='upper right', fontsize=8)            
    axs[2].set_xlabel('time (unit)', fontsize=10)  # Set xlabel for the subplot
    axs[2].set_ylabel('norm. pop.', fontsize=10)
    axs[2].set_xlim(0, 1000)
    axs[2].set_ylim(0, 1.1)

    # plot the residuals

    levels = np.arange(-0.3, 0.31, 0.02)
    res_2D = axs[3].contourf(Time,wl,(data[0,:,:,0] - Fit[:,:]), levels = levels, cmap ="rainbow")
    axs[3].set_xlabel('time (unit)', fontsize=10)
    axs[3].set_ylabel('wavelength (index)', fontsize=10)
    axs[3].set_title("residual = {}".format(np.trunc(score)), fontsize=10)
    fig.colorbar(res_2D, ax=axs[3], ticks=[-0.3, 0, 0.3])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.show(block=False)

    plt.show()

