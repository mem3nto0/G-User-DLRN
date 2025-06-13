#%%
import sys
import numpy as np
import tensorflow as tf
import os
import matplotlib.pyplot as plt
from print_table import Table_print
from skimage.transform import resize
from draw_model import model_drawing
from plot_results import plotting_results, plotting_results_test
from Save_data import Save_analysis
from Analyse_residuals import Check_best_residualts, Check_best_residualts_GEL
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QComboBox, QPushButton, QFileDialog, QLineEdit
from PyQt5.QtCore import Qt  # Import Qt module from PyQt5
import PyQt5

class MyWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        gpus = tf.config.experimental.list_physical_devices('GPU')

        for gpu in gpus:

            tf.config.experimental.set_memory_growth(gpu, True)

        time_1 = np.arange(0,20,0.5)
        time_2 = np.logspace(1, 2.3, num=(256-len(time_1)), base=20)
        self.time =  np.concatenate((time_1, time_2), axis = 0)

        self.wl = np.arange(0,256,1)

        self.setStyleSheet("background-color: #444444; color: #ffffff;")

        self.setWindowTitle("DLRN APP")
        self.resize(400, 300)  # Width, Height

        # Create the main widget and layout
        main_widget = QWidget(self)
        layout = QVBoxLayout(main_widget)

        # Create the second QComboBox and add it to the layout
        self.combo_box_1 = QComboBox(self)
        self.combo_box_1.addItems(["---","Spectra emission","Spectra TA", "Agarose Gel"])
        layout.addWidget(self.combo_box_1)

        # Create the second QComboBox and add it to the layout
        self.combo_box_2 = QComboBox(self)
        self.combo_box_2.addItems(["---","Top1", "Top3"])
        layout.addWidget(self.combo_box_2)

        # Create the label "Adjusting Factor"
        self.adjusting_factor_label = QLabel("Adjusting Factor (recommended:1 max_val:10):")
        layout.addWidget(self.adjusting_factor_label, alignment=Qt.AlignLeft)

        # Create the QLineEdit widget
        self.adjusting_factor_textbox = QLineEdit()
        self.adjusting_factor_textbox.setText("-")
        self.adjusting_factor_textbox.setFixedWidth(100)
        self.adjusting_factor_textbox.textChanged.connect(self.update_id)
        layout.addWidget(self.adjusting_factor_textbox)

        self.load_Files = QPushButton("Load custom time scale", self)
        self.load_Files.clicked.connect(self.load_time_scale)
        layout.addWidget(self.load_Files)

        self.load_Files = QPushButton("Load data", self)
        self.load_Files.clicked.connect(self.load_data)
        layout.addWidget(self.load_Files)

        self.load_Files = QPushButton("test DLRN", self)
        self.load_Files.clicked.connect(self.test_DLRN)
        layout.addWidget(self.load_Files)

        self.load_Files = QPushButton("Data analysis", self)
        self.load_Files.clicked.connect(self.Analyse_data)
        layout.addWidget(self.load_Files)

        # Add stretch to push the label to the left
        layout.addStretch()

        # Set the main widget and layout
        self.setCentralWidget(main_widget)

        # Connect signals (optional)
        self.combo_box_1.currentIndexChanged.connect(self.on_combobox1_changed)
        self.combo_box_2.currentIndexChanged.connect(self.on_combobox2_changed)

    def update_id(self, text):
        try:

            self.adjusting_factor = float(text)

        except ValueError:
            None
            # Handle non-numeric input gracefully
            pass

    def on_combobox1_changed(self, index):

        selected_option = self.combo_box_1.itemText(index)

        def top_3_acc(y_true, y_pred):
            return tf.keras.metrics.top_k_categorical_accuracy(y_true, y_pred, k=3)

        path = os.getcwd()

        if selected_option == "Spectra emission":

            try:
                path_model = path + "/pretrained_model_Spectra/pretrained_model_Spectra/resV2_spectra_model"
                model_model = tf.keras.models.load_model(path_model, custom_objects={'top_3_acc':top_3_acc})
                self.model_model = model_model

                path_amp = path + "/pretrained_model_Spectra/pretrained_model_Spectra/resV2_spectra_amplitude" # 
                model_amp = tf.keras.models.load_model(path_amp)
                self.model_amp = model_amp

                path_tau = path + "/pretrained_model_Spectra/pretrained_model_Spectra/resV2_spectra_Tau"
                model_tau = tf.keras.models.load_model(path_tau)
                self.model_tau = model_tau

            except:
                path_model = path + "/pretrained_model_Spectra/resV2_spectra_model"
                model_model = tf.keras.models.load_model(path_model, custom_objects={'top_3_acc':top_3_acc})
                self.model_model = model_model

                path_amp = path + "/pretrained_model_Spectra/resV2_spectra_amplitude" # 
                model_amp = tf.keras.models.load_model(path_amp)
                self.model_amp = model_amp

                path_tau = path + "/pretrained_model_Spectra/resV2_spectra_Tau"
                model_tau = tf.keras.models.load_model(path_tau)
                self.model_tau = model_tau

        if selected_option == "Spectra TA":

            try:
                path_model = path + "/pretrained_model_TA/resV2_TA_spectra_Model"
                model_model = tf.keras.models.load_model(path_model, custom_objects={'top_3_acc':top_3_acc})
                self.model_model = model_model

                path_amp = path + "/pretrained_model_TA/resV2_TA_spectra_Amplitude" # 
                model_amp = tf.keras.models.load_model(path_amp)
                self.model_amp = model_amp

                path_tau = path + "/pretrained_model_TA/resV2_TA_spectra_cosh_Tau"
                model_tau = tf.keras.models.load_model(path_tau)
                self.model_tau = model_tau

            except:
                path_model = path + "/pretrained_model_TA/pretrained_model_TA/resV2_TA_spectra_Model"
                model_model = tf.keras.models.load_model(path_model, custom_objects={'top_3_acc':top_3_acc})
                self.model_model = model_model

                path_amp = path + "/pretrained_model_TA/pretrained_model_TA/resV2_TA_spectra_Amplitude" # 
                model_amp = tf.keras.models.load_model(path_amp)
                self.model_amp = model_amp

                path_tau = path + "/pretrained_model_TA/pretrained_model_TA/resV2_TA_spectra_cosh_Tau"
                model_tau = tf.keras.models.load_model(path_tau)
                self.model_tau = model_tau

        if selected_option == "Agarose Gel": # to change for the GEL path

            try:
                path_model = path + "/pretrained_model_Gel/pretrained_model_Gel/resV2_GEL_takemodel"
                model_model = tf.keras.models.load_model(path_model, custom_objects={'top_3_acc':top_3_acc})
                self.model_model = model_model

                path_amp = path + "/pretrained_model_Gel/pretrained_model_Gel/resV2_GEL_takespectra+dilat"
                model_amp = tf.keras.models.load_model(path_amp)
                self.model_amp = model_amp

                path_tau = path + "/pretrained_model_Gel/pretrained_model_Gel/resV2_GEL_Tau(MSE)"
                model_tau = tf.keras.models.load_model(path_tau)
                self.model_tau = model_tau

            except:

                path_model = path + "/pretrained_model_Gel/resV2_GEL_takemodel"
                model_model = tf.keras.models.load_model(path_model, custom_objects={'top_3_acc':top_3_acc})
                self.model_model = model_model

                path_amp = path + "/pretrained_model_Gel/resV2_GEL_takespectra+dilat"
                model_amp = tf.keras.models.load_model(path_amp)
                self.model_amp = model_amp

                path_tau = path + "/pretrained_model_Gel/resV2_GEL_Tau(MSE)"
                model_tau = tf.keras.models.load_model(path_tau)
                self.model_tau = model_tau

        self.Kind_analysis = selected_option

    def on_combobox2_changed(self, index):

        selected_option = self.combo_box_2.itemText(index)

        if selected_option == "Top1":

            self.kind_Top = 1

        if selected_option == "Top3":

            self.kind_Top = 3        

    def load_time_scale(self):

        file_path, _ = QFileDialog.getOpenFileName()
        time_original = np.loadtxt(file_path)

        if self.adjusting_factor == -1:

            self.rescale = 1000/np.max(time_original)

        else:
            self.rescale = 1/self.adjusting_factor

        self.time_original = time_original*self.rescale

    def load_data(self):

        file_path, _ = QFileDialog.getOpenFileName()

        try:

            data_general = np.load(file_path)
            data = data_general[data_general.files[0]]

        except:
            
            data = np.loadtxt(file_path)

        time_data = self.time_original # i need few passage to make the time for training and real data the same

        time_training = self.time # right now, time data and time training are the same. to check how to fix it


        if data.shape[0] != 256:

            data_resized = resize(data, (256,data.shape[1]), anti_aliasing= True)
        
        else:

            data_resized = data
                
        Final_data = np.zeros([256,256])

        for j in range(256):

            Final_data[j,:] = np.interp(time_training, time_data, data_resized[j,:])

        norm_value = np.max(np.abs(Final_data))
        Final_data = Final_data / norm_value  

        self.data_to_analyse = Final_data
        self.norm_value = norm_value
        
        if self.Kind_analysis == "Spectra emission":
            
            plt.figure(1)   
            plt.contourf(time_training, self.wl, Final_data)
            plt.show()

        else:

            plt.figure(1)   
            plt.contourf(Final_data)
            plt.show()

    def test_DLRN(self):

        if self.Kind_analysis == "Spectra emission":

            path = os.getcwd()
            path = path + "/Test_data_APP/Spectra"

        else: 

            path = os.getcwd()
            path = path + "/Test_data_APP/Gel"

        list = os.listdir(path)
        rand_data = np.random.choice(list)

        with np.load(path + "/" + rand_data) as data_loaded: #"/" + rand_data

            data = data_loaded["train"]
            #data = data/np.max(data)
            data = np.expand_dims(data,axis = -1)
            data = np.expand_dims(data,axis = 0)

            true_tau = data_loaded["label1"]
            true_model = data_loaded["label2"]
            true_amp = data_loaded["label3"]

        path = os.getcwd()
        dirs = path + "/Model_information"

        with np.load(dirs + "/onehot-reppresentation + binary-matrix.npz") as data_info:

            Binary_matrix = data_info["binary_matrix"]

        #model_model = Global_variable.model_model
        #model_tau = Global_variable.model_tau
        #model_amp = Global_variable.model_amp
        #Time = Global_variable.time
        #kind_top = Global_variable.kind_Top

        pre_model = self.model_model.predict(data)
        top_3_indices = np.argpartition(pre_model[0,:],-int(self.kind_Top))[-int(self.kind_Top):]

        plt.figure(1)
        plt.title("model prediction one-hot encoding")
        plt.plot(pre_model[0,:],marker = "o", color = "red", label = "prediction")
        plt.plot(true_model, marker = "o", color = "black", label = "expectation")
        plt.ylabel('probability confidence', fontsize=12)
        plt.xlabel('model index', fontsize=12)
        plt.legend()
        plt.show()

        Tau_solutions = np.zeros([len(top_3_indices), 7])

        for i in range (len(top_3_indices)):

            model_binary = Binary_matrix[top_3_indices[i]]
            model_binary = np.expand_dims(model_binary, axis = 0)
            X_total = {"input_1": data,"input_2": model_binary}

            print(model_binary)

            pre_tau = self.model_tau.predict(X_total)
            pre_amp = self.model_amp.predict(X_total)

            pre_tau[0,1:][pre_tau[0,1:] < 0.8] = 0
            Tau_solutions[i,:] = pre_tau[0,:]

            if self.Kind_analysis == "Spectra emission" or self.Kind_analysis == "Spectra TA":

                Fit, score, Kinetic_signal = Check_best_residualts(self.time, data[0,:,:,0], model_binary[0,:,:], pre_tau[0,:], pre_amp[0,:,:])
            
            else:

                Fit, score, Kinetic_signal = Check_best_residualts_GEL(self.time, data[0,:,:,0], model_binary[0,:,:], pre_tau[0,:], pre_amp[0,:,:])

            plotting_results_test(model_binary, 
                            pre_model, 
                            i, 
                            top_3_indices, 
                            pre_amp, 
                            self.time, 
                            Kinetic_signal, data, Fit, score,
                            true_amp, self.Kind_analysis)

        Table_print(Tau_solutions, true_tau, self.kind_Top)


    def Analyse_data(self):

        data = self.data_to_analyse
        data = np.expand_dims(data, axis= -1)
        data = np.expand_dims(data, axis= 0)

        path = os.getcwd()
        dirs = path + "/Model_information"

        with np.load(dirs + "/onehot-reppresentation + binary-matrix.npz") as data_info:

            Binary_matrix = data_info["binary_matrix"]

        pre_model = self.model_model.predict(data)
        top_3_indices = np.argpartition(pre_model[0,:],-int(self.kind_Top))[-int(self.kind_Top):]

        plt.figure(1)
        plt.title("model prediction one-hot encoding")
        plt.plot(pre_model[0,:],marker = "o", color = "red", label = "prediction")
        plt.ylabel('probability confidence', fontsize=12)
        plt.xlabel('model index', fontsize=12)
        plt.legend()
        plt.show()

        Tau_solutions = np.zeros([len(top_3_indices), 7])

        # Create a folder to save the analysis
        folder_name = 'Analysis_solutions'
        os.makedirs(folder_name, exist_ok=True)

        for i in range (len(top_3_indices)):

            model_binary = Binary_matrix[top_3_indices[i]]

            N_taus = 0

            for kk in range (5):

                N_taus = N_taus + np.abs(model_binary[kk,kk])

            N_taus =int(N_taus)

            model_binary = np.expand_dims(model_binary, axis = 0)
            X_total = {"input_1": data,"input_2": model_binary}

            pre_tau = self.model_tau.predict(X_total)
            pre_amp = self.model_amp.predict(X_total)

            pre_tau[0,N_taus:] = 0

            Tau_solutions[i,:] = pre_tau[0,:]/self.rescale

            if self.Kind_analysis == "Spectra emission" or self.Kind_analysis == "Spectra TA":

                Fit, score, Kinetic_signal = Check_best_residualts(self.time, data[0,:,:,0], model_binary[0,:,:], pre_tau[0,:], pre_amp[0,:,:])
            
            else:

                Fit, score, Kinetic_signal = Check_best_residualts_GEL(self.time, data[0,:,:,0], model_binary[0,:,:], pre_tau[0,:], pre_amp[0,:,:])

            plotting_results(model_binary, 
                            pre_model, 
                            i, 
                            top_3_indices, 
                            pre_amp, 
                            self.time, 
                            Kinetic_signal, data, Fit, score, self.Kind_analysis)


            residuals = data[0,:,:,0] - Fit[:,:]

            Inner_folder = path +"/"+ folder_name
            Save_analysis(i, Inner_folder, pre_amp[0,:,:], Tau_solutions[i,:], Kinetic_signal, residuals)

        true_tau = np.zeros(7)
        Table_print(Tau_solutions, true_tau, self.kind_Top)
 

if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
