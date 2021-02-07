this is a package for the algorithm brMEGA, a method to remove cardiogenic artifacts in EEG signal.
the data is in the link inside the folder "database"(you need to download it from the provided link in dropbox).
the labels are included in the folder "data_label".
# stage 1
  "step_1_algo_training.m": training the SVM model for detecting cardiogenic artifacts in EEG signal.<br/>
  "step_1_algo_prediction.m": predicting Epochs with/without cardiogenic artifacts in EEG signal.<br/>
  
# stage 2: removing the cardiogenic artifact
 # EAS method
   "Step2_algo_EAS_no_prediction.m" : using the EAS method to remove cardiogenic artifact on the whole EEG signal.<br/>
   "Step2_algo_EAS_with_prediction.m" : using the EAS method to remove cardiogenic artifact only on  epochs that are predicted containing cardiogenic artifacts.<br/>
 # NLEM method
   "Step2_algo_brMEGA_no_prediction.m" : using the brMEGA method to remove cardiogenic artifact on the whole EEG signal.<br/>
   "Step2_algo_brMEGA_with_prediction.m" : using the brMEGA method to remove cardiogenic artifact only on  epochs that are predicted containing cardiogenic artifacts.<br/>
# stage 3:  
  
  
