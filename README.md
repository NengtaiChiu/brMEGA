
# Citation
  Please include the standard citation to brMEGA when using the resources available on this platform:<br/>


  Neng-Tai Chiu, Stephanie Huwiler, M. Laura Ferster, Walter Karlen, Hau-Tieng Wu*, Caroline Lustenberger* (2021). Get rid of the beat in mobile EEG applications: A framework towards automated cardiogenic artifact detection and removal in single-channel EEG. Biomedical Signal Processing and Control: https://doi.org/10.1016/j.bspc.2021.103220<br/>
  
   Note: * Contributed equally to work.<br/>
# Introduction

This is a package for the algorithm brMEGA, a method to remove cardiogenic artifacts in EEG signal.<br/>
the data is in the link inside the folder "database"(you need to download it from the provided link in dropbox).<br/>
the labels are included in the folder "data_label".<br/>

# stage 1
  "step_1_algo_training.m": training the SVM model for detecting cardiogenic artifacts in EEG signal.<br/>
  "step_1_algo_prediction.m": predicting Epochs with/without cardiogenic artifacts in EEG signal.<br/>
  
# stage 2: removing the cardiogenic artifact
 ## EAS method
   "Step2_algo_EAS_no_prediction.m" : using the EAS method to remove cardiogenic artifact on the whole EEG signal.<br/>
   "Step2_algo_EAS_with_prediction.m" : using the EAS method to remove cardiogenic artifact only on  epochs that are predicted containing cardiogenic artifacts.<br/>
 ## NLEM method
   "Step2_algo_brMEGA_no_prediction.m" : using the brMEGA method to remove cardiogenic artifact on the whole EEG signal.<br/>
   "Step2_algo_brMEGA_with_prediction.m" : using the brMEGA method to remove cardiogenic artifact only on  epochs that are predicted containing cardiogenic artifacts.<br/>
# stage 3:  evaluation 
  ## evaluation code: 
  "aR_sC_extraction_combine.m": evaluation code for comparing EAS and brMEGA.<br/>
  use "Step3_evaluation_brMEGA_with_prediction_for_aR_sC.m" and "Step3_evaluation_EAS_with_prediction_for_aR_sC.m" to create dataset for evaluation.

  
  
