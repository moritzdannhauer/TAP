# fMRI-based_TMS_Targeting_Pipeline
Custom MATLAB code to determine an optimal scalp placement of a TMS coil on the computational model of the subject's head. (1) Optimal coil placement is chosen based on maximizing the TMS-induced electric field in a brain region of interest which is derived from individual fMRI peak activity and computed prospectively for neuronavigation-supported TMS experiments. (2) The code's functionality allows also for accuracy assessment of the coil placement recorded during and analyzed after a TMS session.  
  
Dependencies:  
- SimNIBS 3.2.X  
- 'Tools for NIfTI and ANALYZE image' for loading structural MRI data sets (nii files) version 1.27.0.0 (from MATHWORKS MATLAB Central)  
  
For SimNIBS version < 3.2.5 users, the option to control the smoothness of the scalp surface (scalp_normals_smoothing_steps) to allow changing the scalp tangential plane from MATLAB is not available. However, there is a simple patch you can do yourself and described here: https://github.com/simnibs/simnibs/issues/42
and change smooth=1 to smooth=20 for example.
