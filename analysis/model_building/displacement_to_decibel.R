# displacement_to_decibel.R
#
# Conversion between tactor gain (device units, 0-255) and displacement of
# the C-2 tactor at 250 Hz, based on a linear fit to the manufacturer's
# measurements ("Displacement Testing C-2.xlsx", Engineering Acoustics, Inc.).
#
#   gain_to_displacement()    : gain -> displacement (mm)
#   displacement_to_gain()    : displacement (mm) -> gain
#   displacement_to_decibel() : displacement -> dB, 40 * log10(d / d_ref)
#                               with d_ref the displacement at gain 4
#
# The fit for gains 1-255 is stored in
# ../matching_gain_to_displacement_250_Hz.txt (used by gain_to_decibel.R).
# The reference d_ref = 0.001578822 mm and the 40 * log10 definition are
# hard-coded in db_disp() / db_inv_disp() of gpm_generic.stan.
#
# Not sourced by any other script.

gain_to_displacement <- function(gain_val){
  # read matching table 
  # (see "Displacement Testing C-2.xlsx" by the manufacturer; 
  # C-2 @ 250Hz section)
  gain_to_disp <- data.frame(gain = c(35,53,75,115,150,200,250,255),
                             displacement = c(0.0612, 0.0856, 0.13, 0.198, 0.266, 0.356, 0.444, 0.452))  
  # Predict linear relation 
  model_disp_gain <- lm(displacement ~ gain, gain_to_disp)
  # Get predictions from linear model and return
  disp_val <- model_disp_gain$coefficients[1] + model_disp_gain$coefficients[2]*gain_val
  return(disp_val)
}

displacement_to_gain <- function(disp_val){
  # read matching table 
  # (see "Displacement Testing C-2.xlsx" by the manufacturer; 
  # C-2 @ 250Hz section)
  gain_to_disp <- data.frame(gain = c(35,53,75,115,150,200,250,255),
                             displacement = c(0.0612, 0.0856, 0.13, 0.198, 0.266, 0.356, 0.444, 0.452))  
  # Predict linear relation 
  model_disp_gain <- lm(gain ~ displacement, gain_to_disp)
  # Get predictions from linear model and return
  gain_val <- model_disp_gain$coefficients[1] + model_disp_gain$coefficients[2]*disp_val
  return(gain_val)
}

displacement_to_decibel <- function(disp_val){
  library(dplyr)
  
  # Get lowest positively predicted displacement for the tactors as reference
  ref <- gain_to_displacement(4)

  # Create a function to process one value
  process_one <- function(d_val) {
    # Calculate db 
    db_val <- 40 * log10(d_val/ref)
    return(db_val)
  }
  
  # Apply to each element in disp_val
  return(sapply(disp_val, process_one))
}

# Testing function for whole range -------------------------------------------
pred <- gain_to_displacement(4:255)
db_vib <- displacement_to_decibel(pred)
# plot(pred, db_vib)
