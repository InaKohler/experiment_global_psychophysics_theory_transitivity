
gain_to_decibel <- function(gain_val){
  library(dplyr)
  
  # read matching table
  if(!exists("matching_gain_displ")){
    # matching_gain_displ <- read.table(here::here("utils", "matching_gain_to_displacement_250_Hz.txt"))
    matching_gain_displ <- read.table("matching_gain_to_displacement_250_Hz.txt")
  } 
  # Create a function to process one value
  process_one <- function(g_val) {
    # Match displacement from gain
    displ_val <- matching_gain_displ %>%
      filter(gain == g_val) %>%
      pull(displacement)
    
    # Reference value
    ref <- matching_gain_displ %>%
      filter(gain == 9) %>%
      pull(displacement)
    
    # Calculate db 
    db_val <- 50 * log10(displ_val/ref)
    return(db_val)
  }
  
  # Apply to each element in gain_val
  return(sapply(gain_val, process_one))
}
