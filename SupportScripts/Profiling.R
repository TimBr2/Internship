#BiocManager::install("profvis")
# Load library
library(profvis)

# Ensure the path is correct and exists
profile_output_path <- "./ShinyApp/ProfvisOutput.Rprof"

# Run your Shiny app with profiling
p <- profvis({
  runApp('./ShinyApp/BruggemanTim_ShinyApp_NoALK.R')
}, prof_output = profile_output_path)
htmlwidgets::saveWidget(p, "profile_after.html") 
