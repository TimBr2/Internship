#BiocManager::install("profvis")
# Load library
library(profvis)

# Ensure the path is correct and exists
profile_output_path <- "/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/Deploy_ShinyApp_NoAlk/ProfvisOutput.Rprof"

# Run your Shiny app with profiling
p <- profvis({
  runApp('/kyukon/data/gent/vo/000/gvo00027/PPOL/SharedData/2024_TimBruggeman/Deploy_ShinyApp_NoAlk/BruggemanTim_ShinyApp_NoALK.R')
}, prof_output = profile_output_path)
htmlwidgets::saveWidget(p, "profile_after.html") 
