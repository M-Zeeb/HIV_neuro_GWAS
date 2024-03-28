##Required functions
source("PCAs_functions.R")

##File annotation
file_annote = "..."

#tat
main_function(genreg = "tat",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "tat",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#pol
main_function(genreg = "pol",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "pol",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#pol nores
main_function(genreg = "pol",res = "_res",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "pol",res = "_res",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#gag
main_function(genreg = "gag",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "gag",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#env
main_function(genreg = "env",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "env",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#rev
main_function(genreg = "rev",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "rev",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#nef
main_function(genreg = "nef",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "nef",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#vif
main_function(genreg = "vif",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "vif",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#vpu 
main_function(genreg = "vpu",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "vpu",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)
#vpr
main_function(genreg = "vpr",res = "",sequencingtype = "_wga", subtype = "_B", file_annote = file_annote)
main_function(genreg = "vpr",res = "",sequencingtype = "_wga", subtype = "_all", file_annote = file_annote)

##Resistance database(partial pol sequences)
main_function(genreg = "pol",res = "",sequencingtype = "_resdb",subtype = "_B", file_annote = file_annote)
main_function(genreg = "pol",res = "_res",sequencingtype = "_resdb",subtype = "_B", file_annote = file_annote)
main_function(genreg = "pol",res = "",sequencingtype = "_resdb",subtype = "_all", file_annote = file_annote)
main_function(genreg = "pol",res = "_res",sequencingtype = "_resdb",subtype = "_all", file_annote = file_annote)
