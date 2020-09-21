# To prevent roxygen errors because of dplyr stuff (things that look like
# variables but aren't)
utils::globalVariables(c("dimension", 
	                     "perc",
	                     "variable",
	                     "category",
	                     "Freq",
	                     "Freq.y",
	                     "Freq.x"))