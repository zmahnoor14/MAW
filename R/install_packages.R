
install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'uuid', 'digest')) #repos="https://cloud.r-project.org"))

install.packages(c("rvest", "xml2", "dplyr", "stringr", "readr", "remotes", "ncdf4", 
                   "parallel", "doParallel", "foreach", "future", "listenv", "curl", "provViz",
		   "rjson", "jsonlite"))
#install.packages("devtools", repos="http://cloud.r-project.org")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install(c("Spectra", "MsCoreUtils", "mzR"))
BiocManager::install(c("CompoundDb"))
BiocManager::install("CAMERA")
remotes::install_github("rformassspectrometry/MsBackendHmdb")
remotes::install_github("rformassspectrometry/MsBackendMsp")
remotes::install_github("rformassspectrometry/MsBackendMgf")
Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
install.packages("V8", repos="https://cloud.r-project.org")
#devtools::install_github("End-to-end-provenance/provSummarizeR")
#devtools::install_github("End-to-end-provenance/rdt")
#devtools::install_github("End-to-end-provenance/rdtLite")
#devtools::install_github("mmondelli/rdt2repr")
