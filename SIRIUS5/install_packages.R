
install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'uuid', 'digest'))
install.packages(c("rvest", "xml2", "dplyr", "stringr", "readr", "remotes", "ncdf4", 
                   "parallel", "doParallel", "foreach", "future", "listenv", "curl", "provViz"))
Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
install.packages("V8", repos="https://cloud.r-project.org")