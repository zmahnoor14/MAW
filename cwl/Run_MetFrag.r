# provenance library rdtLite integration
library(rdtLite)
options(prov.dir = "./prov_metfrag", snapshot.size = 10000)
prov.init(prov.dir = ".")

met_param <- args[1]
MetFragjarFile <- args[2]

run_metfrag <- function(met_param, MetFragjarFile){
    start_time <- Sys.time()
    filesmet_param <- read.table(met_param)
    for (files in filesmet_param[[1]]){
        system(paste("java -jar",  MetFragjarFile, files))
        Sys.sleep(5)
    }
}