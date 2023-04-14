# provenance library rdtLite integration
library(rdtLite)
options(prov.dir = "./prov_metfrag", snapshot.size = 10000)
prov.init(prov.dir = ".")

#define arguments
args = commandArgs(trailingOnly = TRUE)

# Start time
start.time <- Sys.time()

met_param <- args[1]
MetFragjarFile <- args[2]

print(met_param)
print(MetFragjarFile)

run_metfrag <- function(met_param, MetFragjarFile){
    filesmet_param <- read.table(met_param)
    for (files in filesmet_param[[1]]){
        print(paste("java -jar",  MetFragjarFile, files))
        system(paste("java -jar",  MetFragjarFile, files))
        Sys.sleep(5)
    }
}

run_metfrag(met_param, MetFragjarFile)

end.time <- Sys.time()

# Time taken to run the analysis for MAW-R
time.taken <- end.time - start.time
print(time.taken)

prov.save()
prov.quit()