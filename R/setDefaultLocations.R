uniprobeDirectory = "/p/keles/CAGI2015/volumeB/Data/UniProbe"
chipLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/ChIP-seq/metadata.tsv"
dnaseLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/DNase-seq/metadata.tsv"
dnamethLocation = "/p/keles/CAGI2015/volumeB/Data/ENCODE/List/DNAmethylation/metadata.tsv"
peakDir = "/p/keles/CAGI2015/volumeB/Data/ENCODE/peaks"
inputDTLocation = "/p/keles/CAGI2015/volumeB/ProcessedData/SampleDT.RData"
outputDir = "/p/keles/CAGI2015/volumeB/ProcessedData"
# saving to package
save(uniprobeDirectory, 
     chipLocation,
     dnaseLocation,
     dnamethLocation,
     peakDir,
     inputDTLocation,
     outputDir,
     file = "~/keles/eqtl.pipeline/data/defaultFileLocations.Rda")


