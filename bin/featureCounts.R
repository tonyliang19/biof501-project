library(Rsubread)

# first parse cli
cli <- list()
bam <- cli$bam
gtf <- cli$gtf


featureCounts(files = bam, annot.ext = gtf, isGTAnnotationFile = TRUE, isPairedEnd = TRUE,
GTF.featureType = "gene")