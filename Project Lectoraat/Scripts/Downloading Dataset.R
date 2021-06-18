fn <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152439&format=file"

het bestand is gedownload als tar.gz
download.file(fn,destfile="GSE152439.gz")

untar gebruikt om de GSM bestanden te verkrijgen in R
untar("GSE152439.gz",list=TRUE)
untar("GSE152439.gz")

ook heb ik het bestand ge unzipped. dit is dus een bestand met alle data.
gunzip("GSE152439.gz")
