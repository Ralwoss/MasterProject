library(rGMAP)
path <- "/home/ralwoss/Uni/MasterProject/data/intraChromosomalMatrices"
# files <- list.files(path=path)
# print(files)
# for (f in files) {
#   if(endsWith(f,".txt") == FALSE) {
#     next
#   }
#   print(f)
#   t = read.table(paste(path, f, sep="/"))
#   res = rGMAP(t, resl = 10*1000)
#   remove(t)
#   print(res)
#   write.table(res$tads, paste("/home/ralwoss/Uni/MasterProject/data/rGMAPTADs", f, sep = "/"), sep="\t")
#   res <- NULL
# }
f = "chrY.txt"

t = read.table(paste(path, f, sep="/"))
res = rGMAP(t, resl = 10*1000)
write.table(res$tads, paste("/home/ralwoss/Uni/MasterProject/data/rGMAPTADs", f, sep = "/"), sep="\t")