args <- commandArgs(trailingOnly=TRUE)
env <- new.env()
load(args[1], envir=env)
d <- env$count.data.frame
if (!"Name" %in% names(d)) {
  inf <- vapply(env$count.list, function(x) x$inf, character(1))
  stopifnot(all(grepl("-Fomula-",inf,fixed=TRUE)))
  d <- data.frame(
    Name=sub("^DB#: ","",sub("-Fomula-.*$","",inf)),
    Formula=sub("-CCS-.*$","",sub("^.*-Fomula-","",inf)),
    Adduct=sub("^.*-Add-","",inf)
  )
  rm(inf)
}
rm(env)
gc()
cat("DATA FRAME COLUMNS\n")
print(names(d))
print(head(d,2))
saveRDS(d, paste0(args[2],"_frame.rds"))
find_column <- function(options) {
  found <- names(d)[tolower(names(d)) %in% tolower(options)]
  if (length(found)!=1) stop(paste("Cannot resolve column",paste(options,collapse=",")))
  found
}
names(d)[names(d)==find_column(c("Name","inf","name2"))] <- "Name"
names(d)[names(d)==find_column(c("Formula","formula2"))] <- "Formula"
names(d)[names(d)==find_column(c("Adduct","add"))] <- "Adduct"
# No EQ, LCI, spectrum matching or sample files are used.
des <- sub("^LEVEL[0-9]+_DES_", "", d$Name)
heads <- sub(" .*", "", des)
counts <- sort(table(heads), decreasing=TRUE)
print(counts)
sp <- grepl("Cer|^SM$|^ASM$|^SPB$|^SL($|[+])|^GD|^GM|^GQ|^GT|cysteinolide|SAL", heads, ignore.case=TRUE)
selected <- d[sp,c("Name", "Formula", "Adduct")]
cat("SP NAME EXAMPLES\n")
print(head(selected,8))
# Keep both precursor description and chain description, remove only LEVEL tag.
selected$structure_key <- sub("^LEVEL[0-9]+_DES_", "", selected$Name)
unique_structures <- unique(selected[c("structure_key", "Formula")])
write.table(unique_structures, paste0(args[2],"_structures.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
write.table(data.frame(class=names(counts),records=as.integer(counts)),paste0(args[2],"_classes.tsv"),sep="\t",row.names=FALSE,quote=TRUE)
cat("SP_RECORDS",nrow(selected),"SP_UNIQUE_STRUCTURES",nrow(unique_structures),"\n")
