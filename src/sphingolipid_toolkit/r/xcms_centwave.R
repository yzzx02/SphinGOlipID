# Independent per-file XCMS peak detection. No targets, alignment or filling.
# CentWaveParam units and output: https://sneumann.github.io/xcms/
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(MSnbase))
suppressPackageStartupMessages(library(xcms))
suppressPackageStartupMessages(library(BiocParallel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected --stdin or one JSON configuration path")
connection <- if (args[[1]] == "--stdin") file("stdin", encoding = "UTF-8") else file(args[[1]], encoding = "UTF-8")
cfg <- fromJSON(paste(readLines(connection, warn = FALSE), collapse = "\n"))
close(connection)
options(mc.cores = 1L)
register(SerialParam())
cat("XCMS ", as.character(packageVersion("xcms")), ": ", cfg$input, "\n", sep = "")
raw <- readMSData(files = cfg$input, mode = "onDisk", msLevel. = 1L)
raw <- filterPolarity(raw, polarity. = if (cfg$polarity == "+") 1L else 0L)
raw <- filterEmptySpectra(raw)
if (!length(raw)) stop("No nonempty MS1 spectra with the selected polarity")
centroided <- centroided(raw)
if (anyNA(centroided) || any(!centroided)) stop("centWave requires centroided MS1 spectra; centroid mzML before analysis")

param <- CentWaveParam(
    ppm = cfg$precursor_ppm,
    peakwidth = c(cfg$peakwidth_min, cfg$peakwidth_max),
    noise = cfg$noise,
    snthresh = cfg$snthresh,
    prefilter = c(3L, cfg$noise),
    mzdiff = -0.001,
    mzCenterFun = "wMean",
    integrate = 1L,
    fitgauss = FALSE
)
xdata <- findChromPeaks(raw, param = param, BPPARAM = SerialParam(),
                        return.type = "XCMSnExp", msLevel = 1L)
pks <- chromPeaks(xdata)
columns <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo", "sn")
if (nrow(pks)) {
    peaks <- data.frame(xcms_peak_id = rownames(pks), pks[, columns, drop = FALSE],
                        check.names = FALSE, row.names = NULL)
} else {
    peaks <- data.frame(xcms_peak_id = character())
    for (key in columns) peaks[[key]] <- numeric()
}
write.csv(peaks, cfg$output, row.names = FALSE, fileEncoding = "UTF-8")
write_json(list(R = as.character(getRversion()), xcms = as.character(packageVersion("xcms")),
                MSnbase = as.character(packageVersion("MSnbase")), mzR = as.character(packageVersion("mzR")),
                ms1_scan_count = length(raw), detected_peak_count = nrow(peaks),
                algorithm = "centWave", parallel_workers = 1L,
                centwave_parameters_json = toJSON(as.list(param), auto_unbox = TRUE)),
           cfg$metadata, auto_unbox = TRUE, pretty = TRUE)
cat("XCMS detected ", nrow(peaks), " local chromatographic peaks\n", sep = "")
