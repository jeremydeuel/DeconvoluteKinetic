# Deconvolute scanning kinetics
# GPL 3.0 
# by Jeremy Deuel <jeremy@deuel.ch>
# 21.01.2020
# Dependencies: NNLS package (autoinstalled by this script if necessary)
# Description: This R script deconvolutes a scanning kinetics experiment to the individual spectra.
# Input: Reference Spectra, scanning kinetics experiment (Cary 60 .csv file)
# Output: Deconvoluted .csv files with concentration of substances over time.

# Settings
# set working directory
setwd(getSrcDirectory()[1])
# define the variable referenceSpectraFolder with the path to the folder containing the reference spectra. These spectra have to be saved in files ending on .csv and containing on ly the characters A-Za-z0-9_ in the file name. All spectra in this folder will be used for the deconvolution
referenceSpectraFolder <- "reference spectra"
# define the variable kineticsFile is set to the path with the csv file containing the scanning kinetics. This MUST be an export csv from a Cary 60 scanning kinetics experiments. For other file formats, the importer has to be adapted.
kineticsFile <- "experiments/20151221 transfer.csv"
# define the variable outputFileStem to the path where the output csv file containing the deconvoluted scanning kinetics will be output. The file name will be appended by the sample number and .csv
outputFileStem <- "experiments/example.out."
# check if nnls package is installed, if not, do so.
if (!("nnls" %in% rownames(installed.packages())))
	install.packages("nnls")

# -- start of script, you do not have to modify anything after this line in order to get this script working.

# Function definition: read single spectrum (used for reference spectra)
read.single.spectrum <- function(filename) {
    cat("Reading File ",filename,"\n")
   d <- scan(filename,"",sep="\n")
    x <- d[2]
    m <- regmatches(x,regexec("Conc: (\\d+\\.?\\d*)",x))
    if (length(m)&length(m[[1]][2])) {
        conc <- as.numeric(m[[1]][2])
        if (is.na(conc)) conc <- 1
        else cat("  --> Found Reference Spectrum with Concentration",conc,"µM\n")
    }
    else conc <- 1
    
    d <- na.omit(as.numeric(unlist(lapply(d, function(x) {
        m <- regmatches(x,regexec("^(\\d+\\.?\\d*),(\\d+\\.?\\d*),?",x))
        if (length(m)) return(m[[1]][2:3])
        else return(c(NA,NA))
    }))))
    if (length(d)<2) {
        print("!!! Length of file after match < 0")
        return(NA)
    }
    d <- data.frame(wl=round(d[(1:(length(d)/2))*2-1]),abs=d[(1:(length(d)/2))*2]/conc)
    d <- na.omit(d)
    select <- c()
    wl_done <- c()
    for (wl in d$wl) {
        if (wl%in%wl_done) next
        select <- c(select,which(d$wl == wl)[1])
        wl_done <- c(wl, wl_done)
    }
    d <- d[select, ]
    d <- d[order(d$wl),]
    cat("  -->",dim(d)[1]," points, ",min(d$wl),'-',max(d$wl),'nm with max abs=',max(d$abs),"\n")
    return(d)
}


# Gather the reference spectra
references <- data.frame(url=paste(referenceSpectraFolder,list.files(referenceSpectraFolder, pattern="\\.csv$"),sep="/"),stringsAsFactors=F)
references$name <- unlist(lapply(references$url, function(x) {
    m <- regmatches(x,regexec("/([A-Za-z0-9_ ]+)\\.csv",x))
    if (length(m)) return(m[[1]][2])
    else return(x)
}))
# if necessary, here the references could be sorted, eg. 
# references <- references[references$name %in% c("HbOxy","HbMet","Hpx")]

# Read the reference spectra files
for (r in seq(dim(references)[1])) {
    if (r==1) {
        ref <- read.single.spectrum(references$url[r])
        colnames(ref)[2] <- references$name[r]
    } else {
        a <- read.single.spectrum(references$url[r])
        ref[[references$name[r]]] <- NA
        for (i in seq(dim(a)[1])) {
            w <- which(ref$wl==a$wl[i])
            if (length(w)) ref[[references$name[r]]][w] <- a$abs[i]
        }
    }
}

#include a precipitation spectrum
ref$precip <- 0.1
ref_name <- c(references$name, 'Precipitate')
ref <- na.omit(ref)

# Read the kinetics file
f <- readLines(kineticsFile)
cols <- as.vector(strsplit(f[[1]],",")[[1]])
measurement_id <- sapply(cols,function(x) {
	x <- regmatches(x,regexec('^(.+)_(\\d+)$',x))
	if (length(x)) return(as.numeric(x[[1]][3]))
	else return('')
})
col_id <- sapply(cols,function(x) {
	x <- regmatches(x,regexec('^(.+)_(\\d+)$',x))
	if (length(x)) return(x[[1]][2])
	else return('')
})
#Calculate number of samples, number of cycles from first line.
nSamples <- min(table(measurement_id))
nCycles <- floor(length(cols)/2/nSamples)
#skip the second line, then generate a empty list of matrices.
line <- 2

d <- list()
for (i in 1:nSamples) {
	d[[i]] <- matrix(ncol=nCycles)
}
names(d) <- col_id[2*(1:nSamples)-1]
wl <- c()
#iterate through measurements and load file
while(line<length(f)) {
	line <- line+1
	t <- f[[line]]
	if (nchar(t)<2) break
	t <- as.numeric(strsplit(t,",")[[1]])
	wl <- c(wl,t[1])
	for (i in 1:nSamples)
		d[[i]] <- rbind(d[[i]],t[0:(nCycles-1)*nSamples*2+(i)*2])
}
timepoints <- matrix(nrow=nSamples,ncol=nCycles)
#load timepoint from .csv appendix
count <- 0
while(line<length(f)) {
	line <- line+1
	t <- f[[line]]
	t <- regmatches(t,regexec('^\\[Time\\]\\s,\\s(\\S+)$',t))
	if (length(t[[1]])) {
		timepoints[1+count %% nSamples,1+ count %/% nSamples] <- as.numeric(t[[1]][2])
		count <- count+1
	}
}
rm(f)

# match recorded wl to reference wl
wl <- round(wl)
wlMatch <- match(wl,ref$wl)
#fall back for rounding errors
wlMatch[is.na(wlMatch)] <- match(wl[is.na(wlMatch)],ref$wl-1)
wlMatch[is.na(wlMatch)] <- match(wl[is.na(wlMatch)],ref$wl+1)
if (sum(is.na(wlMatch))>0) stop("Could not match all measured wavelengths to the reference spectra.")
# truncate references to selected references
ref <- as.matrix(ref[wlMatch,2:(dim(ref)[2])])
#perform nnls
#initialize result list
result <- list()
for (s in 1:nSamples) {
	m <- matrix(c(timepoints[,s],rep(NA,dim(ref)[2]*nCycles+nCycles)),byrow=F,ncol=dim(ref)[2]+2,nrow=nCycles)
	colnames(m) <- c("time",colnames(ref),"residue")
	result[[s]] <- m
}
names(result) <- names(d)
	
for (s in 1:nSamples)
	for (i in 1:nCycles) {
		r <- nnls::nnls(ref,as.vector(na.omit(d[[s]][,i])))
		result[[s]][i,2:(length(r$x)+2)] <- c(r$x,r$deviance)
	}
for (s in 1:nSamples) 
	write.csv(result[[s]],paste0(outputFileStem,s,".",names(result)[s],".csv"))