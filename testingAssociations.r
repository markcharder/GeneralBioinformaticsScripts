## First is for repeats.
readFilesRepeast <- function(){
  effectorClosest<<-list()
  secretedClosest<<-list()
  noEffClosest<<-list()
  noSecClosest<<-list()
  fileArray <<- c()
  noeffFiles <- list.files(pattern = "*noEffectors.closest")
  effectorFiles <- list.files(pattern = "*effectors.closest")
  nosecFiles <- list.files(pattern = "*noSecreted.closest")
  secretedFiles <- list.files(pattern = "*secreted.closest")
  for (i in 1:length(secretedFiles)){
    fileArray[i] <<- secretedFiles[i]
    print(paste("Processing file: ", secretedFiles[i]))
    secreted <- read.table(secretedFiles[i], header=F)
    secretedClosest[[i]] <<- secreted
    noSec <- read.table(nosecFiles[i], header=F)
    noSec <- sample(noSec$V1, length(secreted$V1))
    noSecClosest[[i]] <<- noSec
    }
  for (i in 1:length(effectorFiles)){
    print(paste("Processing file: ", effectorFiles[i]))
    effectors <- read.table(effectorFiles[i], header=F)
    effectorClosest[[i]] <<- effectors
    noEff <- read.table(noeffFiles[i], header=F)
    noEff <- sample(noEff$V1, length(effectors$V1))
    noEffClosest[[i]] <<- noEff
    }
  }

## Second is for RIPs.
readFilesRIPs <- function(){
  effectorClosest<<-list()
  secretedClosest<<-list()
  noEffClosest<<-list()
  noSecClosest<<-list()
  fileArray <<- c()
  noeffFiles <- list.files(pattern = "*noEffectors.rips.closest")
  effectorFiles <- list.files(pattern = "*effectors.rips.closest")
  nosecFiles <- list.files(pattern = "*noSecreted.rips.closest")
  secretedFiles <- list.files(pattern = "*secreted.rips.closest")
  for (i in 1:length(secretedFiles)){
    fileArray[i] <<- secretedFiles[i]
    print(paste("Processing file: ", secretedFiles[i]))
    secreted <- read.table(secretedFiles[i], header=F)
    secretedClosest[[i]] <<- secreted
    noSec <- read.table(nosecFiles[i], header=F)
    noSec <- sample(noSec$V1, length(secreted$V1))
    noSecClosest[[i]] <<- noSec
    }
  for (i in 1:length(effectorFiles)){
    print(paste("Processing file: ", effectorFiles[i]))
    effectors <- read.table(effectorFiles[i], header=F)
    effectorClosest[[i]] <<- effectors
    noEff <- read.table(noeffFiles[i], header=F)
    noEff <- sample(noEff$V1, length(effectors$V1))
    noEffClosest[[i]] <<- noEff
    }
  }
require("ggplot2")

##Plot data using ggplot2.
plotData<-function(number){
  nosec<-data.frame(unlist(noSecClosest[[number]]),
    rep("Non-secreted",length(unlist(noSecClosest[[number]]))))
  sec<-data.frame(unlist(secretedClosest[[number]]),
    rep("Secreted",length(unlist(secretedClosest[[number]]))))
  noef<-data.frame(unlist(noEffClosest[[number]]),
    rep("Non-effector",length(unlist(noEffClosest[[number]]))))
  ef<-data.frame(unlist(effectorClosest[[number]]),
    rep("Effector",length(unlist(effectorClosest[[number]]))))
  colnames(nosec)<-c("a","b")
  colnames(sec)<-c("a","b")
  colnames(noef)<-c("a","b")
  colnames(ef)<-c("a","b")
  datagg<<-rbind(nosec,sec,noef,ef)

  ggplot(datagg,aes(b,a)) + xlab("Gene type") + ylab("Distance from nearest repeat") + geom_violin()
  }
