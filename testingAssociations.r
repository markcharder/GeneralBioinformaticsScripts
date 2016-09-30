readFiles <- function(){
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
readFiles <- function(){
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

barplotData <- function(number){
  x<-barplot(
    c(mean(unlist(secretedClosest[[number]])),
    mean(unlist(noSecClosest[[number]])),
    mean(unlist(effectorClosest[[number]])),
    mean(unlist(noEffClosest[[number]]))),
    beside=TRUE,
    ylim=c(0,400))

  arrows(
    x[,number],
    c(mean(unlist(secretedClosest[[number]]))-((sd(unlist(secretedClosest[[number]])))/2),
    mean(unlist(noSecClosest[[number]]))-((sd(unlist(noSecClosest[[number]])))/2),
    mean(unlist(effectorClosest[[number]]))-((sd(unlist(effectorClosest[[number]])))/2),
    mean(unlist(noEffClosest[[number]]))-((sd(unlist(noEffClosest[[number]])))/2)),
    x[,number],c(mean(unlist(secretedClosest[[number]]))+((sd(unlist(secretedClosest[[number]])))/2),
    mean(unlist(noSecClosest[[number]]))+((sd(unlist(noSecClosest[[number]])))/2),
    mean(unlist(effectorClosest[[number]]))+((sd(unlist(effectorClosest[[number]])))/2),
    mean(unlist(noEffClosest[[number]]))+((sd(unlist(noEffClosest[[number]])))/2)),
    code=3)
  }


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
