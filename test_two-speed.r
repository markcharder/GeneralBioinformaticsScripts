require("seqinr")
require("Hmisc")
testtwospeed <- function(windowsize, sequence, genes, CDSs, repeats, effectors){

    if (windowsize < length(sequence)) {
    
        starts <- seq(1, length(sequence) - windowsize, by = windowsize)


        n <- length(starts)

        chunkGCs <<- numeric(n)
        chunkRepeats <<- numeric(n)
        chunkGenes <<- numeric(n)
        chunkEffectors <<- numeric(n)
        chunkCDSs <<- numeric(n)

        if (!exists("scaffolds")){

            subgenes <- genes
            subrepeats <- repeats
            subcds <- CDSs
            subeffectors <- effectors
        }

        else {

            subgenes <- grep(paste(currentscaffold, "\t", sep = "", collapse = ""), 
            genes, value = TRUE)
            print(length(subgenes))
            subcds <- grep(paste(currentscaffold, "\t", sep = "", collapse = ""), 
            CDSs, value = TRUE)
            subrepeats <- grep(paste(currentscaffold, "\t", sep = "", collapse = ""), 
            repeats, value = TRUE)
            subeffectors <- grep(paste(currentscaffold, "\t", sep = "", collapse = ""), 
            effectors, value = TRUE)

        }

        for (i in 1:n){

            chunk <- sequence[starts[i]:(starts[i] + windowsize - 1)]
            start <- starts[i]
            end <- starts[i] + windowsize - 1
            size <- 0
            genecount <- 0
            effectorcount <- 0
            totalCDSs <- 0
            totalrepeats <- 0

            for (j in 1:length(subgenes)){

                fields <- unlist(strsplit(subgenes[j], split = "\\t"))
                featurestart <- as.numeric(fields[4])
                featureend <- as.numeric(fields[5])

                if ((featurestart > start) && (featureend < end)){

                genecount <- genecount + 1

                }

            }

            for (j in 1:length(subcds)){

                fields <- unlist(strsplit(subcds[j], split = "\\t"))
                featurestart <- as.numeric(fields[4])
                featureend <- as.numeric(fields[5])

                if ((featurestart > start) && (featureend < end)){

                    size <- featureend - featurestart
                    totalCDSs <- totalCDSs + size

                }

                if ((featurestart > start) && (featurestart < end) && (featureend > end)){

                    size <- end - featurestart
                    totalCDSs <- totalCDSs + size

                }

                if ((featurestart < start) && (featureend < end) && (featureend > start)){

                    size <- featureend - start
                    totalCDSs <- totalCDSs + size

                }

            }


            for (j in 1:length(subrepeats)){

                fields <- unlist(strsplit(subrepeats[j], split = "\\t"))
                featurestart <- as.numeric(fields[4])
                featureend <- as.numeric(fields[5])

                if ((featurestart > start) && (featureend < end)){

                    size <- featureend - featurestart
                    totalrepeats <- totalrepeats + size

                }

                if ((featurestart > start) && (featurestart < end) && (featureend > end)){

                    size <- end - featurestart
                    totalrepeats <- totalrepeats + size

                }

                if ((featurestart < start) && (featureend < end) && (featureend > start)){

                    size <- featureend - start
                    totalrepeats <- totalrepeats + size

                }

            }


            for (j in 1:length(subeffectors)){

                fields <- unlist(strsplit(subeffectors[j], split = "\\t"))
                featurestart <- as.numeric(fields[4])
                featureend <- as.numeric(fields[5])

                if ((featurestart > start) && (featureend < end)){

                effectorcount <- effectorcount + 1

                }

            }

            chunkGCs[i] <<- GC(chunk)
            chunkGenes[i] <<- genecount
            chunkCDSs[i] <<- totalCDSs / windowsize
            chunkRepeats[i] <<- totalrepeats / windowsize
            chunkEffectors[i] <<- effectorcount

            print(currentscaffold)
            print(paste("Adding genes :", chunkGenes[i]))   
            print(paste("Adding gc :", chunkGCs[i]))   
            print(paste("Adding CDSs :", chunkCDSs[i]))   
            print(paste("Adding effectors :", chunkEffectors[i]))   
            print(paste("Adding repeats :", chunkRepeats[i]))   

        }

    }

    else{

        print("Warning: window size larger than or equal to at least one scaffold")

    }

}

testtwospeedmultiplescaffolds <- function(prefix, windowsize){

    fasta <- read.fasta(paste(prefix, ".fa", sep = "", collapse = ""))
    a <- file(paste(prefix, ".repeats", sep = "", collapse = ""), open = "r")
    b <- file(paste(prefix, ".CDSs", sep = "", collapse = ""), open = "r")
    c <- file(paste(prefix, ".genes", sep = "", collapse = ""), open = "r")
    d <- file(paste(prefix, ".effectors", sep = "", collapse = ""), open = "r")

    repeats <- readLines(a)
    cds <- readLines(b)
    genes <- readLines(c)
    effectors <- readLines(d)

    multirepeats <<- c()
    multiCDSs <<- c()
    multigenes <<- c()
    multieffectors <<- c()
    multiGC <<- c()

    plainfasta <- read.table(paste(prefix,".fa", sep = "", collapse = ""))
    scaffolds <<- grep(">", plainfasta$V1,value=TRUE)
    scaffolds <<- gsub(">", "", scaffolds)

    for (i in 1:length(fasta)){

        currentscaffold <<- scaffolds[i]
        testtwospeed(windowsize, fasta[[i]], genes, cds, repeats, effectors)
        print(paste(i, "done"))
        multirepeats <<- c(multirepeats, chunkRepeats)
        multiCDSs <<- c(multiCDSs, chunkCDSs)
        multigenes <<- c(multigenes, chunkGenes)
        multieffectors <<- c(multieffectors, chunkEffectors)
        multiGC <<- c(multiGC, chunkGCs)

    }

    close(a)
    close(b)
    close(c)
    close(d)

}

testtwospeedmultiplegenomes <- function(genomes, windowsize){

        comparativegenomicsrepeats <<- list()
        comparativegenomicsCDSs <<- list()
        comparativegenomicsgenes <<- list()
        comparativegenomicseffectors <<- list()
        comparativegenomicGC <<- list()

    for (i in 1:length(genomes)){

        testtwospeedmultiplescaffolds(genomes[i], windowsize)
        comparativegenomicsrepeats[[i]] <<- multirepeats
        comparativegenomicsCDSs[[i]] <<- multiCDSs
        comparativegenomicsgenes[[i]] <<- multigenes
        comparativegenomicseffectors[[i]] <<- multieffectors
        comparativegenomicGC[[i]] <<- multiGC

    }

}

spearmancorrelation <- function(number){

    for (i in 1:number){

        effectors <- comparativegenomicseffectors[[i]]
        genes <- comparativegenomicsgenes[[i]]
        repeats <- comparativegenomicsrepeats[[i]]
        cds <- comparativegenomicsCDSs[[i]]
        gc <- comparativegenomicGC[[i]]

        effectors[effectors == 0] <- NA
        genes[genes == 0] <- NA
        repeats[repeats == 0] <- NA
        cds[cds == 0] <- NA
        gc[gc == 0] <- NA
        ratio <- effectors / genes

        setEPS()
        postscript(paste(i, "_cds.eps",sep = "", collapse = ""))
        plot(ratio, cds, xlab = "Secreted / non-secreted", 
        ylab = "CDS content", bty = "l", pch = 20)
        abline(lm(cds ~ ratio), col = "blue")
        dev.off()
        setEPS()

        postscript(paste(i, "_repeats.eps",sep = "", collapse = ""))
        plot(ratio, repeats, xlab = "Secreted / non-secreted", 
        ylab = "Repeat content", bty = "l", pch = 20)
        abline(lm(repeats ~ ratio), col = "blue")
        dev.off()
        setEPS()

        postscript(paste(i, "_gc.eps",sep = "", collapse = ""))
        plot(ratio, gc, xlab = "Secreted / non-secreted", 
        ylab = "GC content", bty = "l", pch = 20)
        abline(lm(gc ~ ratio), col = "blue")
        dev.off()
        setEPS()

        sink(file = "stats.txt", append = TRUE)
        print(paste(i, "Results for CDS"))
        print(rcorr(ratio, cds, type = "spearman"))
        print(paste(i, "Results for repeats"))
        print(rcorr(ratio, repeats, type = "spearman"))
        print(paste(i, "Results for GC"))
        print(rcorr(ratio, gc, type = "spearman"))
        sink(file = NULL)

    }

}
