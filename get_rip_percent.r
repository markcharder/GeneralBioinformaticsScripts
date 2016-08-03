require("seqinr")
require("Hmisc")
testtwospeed <- function(windowsize, sequence, rips){

    if (windowsize < length(sequence)) {

        starts <<- seq(1, length(sequence) - windowsize, by = 1000)

        if (exists("seqsize")){

            positions <<- starts + seqsize

        }

        else {

            positions <<- starts

        }

        if (exists("seqsize")){

            seqsize <<- seqsize + length(sequence)

        }

        else{

            seqsize <<- length(sequence)

        }

        n <- length(starts)

        chunkRIPs <<- numeric(n)

        if (!exists("scaffolds")){

            subRIPs <- rips
        }

        else {

            subRIPs <- grep(paste(currentscaffold, "\t", sep = "", collapse = ""), 
            rips, value = TRUE)
            print(length(subRIPs))

        }

        for (i in 1:n){

            chunk <- sequence[starts[i]:(starts[i] + windowsize - 1)]
            start <- starts[i]
            end <- starts[i] + windowsize - 1
            size <- 0
            totalRIPs <- 0

            for (j in 1:length(subRIPs)){

                fields <- unlist(strsplit(subRIPs[j], split = "\\t"))
                featurestart <- as.numeric(fields[4])
                featureend <- as.numeric(fields[5])


                if (length(fields) > 1){

                    if ((featurestart > start) && (featureend < end)){

                        size <- featureend - featurestart
                        totalRIPs <- totalRIPs + size

                    }

                    if ((featurestart > start) && (featurestart < end) && (featureend > end)){

                        size <- end - featurestart
                        totalRIPs <- totalRIPs + size

                    }

                    if ((featurestart < start) && (featureend < end) && (featureend > start)){

                        size <- featureend - start
                        totalRIPs <- totalRIPs + size

                    }

                }

            }


            chunkRIPs[i] <<- totalRIPs / windowsize
            print(paste("Adding rip :", chunkRIPs[i]))


        }

    }


    else{

        print("Warning: window size larger than or equal to at least one scaffold")

    }

}

testtwospeedmultiplescaffolds <- function(prefix, windowsize){

    fasta <- read.fasta(paste(prefix, ".fa", sep = "", collapse = ""))
    a <- file(paste(prefix, ".gff3", sep = "", collapse = ""), open = "r")

    rips <- readLines(a)

    multirips <<- c()
    multistarts <<- c()

    plainfasta <- read.table(paste(prefix,".fa", sep = "", collapse = ""))
    scaffolds <<- grep(">", plainfasta$V1,value=TRUE)
    scaffolds <<- gsub(">", "", scaffolds)

    for (i in 1:length(fasta)){

        currentscaffold <<- scaffolds[i]
        testtwospeed(windowsize, fasta[[i]], rips)
        print(paste(i, "done"))
        multirips <<- c(multirips, chunkRIPs)
        multistarts <<- c(multistarts, positions)

    }

    close(a)

}

testtwospeedmultiplegenomes <- function(genomes, windowsize){

        comparativegenomicsrips <<- list()
        comparativegenomicstarts <<- list()

    for (i in 1:length(genomes)){

        testtwospeedmultiplescaffolds(genomes[i], windowsize)
        comparativegenomicsrips[[i]] <<- multirips
        comparativegenomicstarts[[i]] <<- multistarts

    }

}

