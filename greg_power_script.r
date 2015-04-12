#!/usr/bin/env Rscript

#grab the bottomly dataset
#http://bowtie-bio.sourceforge.net/recount/
#do this on a 5x5 dataset for good power


#this script was originally written by Dr. Greg Gloor, but has been completely munged by Ruth.

options(error=recover)


getHeatmapAndFDR <- function (nSamplesPerCondition) {




    e <- rea
d.table("data/bottomly_count_table.txt", sep="\t", row.names=1, header=T)

    #remove rows that contain less than 1 counts total
    e.rsum <- apply(e, 1, sum)
    d <- e[which(e.rsum >= 10),] #reduces to 11870 features from 36536

    library(ALDEx2)

    #library(DESeq)

    #minimum number
    #f <- c(1,2,4,8,16,32,64,128,256,512,1024)
    #base abundance
    f <- c(1,4,16,64,256,1024)
    foldf <- c(16,8,4,2,1)

    logBaseAbundance <- as.character(c(0:10)[c(TRUE,FALSE)])
    powerTestConditions <- c(rep("c",length(logBaseAbundance)),rep("e",length(logBaseAbundance)))
    foldChange <- rep(4:0,each=length(powerTestConditions))
    testNames <- paste(powerTestConditions,logBaseAbundance,foldChange,sep="|")



    #foldf <- c(1, 2, 4, 8, 16)
    "rdirichlet" <-
      function(n, alpha) {
        l <- length(alpha)
        x <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
        sm <- x%*%rep(1,l)
        return(x/as.vector(sm))
      }

    cnames.e <- colnames(e[,11:21])

    #for each fold difference
    # for( j in 1:length(foldf) ){
    #         f1 <- f * foldf[j] #make the fold difference table

    ##use this when you want the two test samples to be different conditions

    #       h1 <- c(d$B6033480, f, f1)
    #       h2 <- c(d$B6033488, f, f1)
    #       h3 <- c(d$B6033481, f, f1)
    #       h4 <- c(d$B6033489, f, f1)
    #       h5 <- c(d$B6033482, f, f1)
    #       b1 <- c(d$D2033485, f1, f)
    #       b2 <- c(d$D2033493, f1, f)
    #       b3 <- c(d$D2033486, f1, f)
    #       b4 <- c(d$D2033484, f1, f)
    #       b5 <- c(d$D2033492, f1, f)
    #
    #generate 10 sets for testing, as this is big


    #consider doing replicates later
    #        for(i in 1:2){
                    ##use this when you want the two test samples to be randomly chosen from the same condition
                    a <- sample(cnames.e, 2*nSamplesPerCondition)

                    #each column is a sample, each row is a feature
                    h <- matrix(ncol=nSamplesPerCondition,nrow=(length(d[,a[1]])+length(testNames)))
                    b <- matrix(ncol=nSamplesPerCondition,nrow=(length(d[,a[1]])+length(testNames)))
                    colnames(h) <- paste("h",c(1:nSamplesPerCondition),sep="")
                    colnames(b) <- paste("b",c(1:nSamplesPerCondition),sep="")
                    rownames(h) <- rownames(b) <- c(rownames(d),testNames)
                    for (j in 1:nSamplesPerCondition) {
                        h[,j] <- c(d[,a[j]], f, f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5])
                        b[,j] <- c(d[,a[(nSamplesPerCondition+j)]], f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5], f)
                    }


                    # h1 <- c(d[,a[1]], f, f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5])
                    # h2 <- c(d[,a[2]], f, f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5])
                    # h3 <- c(d[,a[3]], f, f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5])
                    # h4 <- c(d[,a[4]], f, f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5])
                    # h5 <- c(d[,a[5]], f, f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5])

                    # b1 <- c(d[,a[6]], f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5], f)
                    # b2 <- c(d[,a[7]], f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5], f)
                    # b3 <- c(d[,a[8]], f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5], f)
                    # b4 <- c(d[,a[9]], f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5], f)
                    # b5 <- c(d[,a[10]], f*foldf[1], f, f*foldf[2], f, f*foldf[3], f, f*foldf[4], f, f*foldf[5], f)

                    #generate the sampled dataset
                    #the 11870 genes are constant across all samples
                    #so any found are by definition false positives if they occur in the first
                    #the last 22 genes differ by x fold, with expressions that differ as in f1,f2
                    #so any in the last 22 that are found are true positives by definition

                    c <- matrix(ncol=nSamplesPerCondition,nrow=(length(d[,a[1]])+length(testNames)))
                    e <- matrix(ncol=nSamplesPerCondition,nrow=(length(d[,a[1]])+length(testNames)))
                    colnames(c) <- paste("c",c(1:nSamplesPerCondition),sep="")
                    colnames(e) <- paste("e",c(1:nSamplesPerCondition),sep="")
                    rownames(c) <- rownames(e) <- c(rownames(d),testNames)

                    for (j in 1:nSamplesPerCondition) {
                        c[,j] <- as.integer(sum(h[,j]) * (rdirichlet(1, (h[,j] + 0.5))))
                        e[,j] <- as.integer(sum(b[,j]) * (rdirichlet(1, (b[,j] + 0.5))))
                    }


                    # c1 <- as.integer(sum(h1) * (rdirichlet(1, (h1 + 0.5))))
                    # c2 <- as.integer(sum(h2) * (rdirichlet(1, (h2 + 0.5))))
                    # c3 <- as.integer(sum(h3) * (rdirichlet(1, (h3 + 0.5))))
                    # c4 <- as.integer(sum(h4) * (rdirichlet(1, (h4 + 0.5))))
                    # c5 <- as.integer(sum(h5) * (rdirichlet(1, (h5 + 0.5))))

                    # e1 <- as.integer(sum(b1) * (rdirichlet(1, (b1 + 0.5))))
                    # e2 <- as.integer(sum(b2) * (rdirichlet(1, (b2 + 0.5))))
                    # e3 <- as.integer(sum(b3) * (rdirichlet(1, (b3 + 0.5))))
                    # e4 <- as.integer(sum(b4) * (rdirichlet(1, (b4 + 0.5))))
                    # e5 <- as.integer(sum(b5) * (rdirichlet(1, (b5 + 0.5))))

            #ALDEx

                    aldex.in <- data.frame(c,e)
                    names(aldex.in) <- c(colnames(c),colnames(e))
                    rownames(aldex.in) <- c(rownames(d),testNames)
                    conds <- c(rep("C",nSamplesPerCondition),rep("E",nSamplesPerCondition))


                    # aldex.in <- data.frame(c1,c2,c3,c4,c5,e1,e2,e3,e4,e5)
                    # names(aldex.in) <- c("c1","c2","c3","c4","c5","e1","e2","e3","e4","e5")
                    # rownames(aldex.in) <- c(rownames(d), testNames)
                    # conds <- c("C","C","C","C","C","E","E","E","E","E")

                    x <- aldex(aldex.in, conds)
                    # y <- summary.aldex(x, digits=6, median.only=TRUE)

                    # #sig <- which(y$criteria.significant == TRUE)
                    # mean2 <-  c(foldf[j], rownames(y)[which(y$criteria.we.lfdr <= 0.1)] )
                    # mean3 <-  c(foldf[j], rownames(y)[which(y$criteria.we.BH <= 0.1)] )
                    # mean4 <-  c(foldf[j], y[nf,"criteria.we.pval"] )

                    #write(sig, file = "ALDEX_sig.txt", ncolumns = length(sig), append = TRUE, sep = " ")
                    #write(mean2, file = "../analysis/ALDEx_lfdr_bottomly_SAME_COND_0.1_100.txt", ncolumns = length(mean2), append = TRUE, sep = " ")

                    testGenes <- x[testNames,]

                    #write.table(testGenes, file = "analysis/ALDEx_wepval_bottomly_SAME_COND_0.1_100.txt", sep="\t")
                    write.table(testGenes, file = "analysis/ALDEx_wepval_bottomly_SAME_COND_test_genes.txt", sep="\t")


                    #get false positives
                    realGenes <- x[c(1:(nrow(x)-length(testNames))),]

                    falsePositives005 <- realGenes[which(realGenes$we.eBH < 0.05),]
                    falsePositives01 <- realGenes[which(realGenes$we.eBH < 0.1),]

                    write.table(falsePositives005,file="analysis/ALDEx_wepval_bottomly_SAME_COND_false_positives_005.txt",sep="\t")
                    write.table(falsePositives01,file="analysis/ALDEx_wepval_bottomly_SAME_COND_false_positives_01.txt",sep="\t")

                    # #examine false positives
                    # print(falsePositives)

                    #output p-values in table
                    testPValues <- matrix(nrow=length(f),ncol=length(foldf),testGenes$we.eBH)
                    rownames(testPValues) <- f
                    colnames(testPValues) <- foldf

                    write.table(testPValues,file="analysis/pvalues_for_base_abundance_rows_vs_fold_change_cols.txt",sep="\t")

                    #output pvalues as heatmap per fold change per abundance
                    library(gplots)
                    library(RColorBrewer)
                    heatmapPalette <- colorRampPalette(c("deeppink3", "darkorchid1", "gray66"))(n = 299)
                    col_breaks = c(seq(0,0.05,length=100),  # blue is significant at 0.05
                      seq(0.05,0.1,length=100), #cyan is significant at 0.1
                      seq(0.1,1,length=100)) #the rest is yellow

                    # creates a 5 x 5 inch image
                    png("analysis/pvalue_heatmap.png",    # create PNG for the heat map        
                      width = 5*300,        # 5 x 300 pixels
                      height = 5*300,
                      res = 300,            # 300 pixels per inch
                      pointsize = 8)        # smaller font size




                    heatmap.2(testPValues,
                      cellnote = round(testPValues,digits=4),  # same data set for cell labels
                      main = "P Values", # heat map title
                      notecol="black",      # change font color of cell labels to black
                      density.info="none",  # turns off density plot inside color legend
                      trace="none",         # turns off trace lines inside the heat map
                     # margins =c(12,9),     # widens margins around plot
                      col=heatmapPalette,       # use on color palette defined earlier 
                      dendrogram='none',
                      breaks=col_breaks, 
                      Rowv=FALSE,
                      Colv=FALSE,
                      xlab="FOLD CHANGE",
                      ylab="BASE ABUNDANCE")            # turn off column clustering

                    dev.off()               # close the PNG device





            ##
    #       #DESeq
    #               countsTable <- aldex.in
    #               conds <- factor(conds)
    #       #
    #               cds <- newCountDataSet( countsTable, conds )
    #               cds <- estimateSizeFactors( cds )
    #
    #               #corresponds to DESeq when we started
    #               cds <- estimateDispersions( cds, sharingMode="fit-only" , fitType="local")
    #               res <- nbinomTest( cds, "C", "E")
    #
    #               d.sig <-  c(foldf[j], res$id[which(res$padj < 0.1)])
    #               write (d.sig, file="../analysis/DESeq_bottomly_SAME_COND_fdr10_100.txt", ncolumns = length(d.sig), append = TRUE, sep = " ")
    #       #       d.sig <- (res)[which(res$padj < 0.05)]
            #       write (d.sig, file="DESeq_fdr5.txt", ncolumns = length(d.sig), append = TRUE, sep = " ")
            #
    #}


    #}



    returnList <- list()
    returnList$testPValues <- testPValues
    returnList$falsePositives005 <- falsePositives005
    returnList$falsePositives01 <- falsePositives01

    return(returnList)


}