################################################## Part1 GSDataQC ###################################################
#' @export GSDataQC
GSDataQC <- function(markers, phenotype, impute = F, filter = F, NABound = 0.8, MAFBound = 0.05, imputeMethod = "median",hete = 1,silent =F){
  if (!is.matrix(markers) & !is.data.frame(markers)){
    markerFormat <- warning(paste0("The format of markers is wrong! Please check!","\n"))
  }else{
    if (is.data.frame(markers)) {
      markers <- as.matrix(markers)
    }
    markerFormat <- message("The format of markers is right!")
  }
  if (! is.numeric(markers)){
    markersClass <- message("The markers is not numeric! Please transform first with function 'transLetter2number' before doing genomic prediction.")
    format <- "letter"
  }else{
    markersClass <- message("The data in matrix is numeric!")
    format <- "number"
  }
  if (is.matrix(markers) & is.data.frame(phenotype)){
    if (nrow(markers) == nrow(phenotype)) {
      equal <- all(phenotype[,1] == rownames(markers))
      if (!equal) {
        warning(paste0("markers order does not correspond to phenotype order, please check!")) 
      }
    }else{
      warning(paste0("markers count does not correspond to phenotype count, please check!"))
    }
    
    MAF <- MAFSummary(markers,format = format,hete = hete,silent = silent)
    if (length(which(is.na(markers)) == TRUE) == 0){
      NACount <- 0 
      NAPercentTotal <- 0
      NACountRow <- 0
      NACountCol <- 0
      NACountRowPercent <- 0
      NACountColPercent <- 0
      resImpute <- markers
      NASummary <- message("The marker matrix has no missing value!")
      NAIdx <- NA
      
      if(filter){
        # NAFliIdx <- which(NACountColPercent > NABound)
        if (format == "number") {
          MAFFliIdx <- which(MAF < MAFBound)
        }else if(format == "letter"){
          MAFFliIdx <- which(MAF$MAF < MAFBound)
        }
        # message(paste0(length(NAFliIdx)," marker(s) was(were) filtered by missing rate > ",NABound,"."))
        message(paste0(length(MAFFliIdx)," marker(s) was(were) filtered by maf < ",MAFBound,"."))
        fliterIdx <- MAFFliIdx
        markers <- markers[,-fliterIdx]
      }
      
      if(impute){
        resImpute <- gsImputation(G = markers, imputeMethod = imputeMethod,silent = silent)
      }else{
        resImpute <- markers
      }
      # message("Statistical matrix composition ...")
      # tableRes <- table(markers)
    }else{
      NACount <- round(length(which(is.na(markers)) == TRUE))
      NAPercentTotal <- NACount/(ncol(markers)*nrow(markers))
      
      message("Calculating missing rate (MR) ...")
      if(silent){
        NACountRow <- apply(markers,1,function(x){length(which(is.na(x)) == TRUE)})
        NACountCol <- apply(markers,2,function(x){length(which(is.na(x)) == TRUE)})
      }else{
        NACountRow <- apply_pb(markers,1,function(x){length(which(is.na(x)) == TRUE)})
        NACountCol <- apply_pb(markers,2,function(x){length(which(is.na(x)) == TRUE)})
      }
      NACountRowPercent <- round(NACountRow/ncol(markers),2)
      NACountColPercent <- round(NACountCol/nrow(markers),2)
      
      NAIdx <- which(is.na(markers) == TRUE)
      
      # message("Statistical matrix composition ...")
      # tableRes <- table(markers)
      
      if(filter){
        NAFliIdx <- which(NACountColPercent > NABound)
        if (format == "number") {
          MAFFliIdx <- which(MAF < MAFBound)
        }else if(format == "letter"){
          MAFFliIdx <- which(MAF$MAF < MAFBound)
        }
        message(paste0(length(NAFliIdx)," marker(s) was(were) filtered by missing rate > ",NABound,"."))
        message(paste0(length(MAFFliIdx)," marker(s) was(were) filtered by maf < ",MAFBound,"."))
        fliterIdx <- unique(c(NAFliIdx,MAFFliIdx))
        markers <- markers[,-fliterIdx]
      }
      
      if(impute){
        resImpute <- gsImputation(G = markers, imputeMethod = imputeMethod,silent = silent)
      }else{
        resImpute <- markers
      }
    }
    
    phenotypeNACount <- length(which(is.na(phenotype) == TRUE))
    phenotypeNAIdx <- which(is.na(phenotype) == TRUE)
    
    NASum <- c(NACount,NAPercentTotal)
    names(NASum) <- c("missValueCount","missValuePercent") 
    
    # markersTable <- c(tableRes,NASum)
    markersTable <- NASum
    ## summarize
    #     summarize <- paste0(markerFormat,markersClass,"The data has ",NACount," missing value!","\n","Missing value percent:",NAPercentTotal,"\n",
    #                             "The phenotype data has ",phenotypeNACount," missing value!","\n",
    #                             "The data has ",length(tableRes)," element.","\n","The details in markerReport!",
    #                             "\n")
    markerReport <- list(Global = markersTable, resMarker = resImpute,MAF = MAF, NACountRow = NACountRow,NACountRowPercent = NACountRowPercent,NACountCol = NACountCol,
                         NACountColPercent = NACountColPercent, NAIdx = NAIdx,penotypeNAIdx = phenotypeNAIdx)
    cat(paste0(markerFormat,markersClass,"The data has ",NACount," missing value!","\n","Missing value percent:",round(NAPercentTotal,digits = 2),"\n",
               "The phenotype data has ",phenotypeNACount," missing value!","\n",
               "The details are in markerReport!",
               "\n"))
    markerReport
  }else{
    stop(paste0("The format of data is wrong! There may be the following cases:","\n","The format of marker is not matrix or dataframe;","\n",
                "The format of phenotype is not dataframe."))
  }
}
gsImputation <- function(G,imputeMethod = "median",silent = FALSE){
  if(is.numeric(G)){
    if(imputeMethod == "mean"){
      if(silent){
        res <- apply(G,2,function(x){
          x[which(is.na(x))] <- mean(x,na.rm=TRUE)
          x
        })
      }else{
        res <- apply_pb(G,2,function(x){
          x[which(is.na(x))] <- mean(x,na.rm=TRUE)
          x
        })
      }
    }else if(imputeMethod=="median"){
      if(silent){
        res <- apply(G,2,function(x){
          tt <- table(x)
          x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
          x
        })
      }else{
        res <- apply_pb(G,2,function(x){
          tt <- table(x)
          x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
          x
        })
      }
    }
    class(res) <- "numeric"
  }else{
    if(imputeMethod =="mean"){
      stop("Method 'mean' is not available for non-numeric vectors.",call. = FALSE)
    }else {
      if(silent){
        res <- apply(G,2,function(x){
          tt <- table(x)
          x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
          x
        })
      }else{
        res <- apply_pb(G,2,function(x){
          tt <- table(x)
          x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
          x
        })
      }
    }
  }
  return(res)
}
MAFSummary <- function(G,format = "letter",hete = 1,silent = FALSE){
  siSNP2biSNP <- function(x) {
    y <- rep(NA, length(x))
    y[which(x == "A")] <- "AA"
    y[which(x == "T")] <- "TT"
    y[which(x == "C")] <- "CC"
    y[which(x == "G")] <- "GG"
    y[which(x == "R")] <- "AG"
    y[which(x == "Y")] <- "CT"
    y[which(x == "S")] <- "CG"
    y[which(x == "W")] <- "AT"
    y[which(x == "K")] <- "GT"
    y[which(x == "M")] <- "AC"
    y[which(x == "+")] <- NA
    y[which(x == "0")] <- NA
    y[which(x == "-")] <- NA
    y[which(x == "N")] <- NA
    return(y)
  }
  if (format == "letter") {
    judge <- apply(G[,sample(1:ncol(G),round(ncol(G)*0.02))],2,function(x){
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        y <- 0.5
      }else{
        allele <- strsplit(sample(x,1),split = "")[[1]]
        allele.len <- length(allele)
        if (allele.len == 1) {
          y <- 1
        }else{
          y <- 0
        }
      }
    })
    if (sum(judge)/(round(ncol(G)*0.02)) > 0.9) {
      message("Input genotype is coded by singe-letter, if not, please check your data!")
      G.rownames <- rownames(G)
      G.colnames <- colnames(G)
      message("Input genotype is coded by singe-letter, converting single-letter to bi-letter ...")
      if (silent) {
        G <- apply(G,2,siSNP2biSNP)
      } else {
        G <- apply_pb(G,2,siSNP2biSNP)
      }
      rownames(G) <- G.rownames
      colnames(G) <- G.colnames
    }
    
    ## maf calculation
    message("Calculating minor allele frequency (MAF) ...")
    if (silent) {
      maf.info <- apply(G,2,function(x){
        x <- na.omit(x)
        fasta <- paste0(x,collapse = "")
        allele.array <- strsplit(fasta,split = "")
        allele.info <- sort(table(na.omit(allele.array)),decreasing = T)
        if(length(allele.info) < 2){
          af <- allele.info[1]/sum(allele.info)
          maf <- 1-af
          af.allele <- names(allele.info)[1]
          if (is.null(af.allele)) {
            af.allele <- NA
          }
          maf.allele <- NA
        }else if(length(allele.info) > 2){
          warning("Data have multi-allele loci, these loci will be removed!")
          af <- maf <- maf.allele <- af.allele <- "remove"
        }else{
          af <- allele.info[1]/sum(allele.info)
          maf <- allele.info[2]/sum(allele.info)
          maf.allele <- names(allele.info)[2]
          af.allele <- names(allele.info)[1]
        }
        res <- c(af.allele,maf.allele,af,maf)
        res
      })
    }else{
      maf.info <- apply_pb(G,2,function(x){
        x <- na.omit(x)
        fasta <- paste0(x,collapse = "")
        allele.array <- strsplit(fasta,split = "")
        allele.info <- sort(table(na.omit(allele.array)),decreasing = T)
        if(length(allele.info) < 2){
          af <- allele.info[1]/sum(allele.info)
          maf <- 1-af
          af.allele <- names(allele.info)[1]
          if (is.null(af.allele)) {
            af.allele <- NA
          }
          maf.allele <- NA
        }else if(length(allele.info) > 2){
          warning("Data have multi-allele loci, these loci will be removed!")
          af <- maf <- maf.allele <- af.allele <- "remove"
        }else{
          af <- allele.info[1]/sum(allele.info)
          maf <- allele.info[2]/sum(allele.info)
          maf.allele <- names(allele.info)[2]
          af.allele <- names(allele.info)[1]
        }
        res <- c(af.allele,maf.allele,af,maf)
        res
      })
    }
    ## remove multi-allele loci 
    if(length(which(maf.info == "remove",arr.ind = T)) > 0){
      message("Checking for markers with more than 2 alleles. If found will be removed.")
      remove.idx <- unique(which(maf.info == "remove",arr.ind = T)[,1])
      maf.info <- maf.info[-remove.idx,]
      G <- G[,-remove.idx]
    }
    maf.info <- t(maf.info)
    maf.info <- as.data.frame(maf.info)
    
    class(maf.info$V3) <- class(maf.info$V4) <- "numeric"
    names(maf.info) <- c("Ref","Alt","AF","MAF")
    
  } else if (format == "number"){
    if(!is.numeric(G)){
      stop("The input format is not equal with defined! please check your data!")
    }else{
      message("Calculating minor allele frequency (MAF) ...")
      if (silent) {
        maf.info <- apply(G,2,function(x){
          x <- na.omit(x)
          if (hete %in% x) {
            if (length(table(x)) == 2) {
              maf <- length(which(x == hete))/(2*length(x))
            }else{
              a <- length(which(x == hete))
              b <- sort(table(x[-which(x == hete)]))
              maf <- (b[[1]]*2+a)/(length(x)*2)
              maf
            }
          }else{
            a <- table(x)
            b <- sort(a)
            if (length(a) == 0) {
              maf <- NA
            }else{
              maf <- b[1]/length(x)
              if (maf > 0.5) {
                maf <- 1- maf
              }
            }
          }
          maf.info <- maf
        })
      }else{
        maf.info <- apply_pb(G,2,function(x){
          x <- na.omit(x)
          if (hete %in% x) {
            if (length(table(x)) == 2) {
              maf <- length(which(x == hete))/(2*length(x))
            }else{
              a <- length(which(x == hete))
              b <- sort(table(x[-which(x == hete)]))
              maf <- (b[[1]]*2+a)/(length(x)*2)
              maf
            }
          }else{
            a <- table(x)
            b <- sort(a)
            if (length(b) == 0) {
              maf <- NA
            }else{
              maf <- b[1]/length(x)
              if (maf > 0.5) {
                maf <- 1- maf
              }
            }
          }
          maf.info <- maf
        })
      }
    }
  }
  maf.info
}

apply_pb <- function(X, MARGIN, FUN, ...) {
  env <- environment()
  pb_Total <- sum(dim(X)[MARGIN])
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  wrapper <- function(...) {
    curVal <- get("counter", envir = env)
    assign("counter", curVal + 1, envir = env)
    setTxtProgressBar(get("pb", envir = env), curVal + 1)
    FUN(...)
  }
  res <- apply(X, MARGIN, wrapper, ...)
  close(pb)
  res
}