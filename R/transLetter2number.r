####################################################### Trans letter to number #################################################
#'@export transLetter2number 
transLetter2number <- function(G,maf=0,mr=0,silent = FALSE,ref.alt = NULL,impute = TRUE,imputeMethod = "median"){
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
  
  G <- as.matrix(G)
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
  
  ## missing rate
  if (mr > 0) {
    message("Calculating missing rate (MR) ...")
    if (silent) {
      missing.rate <- apply(G,2,function(x){
        NA.idx <- which(is.na(x))
        NA.rate <- length(NA.idx)/length(x)
        NA.rate
      })
      if (sum(missing.rate > mr) > 0) {
        G <- G[,-which(missing.rate > mr)]
        message(paste0(sum(missing.rate > mr)," marker(s) was(were) filtered by missing rate > ",mr,"."))
      }else{
        message("No marker was removed by missing rate with ",mr,".")
      }
    }else{
      missing.rate <- apply_pb(G,2,function(x){
        NA.idx <- which(is.na(x))
        NA.rate <- length(NA.idx)/length(x)
        NA.rate
      })
      if (sum(missing.rate > mr) > 0) {
        G <- G[,-which(missing.rate > mr)]
        message(paste0(sum(missing.rate > mr)," marker(s) was(were) filtered by missing rate > ",mr,"."))
      }else{
        message("None marker was removed by missing rate with ",mr,".")
      }
    }
  }
  
  ## MAF
  message("Calculating minor allele frequency (MAF) ...")
  if (silent) {
    maf.info <- apply(G,2,function(x){
      x <- na.omit(x)
      if (length(x) == 0) {
        warning("Data have column(s) with NA completely, these loci will be removed!")
        af <- maf <- maf.allele <- af.allele <- "remove"
      }else{
        fasta <- paste0(x,collapse = "")
        allele.array <- strsplit(fasta,split = "")
        allele.info <- sort(table(allele.array),decreasing = T)
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
      }
      res <- c(af.allele,maf.allele,af,maf)
      res
    })
  }else{
    maf.info <- apply_pb(G,2,function(x){
      x <- na.omit(x)
      if (length(x) == 0) {
        warning("Data have column(s) with NA completely, these loci will be removed!")
        af <- maf <- maf.allele <- af.allele <- "remove"
      }else{
        fasta <- paste0(x,collapse = "")
        allele.array <- strsplit(fasta,split = "")
        allele.info <- sort(table(allele.array),decreasing = T)
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
      }
      res <- c(af.allele,maf.allele,af,maf)
      res
    })
  }
  
  if(!is.null(ref.alt)){
    if (all(colnames(G) %in% rownames(ref.alt))) {
      ref.alt <- ref.alt[match(colnames(G),rownames(ref.alt)),]
      maf.info[1:2,] <- t(ref.alt)
    } else {
      stop("Please check the ref.alt input, some makers in G are not included in ref.alt!")
    }
    message("The ref and alt infomation which user defined is used in transformation ...")
  }
  
  maf.info <- t(maf.info) 
  
  ## remove multi-allele loci 
  if(length(which(maf.info == "remove",arr.ind = T)) > 0){
    message("Checking for markers with more than 2 alleles. If found will be removed.")
    remove.idx <- unique(which(maf.info == "remove",arr.ind = T)[,1])
    maf.info <- maf.info[-remove.idx,]
    G <- G[,-remove.idx]
  }
  
  maf.info <- as.data.frame(maf.info)
  class(maf.info$V3) <- class(maf.info$V4) <- "numeric"
  names(maf.info) <- c("Ref","Alt","AF","MAF")
  
  ## filter G by maf
  if(maf > 0){
    maf.filter.idx <- which(maf.info$MAF < maf)
    if (length(maf.filter.idx) > 0) {
      G <- G[,-c(which(maf.info$MAF < maf))]
      maf.info <- maf.info[-c(which(maf.info$MAF < maf)),]
      message(paste0(length(maf.filter.idx)," marker(s) was(were) filtered by maf < ",maf,"."))
    }else{
      message("None marker was removed by maf with ",maf,".")
    }
  }
  
  ## Transformation of letter to number
  message("Transformation is begainning ...")
  Ref <- maf.info$Ref
  lines.ID <- rownames(G)
  markers.ID <- colnames(G)
  
  G <- t(G)
  if (silent) {
    G.number <- apply(cbind(Ref, G), 1, function(x) {
      tmp <- gregexpr(pattern = x[1], text = x[-1], 
                      fixed = T)
      res <- as.integer(lapply(tmp, function(z) {
        ifelse(z[1] < 0,2,2 - length(z))
      }))
      return(res)
    })
  }else{
    G.number <- apply_pb(cbind(Ref, G), 1, function(x) {
      tmp <- gregexpr(pattern = x[1], text = x[-1], 
                      fixed = T)
      res <- as.integer(lapply(tmp, function(z) {
        ifelse(z[1] < 0,2,2 - length(z))
      }))
      return(res)
    })
  }
  
  ## impute 
  if (impute) {
    missing <- which(is.na(G.number))
    if (length(missing) > 0) {
      message(paste0("Imputing missing data with mode ",imputeMethod,"."))
      G.number <- gsImputation(G.number,imputeMethod = imputeMethod,silent = silent)
    }
  }else {
    message("Imputation not required. Be careful using non-imputed matrices in mixed model solvers.")
  }
  
  rownames(G.number) <- lines.ID
  colnames(G.number) <- markers.ID
  final.res <- list(G.n = G.number,maf.info = maf.info)
  final.res
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