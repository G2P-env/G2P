multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

.sim <- function(markers,mark6.centered,mark.svd,tp.size,vp.names,method) {
  mark6.centered  <- mark6.centered
  mark.svd <- mark.svd
  po <- mark.svd$d^2/sum(mark.svd$d^2)
  ix <- which(mark.svd$d < 1e-09)
  if (length(ix) > 0) {
    mark6.PCbasis <- mark6.centered %*% mark.svd$v[, -ix]
  }else {
    mark6.PCbasis <- mark6.centered %*% mark.svd$v
  }
  VP <- which(rownames(markers) %in% vp.names)
  names(VP) <- rownames(mark6.PCbasis)[VP]
  ALLDIST <- dist(mark6.PCbasis[, 1:2])
  ALLDIST2 <- as.matrix(ALLDIST)
  ALLDIST2[1:5, 1:5]
  remove.list <- list(NA)
  for (k in 1:length(tp.size)) {
    fff <- tp.size[k]
    remove <- 1
    starto <- 1
    is.even <- function(x) {
      x%%2 == 0
    }
    while (length(remove) < fff) {
      potential <- apply(as.data.frame(VP), 1, function(x,y, z) {
        if (method == "sim") {
          closer <- (sort(y[x, ], decreasing = FALSE))[-1]
          best <- closer[1:z]
        }
        else if (method == "sim-dissim") {
          closer1 <- (sort(y[x, ], decreasing = FALSE))[-1]
          closer2 <- (sort(y[x, ], decreasing = TRUE))[-1]
          best <- c(closer1[1:z], closer2[1:z])
        }
        return(names(best))
      }, y = ALLDIST2, z = starto)
      TP <- unique(as.vector(unlist(potential)))
      VPP <- names(VP)
      remove <- setdiff(TP, VPP)
      starto <- starto + 1
    }
    actual <- length(remove) - tp.size[k]
    if (actual > 0) {
      xxx <- tp.size[k]
      remove <- remove[1:xxx]
    }
    remove.list[[k]] <- remove
  }
  tp.list <- remove.list
}

########################################### SPTPSel function ##################################################
#' @export TSRefine
#' @import STPGA qtlDesign reshape pbapply
#' @importFrom ggplot2 margin

TSRefine <- function(markers,candidates,test = NULL,ntosel,method = "PEVmean",npop = 100, nelite =5 ,mutprob = .8,niterations = 500,
                    cores = 1,lambda = NULL,sequent = FALSE,visulization = F){
  suppressMessages(require("STPGA"))
  suppressMessages(require("qtlDesign"))
  suppressMessages(require("ggplot2"))
  suppressMessages(require("reshape"))
  suppressMessages(require("pbapply"))
  pboptions(type = "timer")
  
  scaleRes <- scale(markers, center = T, scale = F)
  svdRes <- svd(scaleRes, nu=50, nv=50)
  idx <- which(svdRes$d < 1e-09)
  PCbasis <- scaleRes %*% svdRes$v[,-idx]
  po <- svdRes$d^2/sum(svdRes$d^2)
  
  if(is.null(test)){
    res <- list()
    length(res) <- length(ntosel)
    cat("Phenotyping selection is doing ... \n")
    cat(ntosel, "samples will be selected from markers ... \n")
    res <- pblapply(1:length(ntosel),function(i){
      resIdx <- mma(markers[candidates,],p = ntosel[i],sequent =sequent)[[1]]
      candidates[resIdx]
    })
  }else{
    if(method %in% c("PEVmean","RD","CDmean","sim")){
      method <- method
      cat("Selection of trainning population is doing ... \n")
      cat(ntosel, " samples will be selected from candidates with method", method,"... \n")
    }else{
      stop(cat("The methods selection is wrong!"))
    }
    res <- list()
    length(res) <- length(ntosel)
    
    if(method == "RD"){
      res <- pblapply(1:length(ntosel),function(i){
        resRD <- sample(candidates,ntosel[i])
        resRD
      })
      #     }else if(method == "sim"){
      #       res <- .TPO(markers = markers, vp.names = test, tp.size = ntosel, method = "sim", 
      #            npop = npop, nelite = nelite, mutprob = mutprob, niterations = niterations, 
      #            lambda = lambda)
    }else{
      if(is.null(lambda)){
        LambdaTrait <- 1/ncol(markers)
      }else{
        LambdaTrait <- lambda
      }
      if(method == "sim"){
        res <- pblapply(1:length(ntosel),function(i){
          .sim(markers = markers,mark6.centered = scaleRes,mark.svd =svdRes,tp.size =ntosel[i],vp.names = test,method = method)[[1]]
        })
      }else{
        rm(markers)
        if(method == "PEVmean"){
          errorstat = "PEVMEAN"
        }else if(method == "CDmean"){
          errorstat = "CDMEAN"
        }
        
        res <- pblapply(1:length(ntosel),function(i){
          suppressWarnings(res1 <- GenAlgForSubsetSelection(P = PCbasis, Candidates = candidates, Test = test, ntoselect = ntosel[i], npop = npop, nelite =nelite, keepbest = TRUE, tabu = T, tabumemsize = 1, mutprob = mutprob,
                                           mutintensity = 1, niterations = niterations,
                                           minitbefstop = 200, niterreg = 5, lambda = LambdaTrait,
                                           plotiters = FALSE, plottype=1,errorstat = errorstat, C = NULL,
                                           mc.cores = cores, InitPop = NULL, tolconv = 1e-07, Vg = NULL, Ve = NULL))
          res1[[1]]
        })
      }
    }
    res
  }
  if(visulization){
    cat("Plotting ...")
    po_plot <- melt(po)
    g0 <- ggplot(po_plot,aes(x = 1:length(po),y = value)) + geom_point(size =1 ,color = "red") + 
      labs(x = "Principal components",y = "Percent variation explained",title =  paste("Var explained 3 PCs = ", round(sum(po[1:3]), 3))) +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",color = "black"),
            plot.title = element_text(hjust = 1,size = 14),axis.text = element_text(size = 14),
            axis.title= element_text(size=14,face="bold")
      )
    PCA_plot_backGround <- as.data.frame(PCbasis[, 1:2])
    if(is.null(test)){
      pblapply(1:length(res),function(i){
        p <- geom_point(data = as.data.frame(PCbasis[res[[i]],1:2]) ,aes(x = V1,y = V2,colour = "Phenotyping selected"),size =1,alpha = 1)
        g00 <- ggplot() + 
          geom_point(data = PCA_plot_backGround,aes(x = V1,y = V2,colour = "Background"),size =1,alpha = 1) + 
          labs(x = "PC1",y = "PC2",title= paste0("TP size = ",ntosel[i]))
        g01 <- g00 + p + scale_colour_manual(values = c("cadetblue","red"))+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "transparent",color = "black"),
                axis.text = element_text(size = 12),axis.title=element_text(size=14,face = "bold"),
                legend.position = "right",legend.text = element_text(size = 12),
                legend.title = element_blank(),legend.background = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 14,face = "bold")
          )
        pdf(paste0("TST_SP","_",ntosel[i],".pdf"),width = 12,height = 5)
        multiplot(g0,g01,cols = 2,layout = matrix(c(1,1,1,2,2,2,2,2),nrow = 1))
        dev.off()
        multiplot(g0,g01,cols = 2,layout = matrix(c(1,1,1,2,2,2,2,2),nrow = 1))
      })
    }else{
      pblapply(1:length(res),function(i){
        p <- geom_point(data = as.data.frame(PCbasis[res[[i]],1:2]) ,aes(x = V1,y = V2,colour = "Trainnig pop selected"),size =1,alpha = 1)
        eval(parse(text = paste0("p",i,"<-p")))
        g00 <- ggplot() + 
          geom_point(data = PCA_plot_backGround,aes(x = V1,y = V2,colour = "Background"),size =1,alpha = 1) + 
          labs(x = "PC1",y = "PC2",title= paste0("TP size = ",ntosel[i])) +
          geom_point(data = as.data.frame(PCbasis[test,1:2]),aes(x = V1,y = V2,colour = "Testing pop"),size =1,alpha = 1)
        
        g01 <- g00 + p + scale_colour_manual(values = c("cadetblue","blue","red"))+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "transparent",color = "black"),
                axis.text = element_text(size = 12),axis.title=element_text(size=14,face = "bold"),
                legend.position = "right",legend.text = element_text(size = 12),
                legend.title = element_blank(),legend.background = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 14,face = "bold")
          )
        pdf(paste0("TSR_",method,"_",ntosel[i],".pdf"),width = 12,height = 5)
        multiplot(g0,g01,cols = 2,layout = matrix(c(1,1,1,2,2,2,2,2),nrow = 1))
        dev.off()
        multiplot(g0,g01,cols = 2,layout = matrix(c(1,1,1,2,2,2,2,2),nrow = 1))
      })
    }
  }
  return(res)
}