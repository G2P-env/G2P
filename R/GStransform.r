############################### big data readin and format transformation ################################
#' @export GSTransForm
GSTransForm <- function(pattern, input, output, workdir = "./", memory = 8 ){
  oriPath <- getwd()
  setwd(workdir)
  checkPattern <-  pattern %in% c("hmp2plink","plink2vcf","vcf2hmp","vcf2plink","plink2hmp","hmp2vcf")
  
  oriFilePath <- input
  fileName <- output
  if (! checkPattern ){
    stop("Error: not defined parttern for genotype transformation!")
  }
  suppressMessages(switch(pattern,
         hmp2plink = {system(paste0("run_pipeline.pl -Xmx",memory,"g -SortGenotypeFilePlugin -inputFile ",oriFilePath," -outputFile ",fileName,"_sort -fileType Hapmap"));
           system(paste0("run_pipeline.pl -Xmx",memory,"g -fork1 -h ",fileName,"_sort.hmp.txt -export -exportType VCF"));
           system(paste0("plink --allow-extra-chr --keep-allele-order --vcf ",fileName,"_sort.vcf --out ",fileName," --double-id --biallelic-only"));
           system(paste0("rm ",fileName,"_sort.hmp.txt ",fileName,"_sort.vcf"))},
         plink2vcf = system(paste0("plink --allow-extra-chr --keep-allele-order --bfile ",oriFilePath," --recode vcf-iid --out ",fileName)),
         vcf2hmp = system(paste0("run_pipeline.pl -Xmx",memory,"g -importGuess ",oriFilePath," -export ",fileName," -exportType HapmapDiploid")),
         vcf2plink = system(paste0("plink --allow-extra-chr --keep-allele-order --vcf ",oriFilePath," --out ",fileName," --double-id --biallelic-only")),
         plink2hmp = {system(paste0("plink --allow-extra-chr --keep-allele-order --bfile ",oriFilePath," --recode vcf-iid --out ",fileName));
           system(paste0("run_pipeline.pl -Xmx",memory,"g -importGuess ",fileName,".vcf -export ",fileName," -exportType HapmapDiploid"));
           system(paste0("rm ",fileName,".vcf ",fileName,".nosex"))},
         hmp2vcf = {
           system(paste0("run_pipeline.pl -Xmx",memory,"g -SortGenotypeFilePlugin -inputFile ",oriFilePath," -outputFile ",fileName,"_sort -fileType Hapmap"));
           system(paste0("run_pipeline.pl -Xmx",memory,"g -fork1 -h ",fileName,"_sort.hmp.txt -export ",fileName," -exportType VCF"));
           system(paste0("rm ",fileName,"_sort.hmp.txt "))
         }
                     
  ))
  setwd(oriPath)
}

## data filtering ###############
#' @export GSFiltering 
GSFiltering <- function(input, output, inputformat, missingRate = 0.2, maf = 0.05, window = 100, step =5, r2 = 0.1, workdir = "./", memory = 8){
  oriPath <- getwd()
  setwd(workdir)
  checkPattern <-  inputformat %in% c("hmp","vcf","plink")
  
  oriFilePath <- input
  fileName <- output
  if (! checkPattern ){
    stop("Error: not defined parttern for genotype transformation!")
  }
  
  if (inputformat %in% c("hmp","vcf")) {
    if (inputformat == "hmp"){
      GSTransForm("hmp2plink",input, output, workdir = workdir, memory = memory )
    } else if (inputformat == "vcf"){
      GSTransForm("vcf2plink",input, output, workdir = workdir, memory = memory )
    }
    system(paste0("plink --bfile ",output," --geno ",missingRate," --maf ",maf," --out ",output," --write-snplist -indep-pairwise ",window," ",step," ",r2))
  } else {
    system(paste0("plink --bfile ",input," --geno ",missingRate," --maf ",maf," --out ",output," --write-snplist -indep-pairwise ",window," ",step," ",r2))
  }
  
  ## trans and extract and trans
  system(paste0("plink --bfile ",output," --extract ",output,".prune.in --out ",output,"_prune --make-bed"))
  rmfile <- paste0(output,"_prune")
  
  GSTransForm("plink2hmp",paste0(output,"_prune"), output, workdir = workdir, memory = memory )
  system(paste0("rm ",rmfile,".bed ",rmfile,".nosex ",rmfile,".bim ",rmfile,".fam "))
  system(paste0("rm ",output,".bed ",output,".snplist ",output,".bim ",output,".fam "))
}

## GS imputation ##
#' @export GSImputation
GSImputation <- function(input, output, inputformat, workdir = "./", memory = 8,threads = 1,window = 40,overlap = 2){
  oriPath <- getwd()
  setwd(workdir)
  checkPattern <-  inputformat %in% c("hmp","vcf","plink")
  
  # oriFilePath <- input
  # fileName <- output
  if (! checkPattern ){
    stop("Error: Unsupported genotypic data type!")
  }
  
  suppressMessages(switch(inputformat,
                          hmp = {GSTransForm(pattern="hmp2vcf",input=input,output=paste0(output,"_tmp"),workdir=workdir,memory = memory);
                                 system(paste0("java -Xmx",memory,"G -jar /opt/beagle.jar gt=",output,"_tmp.vcf"," out=",output," nthreads=",threads," window=",window," overlap=",overlap));
                                 system(paste0("gzip -d ",output,".vcf.gz"))
                            },
                          plink = {GSTransForm(pattern="plink2vcf",input=input,output=paste0(output,"_tmp"),workdir=workdir,memory = memory);
                                   system(paste0("java -Xmx",memory,"G -jar /opt/beagle.jar gt=",output,"_tmp.vcf"," out=",output," nthreads=",threads," window=",window," overlap=",overlap));
                                   system(paste0("gzip -d ",output,".vcf.gz"))},
                          vcf = {system(paste0("java -Xmx",memory,"G -jar /opt/beagle.jar gt=",input," out=",output," nthreads=",threads," window=",window," overlap=",overlap));
                                 system(paste0("gzip -d ",output,".vcf.gz"))}
  ))
  setwd(oriPath)
}

## GS read from file ###############
#' @export GSRead
GSRead <- function(input, output, inputformat, workdir = "./", memory = 8){
  require("data.table")
  oriPath <- getwd()
  setwd(workdir)
  checkPattern <-  inputformat %in% c("hmp","vcf","plink")
  
  oriFilePath <- input
  fileName <- output
  if (! checkPattern ){
    stop("Error: Unsupported genotypic data type!")
  }
  
  suppressMessages(switch(inputformat,
                          hmp = {GSTransForm(pattern="hmp2plink",input=input,output=paste0(output,"_tmp"),workdir=workdir,memory = memory);
                            system(paste0("plink --bfile ",output,"_tmp --recode A --output-missing-genotype - --out ",output));
                            G <- as.matrix(fread(paste0(output,".raw")));
                            rownames(G) <- G[,1]
                            G <- G[,-c(1:6)]
                            suppressWarnings(class(G) <- "numeric") 
                          },
                          vcf = {GSTransForm(pattern="vcf2plink",input=input,output=paste0(output,"_tmp"),workdir=workdir,memory = memory);
                            system(paste0("plink --bfile ",output,"_tmp --recode A --output-missing-genotype - --out ",output));
                            G <- as.matrix(fread(paste0(output,".raw")));
                            rownames(G) <- G[,1]
                            G <- G[,-c(1:6)]
                            suppressWarnings(class(G) <- "numeric") 
                          },
                          plink = {system(paste0("plink --bfile ",input," --recode A --output-missing-genotype - --out ",output));
                            G <- as.matrix(fread(paste0(output,".raw")));
                            rownames(G) <- G[,1]
                            G <- G[,-c(1:6)]
                            suppressWarnings(class(G) <- "numeric") }
  ))
  setwd(oriPath)
  G
}