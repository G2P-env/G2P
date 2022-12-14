# G2P
<a href="https://www.r-project.org/" target="_blank"><img src="https://img.shields.io/badge/language-R-orange?style=plastic"></a>
<a href="https://cran.r-project.org/bin/windows/base/old/" target="_blank"><img src="https://img.shields.io/badge/R%20version-%3E%3D%203.6.0-orange?style=plastic"></a>
<a href="https://govs-pack.github.io/" target="_blank"><img src="https://img.shields.io/badge/webpage-ready-green?style=plastic"></a>
<a href="https://sylabs.io/" target="_blank"><img src="https://img.shields.io/badge/Singularity-%3E%3D3.1-orange?style=plastic"></a>
![](https://img.shields.io/badge/platform-Win%20%7C%20Linux%20%7C%20MacOS-lightgrey?style=plastic)<br/>

## Overview
__G2P__ (**G**enotype **to** **P**henotype) is an integrative environment in form of Singularity container, which not only contains a library of 16 state-of-the-art GS models and 13 evaluation metrics for model evaluation and selection but also provides several stand-alone software and easy-to-use functions for data preprocessing, population analysis, integration of prediction results from multiple models, and refinement of training datasets. G2P provides a comprehensive environment for genomic selection to facilitate the comparison and selection of appropriate model. G2P not only was developed and packaged in a Singularity container, allowing easy installation and upgrade, but also was accessible as an stand-alone R package.

<div align="center">
<img src="https://g2p-env.github.io/img/overall.png" width="900"/>
</div>
<div align="center">
  Overview of G2P
</div>

## Installation 
## Installation of G2P container 

#### 1. Installation of Singularity (Linux, assuming ubuntu)
#### 1.1 Install system dependencies
```bash
$ sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup
```
#### 1.2 Installing GO

```bash
$ export VERSION=1.16.4 OS=linux ARCH=amd64 && \  # Replace the values as needed
  wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \ # Downloads the required Go package
  sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \ # Extracts the archive
  rm go$VERSION.$OS-$ARCH.tar.gz    # Deletes the ``tar`` file
```

![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) __Note__: you can vist [Go Downloads page](https://go.dev/dl/) for suitable to the environment you are in
```

#### 1.3 Installing Singularity
You can download SingularityCE from one of the releases. To see a full list, visit the [GitHub release page](https://github.com/sylabs/singularity/releases). After deciding on a release to install, you can run the following commands to proceed with the installation (here we use version 3.8.1 as example).

```bash
$ export VERSION=3.8.1 && # adjust this as necessary \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
    tar -xzf singularity-ce-${VERSION}.tar.gz && \
    cd singularity-ce-${VERSION}
```
#### 1.4 Compile singularity
```bash 
$ ./mconfig && \
    make -C builddir && \
    sudo make -C builddir install
```
![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) __Note__: If you get an error when compile, please try the following command or refer to [Singularity websit for documents](https://sylabs.io/docs/). 

```bash
## Motify the file mlocal/frags/go_common_opts.mk in singularity path???find the following line and change it
GOPROXY := https://goproxy.cn,direct
## At the same time??? motify the global proxy
go env -w GOPROXY=https://goproxy.cn,direct
```

#### 1.5 Installing Singularity for other operating system
Singularity can be installed via Vagrant Boxes for Windows and Intel cores MacOS, the details please refer to [Singularity installation guide](https://docs.sylabs.io/guides/3.8/admin-guide/installation.html#installation-on-windows-or-mac).

#### 2. Pulling G2P container
```bash
singularity pull G2P.sif library://mym89757/repo/g2p:latest
```
## Installation of G2P R package 
![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) __Note__: the functions, like _GSFiltering_, _GSImputation_, _GSTransForm_, _GSRead_ shown in figure ["Overview of G2P" a](https://github.com/G2P-env/G2P/blob/main/figures/Fig%201.png) needs the opearation evironment and stand-alone software on the circle, namely??? these function in G2P pcakage can't be performed normally unless the dependencis shown in figure ["Overview of G2P" a](https://github.com/G2P-env/G2P/blob/main/figures/Fig%201.png)(upper semicircle) have been correctly installed and configuration. G2P package can be installed for MacOS, Linux and windows and the functions shown in figure ["Overview of G2P" b](https://github.com/G2P-env/G2P/blob/main/figures/Fig%201.png) can run normally.
#### 1.  Github install
```R
## install dependencies and GOVS
install.packages(c("ggplot2","brnn","glmnet","spls","pls","e1071","BGLR","rrBLUP","randomForest","hglm","hglm.data","parallel","pROC","PRROC","STPGA","reshape","reshape2","grid","pbapply","pheatmap"))

require("devtools")
install_github("G2P-env/G2P") 
```
#### 2.  Download [.tar.gz package](https://github.com/G2P-env/G2P/raw/main/G2P_1.0.tar.gz) and install <br/>
```R
## install dependencies and GOVS with bult-in vignette
install.packages(c("ggplot2","brnn","glmnet","spls","pls","e1071","BGLR","rrBLUP","randomForest","hglm","hglm.data","parallel","pROC","PRROC","STPGA","reshape","reshape2","grid","pbapply","pheatmap"))

install.packages("DownloadPath/GOVS_1.0.tar.gz")
```
## Links
* G2P Homepage: https://g2p-env.github.io/
* QuickStart: https://g2p-env.github.io/QuickStart/
* Tutorial: https://g2p-env.github.io/Tutorial/
* Singularity container download: https://cloud.sylabs.io/library/mym89757/repo/g2p
* R package download: https://github.com/G2P-env/G2P/raw/main/G2P_1.0.tar.gz
