FROM rocker/tidyverse:3.4.0

MAINTAINER Jason Serviss <jason.serviss@ki.se>

# System dependencies for required R packages
RUN  rm -f /var/lib/dpkg/available \
  && rm -rf  /var/cache/apt/* \
  && apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    git

RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('Rtsne/0.13','ggthemes/3.4.0','liftr/0.7','igraph/1.0.1','magrittr/1.5','knitr/1.17','rmarkdown/1.6','testthat/1.0.2','multipanelfigure/0.9.0','gtable/0.2.0','ggraph/1.0.0','printr/0.1'))"

RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('topGO/2.30.0','biomaRt/2.34.0','GOSim/1.14.0','AnnotationDbi/1.38.1','GO.db/3.4.1'))"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_remotes.R');install_remotes(c('GranderLab/acidAdaptedRNAseq'))"

RUN mkdir /liftrroot/
WORKDIR /liftrroot/
