FROM bioconductor/release_base2:R3.4.0_Bioc3.5

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
    git \
    libjpeg-dev \
    libpng-dev \
    librsvg2-dev \
    libtiff5-dev \
    libudunits2-dev

RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R');install_cran(c('Rtsne/0.13','ggthemes/3.4.0','liftr/0.7','igraph/1.0.1','magrittr/1.5','knitr/1.17','rmarkdown/1.6','testthat/1.0.2','multipanelfigure/0.9.0','gtable/0.2.0','ggraph/1.0.0','printr/0.1','tidyr/0.7.2'))"

RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/biocLite.R');biocLite(c('topGO','biomaRt','GOSim','AnnotationDbi','GO.db', 'SummarizedExperiment'))"

#RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_remotes.R');install_remotes(c('GranderLab/acidAdaptedRNAseq'))"
RUN Rscript -e "devtools::install_github('GranderLab/acidAdaptedRNAseq', dep = FALSE)"

RUN mkdir /liftrroot/
WORKDIR /liftrroot/
