FROM rocker/verse:4.3.3

RUN apt-get update -qq
RUN apt-get -y install cmake
RUN apt-get -y --no-install-recommends install libglu1-mesa-dev
RUN apt-get -y install libxml2-dev
RUN apt-get -y install ghostscript
RUN apt-get -y install libxt-dev
RUN apt-get -y install git
RUN apt-get -y install -y python3
RUN apt-get -y install python3-pip
RUN apt-get -y install python3-venv

RUN install2.r –error \
--deps TRUE \
quarto \
reticulate \
knitr \
osfr \
bookdown \
grImport \
tiff \
rstatix \
ez \
readr \
BayesFactor \
bmp \
pracma \
eyelinker