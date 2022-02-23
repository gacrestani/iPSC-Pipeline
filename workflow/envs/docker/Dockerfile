FROM rocker/r-base

# Adds R packages install script.
RUN mkdir Code
RUN mkdir Software

ADD install_packages.R Code/install_packages.R
ADD fusca.zip Software/fusca.zip

# Install required system dependencies.
RUN apt-get update -qq && apt-get -y install \
libxml2-dev \
libxt-dev \
libcairo2-dev \
libcairo2-doc \
libsqlite3-dev \
libmariadbd-dev \
libpq-dev \
libssh2-1-dev \
unixodbc-dev \
libcurl4-openssl-dev \
libssl-dev \
libharfbuzz-dev \
libfribidi-dev \
libfreetype6-dev \
libpng-dev \
libtiff5-dev \
libjpeg-dev \
default-jre \
default-jdk \
cmake

# Reconfigure java to install rJava
RUN R CMD javareconf

# Install required R pakcages. This may take a while.
RUN Rscript "/Code/install_packages.R"
