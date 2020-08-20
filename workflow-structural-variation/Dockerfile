FROM modulator:latest

MAINTAINER Fenglin Chen <f73chen@uwaterloo.ca>

# packages should already be set up in modulator:latest
USER root

# move in the yaml to build modulefiles from
COPY recipes/delly_recipe.yaml /modulator/code/gsi/recipe.yaml

# build the modules and set folder & file permissions
RUN ./build-local-code /modulator/code/gsi/recipe.yaml --initsh /usr/share/modules/init/sh --output /modules && \
	find /modules -type d -exec chmod 777 {} \; && \
	find /modules -type f -exec chmod 777 {} \;

# install required packages
RUN apt-get -m update && apt-get install -y gzip zip unzip

# add the user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 ubuntu
USER ubuntu

# copy the setup file to load the modules at startup
COPY .bashrc /home/ubuntu/.bashrc

# set environment paths for modules
#ENV PICARD_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/picard-2.19.2"
#ENV JAVA_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/java-8"
#ENV DELLY_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/delly-0.8.1"
#ENV BCFTOOLS_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/bcftools-1.9"
#ENV HTSLIB_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9"
#ENV TABIX_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/tabix-0.2.6"
#ENV VCFTOOLS_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/vcftools-0.1.16"
#ENV PERL_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/perl-5.30"

#ENV PATH="/modules/gsi/modulator/sw/Ubuntu18.04/vcftools-0.1.16/bin:/modules/gsi/modulator/sw/Ubuntu18.04/perl-5.30/bin:/modules/gsi/modulator/sw/Ubuntu18.04/tabix-0.2.6/bin:/modules/gsi/modulator/sw/Ubuntu18.04/bcftools-1.9/bin:/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/bin:/modules/gsi/modulator/sw/Ubuntu18.04/delly-0.8.1/bin:/modules/gsi/modulator/sw/Ubuntu18.04/java-8/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
#ENV MANPATH="/modules/gsi/modulator/sw/Ubuntu18.04/vcftools-0.1.16/share/man:/modules/gsi/modulator/sw/Ubuntu18.04/perl-5.30/share/man:/modules/gsi/modulator/sw/Ubuntu18.04/bcftools-1.9/share/man:/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/share/man:/modules/gsi/modulator/sw/Ubuntu18.04/java-8/man"
#ENV LD_LIBRARY_PATH="/modules/gsi/modulator/sw/Ubuntu18.04/vcftools-0.1.16/lib:/modules/gsi/modulator/sw/Ubuntu18.04/perl-5.30/lib:/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/lib:/modules/gsi/modulator/sw/Ubuntu18.04/java-8/li"
#ENV PERL5LIB="/modules/gsi/modulator/sw/Ubuntu18.04/vcftools-0.1.16/lib/site_perl:/modules/gsi/modulator/sw/Ubuntu18.04/perl-5.30/lib/site_perl"
#ENV PKG_CONFIG_PATH="/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/lib/pkgconfig"
#ENV LD_RUN_PATH="/modules/gsi/modulator/sw/Ubuntu18.04/bcftools-1.9/libexec"

CMD /bin/bash
