FROM rocker/tidyverse
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
libnetcdf-dev \
libxml2-dev
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("plumber")'
RUN R -e 'remotes::install_github("chasemc/IDBacApp@1fd3186")'
RUN R -e 'remotes::install_github("chasemc/maldiblast")'
EXPOSE 8000
CMD R -e 'plumber::plumb(file=system.file("api.R",package = "maldiblast"))$run(host="0.0.0.0", port=8000)'
