FROM rocker/tidyverse
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
libnetcdf-dev \
libxml2-dev
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("plumber")'
RUN R -e 'remotes::install_github("chasemc/IDBacApp@1fd3186")'
RUN R -e 'remotes::install_github("chasemc/maldiblast")'
COPY maldiblast_*.tar.gz /app.tar.gz
RUN R -e 'remotes::install_local("/app.tar.gz")'
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');maldiblast::run_app()"
