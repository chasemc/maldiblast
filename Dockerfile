FROM trestletech/plumber
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_github("r-lib/remotes")'
RUN R -e 'remotes::install_github("chasemc/IDBacApp@1fd3186")'
COPY maldiblast_*.tar.gz /app.tar.gz
RUN R -e 'remotes::install_local("/app.tar.gz")'
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');maldiblast::run_app()"
