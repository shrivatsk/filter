FROM r-base:latest

RUN mkdir /app

COPY . /app

ADD Sift_PP.R /app/Sift_PP.R

ADD dbnsfp_filter.R /app/dbnsfp_filter.R

WORKDIR /app

RUN Rscript -e "install.packages('BiocManager');\
    devtools::install_github('sbg/sevenbridges-r', repos = BiocManager::repositories(), build_vignettes = FALSE, dependencies = TRUE)"

CMD ["Rscript", "/app/dbnsfp_filter.R"]
