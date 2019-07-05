FROM bioconductor/devel_base2:latest

RUN mkdir /app

COPY . /app

ADD Sift_PP.R /app/Sift_PP.R

ADD dbnsfp_filter.R /app/dbnsfp_filter.R

WORKDIR /app

RUN R -e 'install.packages(c("magrittr", "tidyr"))'

CMD ["Rscript", "/app/dbnsfp_filter.R"]




