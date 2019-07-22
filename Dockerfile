FROM bioconductor/devel_base2:latest

RUN mkdir /app

COPY . /app

ADD Sift_PP.R /app/Sift_PP.R

ADD dbnsfp_filter.R /app/dbnsfp_filter.R 

ADD annovarfilter.R /app/annovarfilter.R

WORKDIR /app

RUN R -e 'install.packages("tidyr")'

CMD ["Rscript", "/app/annovarfilter.R"]




