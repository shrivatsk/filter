FROM r-base:latest

RUN mkdir /app

COPY . /app

ADD Sift_PP.R /app/Sift_PP.R

WORKDIR /app

CMD ["Rscript", "/app/Sift_PP.R"]
