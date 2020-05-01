FROM ruby:2.6.5

## System dependencies

RUN apt-get update && apt install -y \
	rsync \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

## Install gems

COPY Gemfile /tmp
COPY Gemfile.lock /tmp
WORKDIR /tmp
## RUN bundle lock
RUN bundle install --jobs=4 && bundle clean --force

## Set up web server

EXPOSE 3000

## Create startup script

RUN echo '#! /bin/bash       \n\
    cd /bioconductor.org     \n\
    rake                     \n\
    cd output                \n\
    adsf\n' > .startup.sh     \\
    && chmod +x .startup.sh

CMD ["./.startup.sh"]
