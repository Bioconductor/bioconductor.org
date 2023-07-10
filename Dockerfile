FROM ruby:2.6.5

## System dependencies
RUN apt-get update && apt-get install -y \
    rsync \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

## Install gems
COPY Gemfile /tmp
COPY Gemfile.lock /tmp

## Set the working directory to /tmp/bioconductor.org
WORKDIR /tmp/bioconductor.org

## Copy the remaining files from bioconductor.org/ to /tmp/bioconductor.org
COPY . .

## Install bundle and dependencies
RUN gem install bundler
RUN bundle install --jobs=4 && bundle clean --force

## Set up web server
EXPOSE 3000

## Create startup script
RUN echo '#! /bin/bash' > .startup.sh \
    && echo 'rake' >> .startup.sh \
    && echo 'cd output' >> .startup.sh \
    && echo 'adsf' >> .startup.sh \
    && chmod +x .startup.sh

CMD ["./.startup.sh"]
