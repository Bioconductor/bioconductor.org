FROM ruby:2.6.5

## System dependencies
RUN apt-get update && apt-get install -y \
    rsync \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

## Install gems
COPY Gemfile /opt
COPY Gemfile.lock /opt

## Set the working directory to /tmp/bioconductor.org
WORKDIR /opt/bioconductor.org

## Copy the remaining files from bioconductor.org/ to /opt/bioconductor.org
COPY . .

## Install bundle and dependencies
RUN gem install bundler
RUN bundle install --jobs=4 && bundle clean --force

## Set up web server
EXPOSE 3000

## Create startup script
RUN echo '#! /bin/bash' > .startup.sh \
    && echo 'cd /opt/bioconductor.org' >> .startup.sh \
    && echo 'rake' >> .startup.sh \
    && echo 'cd output' >> .startup.sh \
    && echo 'adsf' >> .startup.sh \
    && chmod +x .startup.sh

CMD ["/opt/bioconductor.org/.startup.sh"]
