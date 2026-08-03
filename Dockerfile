FROM ubuntu:24.04

# Ruby 3.2 from Ubuntu 24.04 LTS (supported until 2029). The previous image was
# ruby:2.6.5 on Debian Buster, both long EOL — Buster was archived, so that build
# needed apt sources rewritten to archive.debian.org just to install anything.

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    BUNDLE_PATH=/usr/local/bundle

RUN apt-get update && apt-get install -y --no-install-recommends \
      ruby-full ruby-dev \
      build-essential pkg-config \
      libxml2-dev libxslt1-dev zlib1g-dev libffi-dev libyaml-dev \
      libsqlite3-dev libpq-dev \
      git rsync ca-certificates curl \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/bioconductor.org

# Gems first, so dependency installation caches independently of content changes.
COPY Gemfile Gemfile.lock ./
RUN gem install bundler --no-document \
    && bundle install --jobs=4 \
    && bundle clean --force

COPY . .

EXPOSE 3000

# `rake` builds the site into output/; adsf serves it. Same contract as before.
CMD ["bash", "-lc", "bundle exec rake && cd output && bundle exec adsf -p 3000"]
