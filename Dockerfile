FROM ruby:2.7.0

## bioconductor.org website requirements
RUN apt-get update && apt install -y \
	curl \
	libssl-dev \
	zlib1g-dev \
	autoconf \
	bison \
	build-essential \
	libyaml-dev \
	libreadline-dev \
	libncurses5-dev \
	libffi-dev \
	libgdbm-dev \
	libpq-dev \
	sqlite3 \
	libsqlite3-dev \
	rsync \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

## website
RUN mkdir /myapp
WORKDIR /myapp
COPY . /myapp

## Install gems
RUN bundle update --bundler

RUN rake

# Add a script to be executed every time the container starts.
EXPOSE 3000
WORKDIR /myapp/output

CMD ["adsf"]
