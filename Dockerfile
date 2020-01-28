FROM ruby:2.5.1

## bioconductor.org website requirements
RUN apt-get update && apt install -y \
	git \
	curl \
	libssl-dev \
	libreadline-dev \
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

## Install gems
COPY Gemfile /myapp/Gemfile
COPY Gemfile.lock /myapp/Gemfile.lock
RUN bundle install

## Copy website into image
COPY . /myapp

RUN rake

# Add a script to be executed every time the container starts.
EXPOSE 3000
WORKDIR /myapp/output

CMD ["adsf"]
