FROM ruby:2.6.5

## bioconductor.org website requirements
RUN apt-get update && apt install -y \
	rsync \
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
