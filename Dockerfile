FROM ubuntu:bionic

## Update apt
## Source: https://linuxize.com/post/how-to-install-ruby-on-ubuntu-18-04/
## Install rbenv essentials
RUN apt-get update \
	&& apt install -y git \
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
	rsync

## Install rbenv and ruby-build
RUN git clone https://github.com/rbenv/rbenv.git /root/.rbenv
RUN git clone https://github.com/rbenv/ruby-build.git /root/.rbenv/plugins/ruby-build
RUN /root/.rbenv/plugins/ruby-build/install.sh
ENV PATH /root/.rbenv/bin:$PATH
RUN echo 'eval "$(rbenv init -)"' >> /etc/profile.d/rbenv.sh # or /etc/profile
RUN echo 'eval "$(rbenv init -)"' >> .bashrc

## Add $HOME/.rbenv/bin to the user PATH
#RUN echo 'export PATH="$HOME/.rbenv/bin:$PATH"' >> ~/.bashrc \
#    && echo 'eval "$(rbenv init -)"' >> ~/.bashrc \
#    && (source ~/.bashrc)

## Install latest version of ruby
RUN rbenv install 2.5.1 \
    && rbenv global 2.5.1

RUN gem update --system \
    && gem install bundler:1.16.1

## Clone repo
RUN git clone https://github.com/Bioconductor/bioconductor.org \
	&& cd bioconductor.org \
	&& bundle install \
	&& rake
