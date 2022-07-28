FROM alpine:3.14
RUN apk update

ENV SAMTOOLS_VERSION 1.3.1

RUN apk add --update pcre-dev openssl-dev zlib-dev \
 && apk add --virtual build-dependencies build-base curl \
 && curl -L -o samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    http://jaist.dl.sourceforge.net/project/samtools/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
 && cd samtools-${SAMTOOLS_VERSION}/ \
 && ./configure --without-curses \
 && make \
 && make install \
&& apk del build-dependencies \
&& rm -rf /var/chache/apk/*

RUN apk add python3
RUN apk add py3-pip
RUN pip install MAnorm2-utils
RUN apk add bash

RUN mkdir /data
