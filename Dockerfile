FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="ShiyuanChen <schen647@ucr.edu>"
RUN apt update -y && apt intall git python3 python3-pip -y

# cd to /opt and wget https://nodejs.org/dist/v18.15.0/node-v18.15.0-linux-x64.tar.xz
# tar -zxvf node-v18.15.0-linux-x64.tar.xz
# ln -s /opt/node-v18.15.0-linux-x64/bin/node /usr/bin/
# ln -s /opt/node-v18.15.0-linux-x64/bin/npm /usr/bin/
RUN apt install wget -y && cd /opt && wget https://nodejs.org/dist/v18.15.0/node-v18.15.0-linux-x64.tar.xz && tar -zxvf node-v18.15.0-linux-x64.tar.xz && ln -s /opt/node-v18.15.0-linux-x64/bin/node /usr/bin/ && ln -s /opt/node-v18.15.0-linux-x64/bin/npm /usr/bin/

#clone https://github.com/schen647/shiba2Motif into /opt and install requirements with pip3 and install packages with npm
RUN cd /opt && git clone https://github.com/schen647/shiba2Motif && cd shiba2Motif && pip3 install -r requirements.txt && npm install && cd /opt


#wget https://meme-suite.org/meme/meme-software/5.5.0/meme-5.5.0.tar.gz
#tar -zxvf *.gz
#cd meme*
#apt install g++ make perl ghostscript zlib1g-dev libxml-simple-perl && ./configure --prefix=/opt/meme  --enable-build-libxml2 --enable-build-libxslt && make && make test && make install && ln -s /opt/meme/bin/* /usr/bin/

RUN cd /opt && wget https://meme-suite.org/meme/meme-software/5.5.0/meme-5.5.0.tar.gz && tar -zxvf *.gz && cd meme* && apt install g++ make perl ghostscript zlib1g-dev libxml-simple-perl && ./configure --prefix=/opt/meme  --enable-build-libxml2 --enable-build-libxslt && make && make test && make install && ln -s /opt/meme/bin/* /usr/bin/

WORKDIR /opt/shiba2Motif
RUN npm run host &
CMD ["bash"]
