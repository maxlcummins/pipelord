FROM debian:stretch

ENV DEBIAN_FRONTEND noninteractive

### RUN set -ex; \

RUN apt-get update -qq; \
    apt-get install -y -qq git \
    apt-utils \
    wget \
    python3-pip \
    ncbi-blast+ \
    libz-dev \
    ; \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*;
    
ENV DEBIAN_FRONTEND Teletype

# Install python dependencies
RUN pip3 install -U biopython tabulate cgecore==1.3.2;

# Install kma 
RUN git clone --branch 1.0.1 --depth 1 https://bitbucket.org/genomicepidemiology/kma.git; \
    cd kma && make; \
    mv kma* /bin/

COPY pmlst.py /usr/src/pmlst.py 

RUN chmod 755 /usr/src/pmlst.py; 

ENV PATH $PATH:/usr/src
# Setup .bashrc file for convenience during debugging
RUN echo "alias ls='ls -h --color=tty'\n"\
"alias ll='ls -lrt'\n"\
"alias l='less'\n"\
"alias du='du -hP --max-depth=1'\n"\
"alias cwd='readlink -f .'\n"\
"PATH=$PATH\n">> ~/.bashrc

WORKDIR /workdir

# Execute program when running the container
ENTRYPOINT ["/usr/src/pmlst.py"]
