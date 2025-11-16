FROM ubuntu:24.04

# Install system packages
RUN apt-get update && apt-get install -y \
    nano \
    curl \
    wget \
    htop \
    git \
    build-essential \
    libboost1.74-all-dev \
    liblapack-dev \
    csh \
    python3 \
    python3-dev \
    python3-venv \
    python3-pip \
    cmake \
    libx11-dev \
    libxext-dev \
    libxt-dev \
    libxpm-dev \
    libxft-dev \
    libxrandr-dev \
    libxinerama-dev \
    libomp-dev \
    libxi-dev \
    libomp-dev \
    zlib1g-dev \
    libxcursor-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Create Python venv
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade pip inside venv
RUN pip install --upgrade pip setuptools wheel

# Clone Scarabee from source
RUN git clone https://github.com/scarabee-dev/scarabee.git /src/scarabee

# OPTIONAL: checkout development branch
# RUN cd /src/scarabee && git switch develop

# Build and install Scarabee from source (inside venv)
RUN cd /src/scarabee && pip install -v .

# Install ENDFtk from source
RUN git clone https://github.com/njoy/ENDFtk.git /src/ENDFtk && \
    mkdir /src/ENDFtk/build && \
    cd /src/ENDFtk/build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=$(which python3) ../ && \
    make -j8

# Install Papillon from source
RUN git clone https://github.com/HunterBelanger/papillon-ndl.git /src/papillon-ndl && \
    cd /src/papillon-ndl && \
    cmake -E make_directory build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    cmake --build . --target install

# Add build paths to $PYTHONPATH
ENV PYTHONPATH="/src/ENDFtk/build/python:/src/ENDFtk/build/_deps/tools-build/python:/src/papillon-ndl/build:$PYTHONPATH"

# Install h5py inside venv
RUN pip install h5py

# Copy FRENDY install script into the image
COPY install_frendy.sh /tmp/install_frendy.sh

# Make FRENDY install script executable and run it
RUN chmod +x /tmp/install_frendy.sh && \
    /bin/bash /tmp/install_frendy.sh

# Create symbolic link for frendy executable
RUN ln -s /src/frendy/frendy_20241030/frendy/main/frendy.exe /usr/bin/frendy

# Download ENDF/B-VIII.0 library
RUN mkdir -pv /nuclear_data/endf && \
    cd /nuclear_data/endf && \
    wget https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0.zip && \
    unzip ENDF-B-VIII.0.zip 

CMD ["/bin/bash"]
