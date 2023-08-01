FROM openjdk:11

# Install dependencies
RUN apt-get update && \
    apt-get install -y cmake build-essential git && \
    rm -rf /var/lib/apt/lists/*

# Clone and install jnainchi
WORKDIR /usr/src
RUN git clone https://github.com/dan2097/jna-inchi.git && \
    mkdir -p jna-inchi/build && \
    cd jna-inchi/build && \
    cmake .. && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf jna-inchi

# Copy MetFragCommandLine-2.5.0.jar to the container
WORKDIR /usr/src/myapp
COPY . /usr/src/myapp
COPY MetFragCommandLine-2.5.0.jar .

# Set the entrypoint to run MetFrag
#ENTRYPOINT ["java", "-jar", "MetFragCommandLine-2.5.0.jar"]
