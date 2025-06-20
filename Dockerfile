FROM ubuntu:22.04

# Install required packages
RUN apt-get update && \
    apt-get install -y cmake g++ git make && \
    apt-get clean

# Set working directory
WORKDIR /app

# Copy project into container
COPY . .

# Configure and build
RUN cmake -S . -B build && cmake --build build

# Optional: set default command to bash (interactive shell)
CMD ["/bin/bash"]
