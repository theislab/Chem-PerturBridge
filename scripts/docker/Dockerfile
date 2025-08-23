FROM python:3.12

# Define the HOME environment variable
ENV HOME=/home/ubuntu

# Install Python packages
RUN pip install --no-cache-dir lamindb scanpy gcsfs decoupler pubchempy magnitude

# Set working directory
WORKDIR $HOME
