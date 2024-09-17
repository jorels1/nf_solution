# Use an official Python runtime as a parent image
FROM python:3.8-slim-buster

# Set the working directory in the container
WORKDIR /bin

# Copy the current directory contents into the container at /bin
COPY . /bin

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Run script.py when the container launches
CMD ["python"]
