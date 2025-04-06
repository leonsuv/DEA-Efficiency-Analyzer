# Use Python 3.13 as the base image
FROM python:3.13

# Set working directory
WORKDIR /app

# Copy requirements file first for caching
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy all project files to the container
COPY . .

# Make sure the application is accessible from outside the container
EXPOSE 8080

# Command to run the application
CMD ["python", "main.py"]
