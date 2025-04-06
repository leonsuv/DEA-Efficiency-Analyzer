# Remove any existing container with the same name
docker rm -f dea-efficiency-analyzer-container 2>/dev/null || true

# Build the Docker image and run it
docker build -t dea-efficiency-analyzer . && \

# start the container and expose the port 8080
docker run -p 8080:8080 --name dea-efficiency-analyzer-container dea-efficiency-analyzer

echo "Server: http://localhost:8080"
