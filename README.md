# Instructions:
1. Clone or Download the Repository
```bash
git clone https://github.com/arjunpatel96/StressME-DynamicME
```
2. From the root of the cloned/unzipped folder (where the Dockerfile is), build the docker image by running: 
```bash
docker build -t stressme_with_dynamicme .
```
3. Still in the same folder/terminal, start the container with:
```bash
docker run -it -p 5000:5000 -v "$PWD":/app stressme_with_dynamicme
```
4. Open a browser and go to http://localhost:5000