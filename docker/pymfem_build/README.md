# pymfem-docker

## Getting started
Dockerfile to installation of PyMFEM

### Docker build
Put `Dockerfile` and `script.sh` in the same folder and execute
```
docker build -t pymfem .
```

### Run pymfem (delete container after exit)
```
# start pymfem without gpu but with jupyter notebook

docker run -it --name pymfem_nb -p 8888:8888 -v D:\docker_ws:/docker_ws pymfem:snapshot bash

# start with gpu (nvidia-smi will be available)
docker run -it --gpus all --name gpu_test -p 8888:8888 -v D:\docker_ws:/docker_ws pymfem:snapshot bash

docker stop container_name
docker rm -f container_name
docker start container_name
docker exec -it container_id bash
```