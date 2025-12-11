# Nuclear data library generation using containers

To run, install Docker and Docker Compose (or some alternative container runtime),
then execute:

```shell
docker compose up --build -d
docker exec -it scarabee-ubuntu /bin/bash
```

The first command takes a while to complete (in the order of tens of minutes). Once in the container, run:

```shell
cd /data
python3 endf8.py
```

The second command takes a while to complete (in the order of days)

If a crash occurs, it might be due to RAM limitations. 
Try increasing the RAM available to container runtime.