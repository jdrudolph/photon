# Installation

PHOTON runs inside a 'container' and therefore requires docker to run across all platforms.
General information regarding the installation of docker can be
found on the docker [website](https://docs.docker.com/engine/installation/).

# Running PHOTON

## Docker Toolbox (Windows, Mac)

First open the `Docker Quickstart Terminal`. After initialization (can take some time),
denote the IP address of docker (under the whale image).
Now you can run PHOTON by entering:

```bash
docker run -d -p 5000:5000 jdrudolph/photon
```

Now you can access PHOTON from your browser under the IP address
of docker followed by colon and 5000. For example `192.168.0.100:5000`.
Alternatively you can lookup the IP address using `docker-machine ip default`.

## Docker for Windows / Linux

Open a Powershell (terminal in linux) and run:

```bash
docker run -d -p 5000:5000 jdrudolph/photon
```

Now you can access PHOTON from your [browser](http://localhost:5000).

## Native Linux (Experts)

Consult the `Dockerfile` on how to run PHOTON on Linux natively.

# Stopping PHOTON

First list all running containers:

```bash
docker ps
```

Lookup the name of the container (last column) and stop it using

```bash
docker stop name_of_container
```

# Updating PHOTON

New versions of docker can be downloaded by running

```bash
docker pull jdrudolph/photon
```

before issuing the `run` command as shown above.

# Troubleshooting
Please open an issue [here](https://github.com/jdrudolph/photon/issues), or
contact me directly if you have issues with installing or using PHOTON.

# Licence
If the current licencing (see `LICENCE.txt`) does not suit your needs,
please contacte me for alternative licencing options.

# PHOTON data format example

    GeneID,Amino.Acid,Position,avg,Symbol
    5097,S,962,-0.47417884405658556,PCDH1
    5097,S,984,-0.6438018715183813,PCDH1
    4300,S,522,-0.5172999908818245,MLLT3
    81628,S,264,0.8682687511988673,TSC22D4
    81628,S,279,2.684096205546553,TSC22D4
    4034,S,25,-0.593695490783283,LRCH4
    4034,S,380,1.6585730683114723,LRCH4
    6651,S,1822,0.046145203537912814,SON
    6651,S,1737,-0.06153405253207198,SON
 
