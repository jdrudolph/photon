# Running Perseus through Docker

PHOTON can be run in `legacy` mode through a webserver that is packaged as a docker container,
as is described in the original publication.

The following document served mostly as an archive, as running PHOTON through Perseus is
recommended.

## Running PHOTON on Docker

PHOTON runs inside a 'container' and therefore requires docker to run across all platforms.
General information regarding the installation of docker can be
found on the docker [website](https://docs.docker.com/engine/installation/).


### Docker Toolbox (Windows, Mac)

First open the `Docker Quickstart Terminal`. After initialization (can take some time),
denote the IP address of docker (under the whale image).
Now you can run PHOTON by entering:

```bash
docker run -d -p 5000:5000 jdrudolph/photon
```

Now you can access PHOTON from your browser under the IP address
of docker followed by colon and 5000. For example `192.168.0.100:5000`.
Alternatively you can lookup the IP address using `docker-machine ip default`.

### Docker for Windows / Linux

Open a Powershell (terminal in linux) and run:

```bash
docker run -d -p 5000:5000 jdrudolph/photon
```

Now you can access PHOTON from your [browser](http://localhost:5000).

### Native Linux (Experts)

Consult the `Dockerfile` on how to run PHOTON on Linux natively.

## Stopping PHOTON

First list all running containers:

```bash
docker ps
```

Lookup the name of the container (last column) and stop it using

```bash
docker stop name_of_container
```

## Updating PHOTON

New versions of docker can be downloaded by running

```bash
docker pull jdrudolph/photon
```

before issuing the `run` command as shown above.


## Output description
A zip archive containing several output tables can be downloaded from the results page of PHOTON.

1. `subnet.csv`: edges of the reconstructed signaling network. This is a sub-network of the input network (which is located at `/db/anat/H_sapiens.net`). The `subnet.csv` file can be e.g. imported into Cytoscape for further analysis.
2. `scores.csv`: signaling functionality scores devised by PHOTON. The 'Significant' column allows for easy filtering of the table. The other columns provide more detail on each step of the algorithm, i.e. calculating an empirical score, obtaining p-values, and finally FDR-corrected q-values.
3. `go_scores.csv`: details all GO annotation enrichment results for the reconstructed network. Filtering for the 'rejected' column (=1) will yield all sign. enriched categories. The naming of the other columns follows the convention of the [statistical test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html) used.
4. `predictions.csv`: contains information regarding the functional phosphorylation site prediction. The '[proba](http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html#sklearn.linear_model.LogisticRegression.predict_proba)' column is the prediction score for each site. If the phosphorylation site was part of the training set the 'label' column will be 'TRUE'. One should take special care when the training set was small, i.e. few sites with 'label' set to TRUE are in the data set.

## My result graph contains only a sigle node '-1'
Sometimes PHOTON will not yield any results. Known technical sources of failure are:

1. The ANAT service is down. PHOTON relies on the [ANAT web server](http://www.cs.tau.ac.il/~bnet/ANAT/) which can be down for various reasons. Before reporting an error with PHOTON, please confirm that the Cytoscape plugin for ANAT is working as expected. Please note that you might still download the resulting activity scores derived by PHOTON.

## Common errors when running PHOTON

1. Uploading a data file which is not `ASCII` encoded might yield errors such as
    
        'utf-8' codec can't decode byte 0xa0 in position 0: invalid start byte
        
   Converting the file to `ASCII` will fix this issue.
2. There is a popup with a long error message similar to

        <ns2:confidence>0.6912981534952782</ns2:confidence>
        <ns2:fromNodeId>473</ns2:fromNodeId>
        <ns2:toNodeId>126961</ns2:toNodeId>
        </ns2:edgesData>
        <ns2:edgesData>
        ...
        
    This error messages is due to the [ANAT web server](http://www.cs.tau.ac.il/~bnet/ANAT/) not being available.
    Please confirm that the Cytoscape plugin for ANAT is working as expected before reporting an error for PHOTON.
    Please contact me or the authors of ANAT in case of errors.
    
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
 
