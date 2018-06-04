# Signaling functionality scores from phosphoproteomic data (PHOTON)

Make sure to follow the [installation instructions](docs/installation.md)

## Obtaining a high-confidence mouse protein-protein interaction network from STRING
1. We start by loading the STRING network table.

    **Matrix => Load => Raw upload**: Choose file `10090.protein.links.v10.5.txt.gz` and
    separating into columns by 'space'.

2. Next we want to filter out low-confidence interactions.

    First we have to change to column type of the score column.
    **Matrix => Rearrange => Change column types**: Convert `comined_score` from 'Text' to 'Numeric'.

3. We can now filter the edges and observe the changes in number of edges
    and score distribution by looking at the histogram of the scores before and after filtering.

    **Matrix => Filter rows => Filter rows based on numerical/main column**: Set `x >= 900`.

    **Matrix => Visualization => Histogram**: Select `combined_score` and set the y-axis to log-scale.

4. As a matter of convenience we can also scale the score to [0, 1] and rename the column to confidence.

    **Matrix => Basic => Transform**: Select `combined_score` and `x / 1000`.

    **Matrix => Rearrange => Rename columns**: Change `combined_score` to `Confidence`.

5. In order to improve compatability of the network with UniProt annotations we will
    change the format of the protein identifiers and remove the species prefix.

    **Matrix => Rearrange => Process text column**: Select both columns and set 'Regular expression' to `^10090\.(.*)`.

6. The resulting edge table provides the template for our network. We will save this
    intermediate result to file for use in the next step.

    **Matrix => Export => Generic matrix export**: Save as `10090.protein.links.v10.5_high_confidence.txt`.

## Create the PPI network from the edge table and prepare it for data annotation

1. We start by loading the data table from the previous step.

    **Matrix => Load => Generic matrix upload**: Choose file `10090.protein.links.v10.5_high_confidence.txt`.

2. Now we can transform the matrix into the dedicated 'Network collection' data structure.

    **Network => From matrix => Basic => From matrix**: Choose `protein1` and `protein2` as 'Source' and 'Target'.

3. Take your time to browse the different tables of the 'Network collection'.

4. We can explore the topology of the network by calculating node degrees (number of neighbors of any given node).

    **Network => Topology => Node degrees**

5. The node degrees were added to the node table of the network. It is easiest to plot them
    using any of the matrix visualizations in Perseus. We can easily convert network tables
    to matrices.

    **Network => To matrix => Basic => To matrix**: Select 'Nodes'.

6. When we plot the histogram of node distributions we can see that there are a number nodes with
    very large degrees.

    **Matrix => Analysis => Visualization => Histogram**: Select `Degree` and change the y-axis to log-scale.

7. Finding out the identity of the proteins by browsing ensemble protein identifiers is tedious.
    Instead, we would like to map them to gene names. We can't do so directy, but first have to map
    ENSP to UniProt identifiers.

    **Matrix => Processing => Annot. columns => To base identifiers**: Select `mainAnnot.mus_musculus.txt.gz`
    and 'Identifier type' as `ENSP`.

    **Matrix => Processing => Annot. columns => Add annotation**: Select `mainAnnot.mus_musculus.txt.gz`,
    set 'UniProt column' to `UniProt` and add the `Gene name` and `GOBP name` annotation.

8. Sorting by degree will show you the highest connected proteins in the network. Due to their high connectivity,
    These proteins can have a dominating effect on network analyses. Therefore, unless these proteins
    are of interest, removing them from the network collection can be helpful.

    **Network => Processing => Filter nodes => Filter nodes by numerical column**: Enter `x <= 2000`.

## Prepare the phospho-data and map it to the network

1. Upload the provided file.

    **Matrix => Load => Generic matrix upload**: Select `humphrey_phosphoSTY.txt`.

2. Before mapping to the network we need to make sure the data and network are defined
    in terms of the same identifiers. Here we need to map UniProt to ENSP.

    **Matrix => Processin => Annot. columns => Add annotation**: Select `mainAnnot.mus_musculus.txt.gz`
    and annotate with `ENSP`.

3. Now we can finally map the data to the network. And browse the network and the data simultaneously.

    **Network => Merge with matrix => Annotate => Annotate nodes**: Match `Node` with `ENSP` and copy
    all main columns. Copied main value should be kept separate.

## Perform the analysis and visualize the results

1. Perform PHOTON analysis to obtain 'Signaling functionality scores' for all proteins in the network.

    **Network => Modifications => PHOTON**: Leave default parameters.

2. For convenience the result is listed in different data nodes. The big network collection contains
    all results of the analysis. The second network collection, empty is this case, contains reconstructed
    signaling networks, if the box is checked in the previous case. The matrix contains all signaling
    functionality scores which can be analyzed further, e.g. by clustering

    **Matrix => Processing => Annot. columns => To base identifiers**: Select `mainAnnot.mus_musculus.txt.gz`
    and 'Identifier type' as `ENSP`.

    **Matrix => Processing => Annot. columns => Add annotation**: Select `mainAnnot.mus_musculus.txt.gz`,
    set 'UniProt column' to `UniProt` and add the `Gene name` and `GOBP name` annotation.

    **Matrix => Analysis => Hierarchical clustering**: The time points will cluster together. You can
    flip branches of the cluster dendrogram by CTRL+click on the branching points. Can you
    identify some interesting cluster.

This covers the basic analysis. The kinase functionality scores can be analyzed in a similar way
to protein quantification, i.e. via clustering, statistical testing etc.

