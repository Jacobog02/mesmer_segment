# MESMER Segmentation

Purpose of this package is to apply mesmer segmentation to SMI CellComposite images and use those segmentations to recompute downstream outputs like expression matrices.


# INSTALL

This package hasnt been uploaded to pip current (as of 5/5/22) 

To use this function on commandline you want to clone this repository followed by 

```
git clone https://github.com/Jacobog02/mesmer_segment.git # get repository
cd mesmer_segment # go there 
pip install . # install local package to your python
```

# Development

If you are interested in modifying this package and rerunning pipelines live you must start a new python virtualenv, then call the setup script with the `develop` parameter. You can then pip install the package and when you run the command it recompiles on the fly


```
python3 -m venv /path/to/mesmer_env ## make the env
source /path/to/mesmer_env/bin/activate ## Spark up the env

cd /path/to/mesmer ## Go to the github repo locally

python3 setup.py develop ## Config package to recompile on the fly

pip install . ## Install the package 

mesmer_segment --help ## Yields help screen from onfly compliation

```

In order to both segment and generate an expression matrix, polygon vertices, etc. use the following command:
```
mesmer_segment -i /path/to/outs/ --save-npz --build-outs
```

The outputs can then be loaded into Seurat (feat/imaging branch) with:
```
data <- ReadNanostring(
  mtx.file = "/path/to/mesmer_exprMat_file.csv",
  metadata.file = "/path/to/mesmer_metadata_file.csv",
  molecules.file = "/path/to/mesmer_tx_file.csv",
  segmentations.file = "/path/to/mesmer_polygons_file.csv",
  type = c("centroids", "segmentation")
)
segs <- CreateSegmentation(data$segmentations)
cents <- CreateCentroids(data$centroids)
segmentations.data <- list(
  "centroids" = cents,
  "segmentation" = segs
)
coords <- CreateFOV(
  coords = segmentations.data,
  type = c("centroids", "segmentation"),
  molecules = data$pixels,
  assay = 'RNA'
)
obj <- CreateSeuratObject(counts = data$matrix, assay = 'RNA')

# select cells in segmentations and in expression matrix
cells <- intersect(Cells(obj), Cells(x = coords, boundary = "centroids"))
coords <- subset(x = coords, cells = cells)

# Add coords as an Image with name `SMI`
obj[['SMI']] <- coords
```
