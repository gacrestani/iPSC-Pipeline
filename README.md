# iPSC Pipeline

iPSC Pipeline is a pipeline developed with Snakemake to apply several single cell analysis tools into query single cell transcriptomic datasets. The pipeline requires one  reference dataset, with known cell types in its metadata and one query dataset with a metadata file. The output of the pipeline is a series of figures that can be used to generate reports and/or scientific papers.

This pipeline was built for the final undergrad paper (TCC) of Giovanni Alberto Crestani. The paper will be available soon.

## Installation

The only installation needed for this pipeline is Snakemake. To install Snakemake, go to [Snakemake installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Quick start

1. Clone or download this repo and go to the directory.

``` sh
 git clone https://github.com/gacrestani/ipsc-pipeline
 cd ipsc-pipeline
 ```

2. Copy the `singularity.img` Singularity container into [`workflow/envs/`](workflow/envs). You can download the image [here](https://www.dropbox.com/s/680doaawgn7avdb/singularity.img?dl=0).

3. Run the pipeline with
``` sh
snakemake -c4 --use-singularity
```
You can change `4` to be the number of cores used by Snakemake.

The pipeline will run and generate the figures using the [La Manno et al (2016)](https://pubmed.ncbi.nlm.nih.gov/27716510/) dataset.


## Customization

Edit the configuration file [`config/config.yml`](config/config.yml).

  - `datasets`: the directory names of your datasets. There must be a different directory for each query/reference dataset pair. The directory must contain at the following four files:

    - `query_exp.rda` is the `rda` file of the expression matrix of your query dataset. This should be a [dgCMatrix](https://cran.r-project.org/web/packages/Matrix/Matrix.pdf) where columns are samples and rows are genes.

    - `query_meta.rda` is the `rda` file of the metadata of your query dataset. The `rownames` must be samples.

    - `train_exp.rda` is the `rda` file of the expression matrix of your reference dataset. Its structure should be similar to `query_exp.rda`.

    - `train_meta.rda` is the `rda` file of the metadata of your reference dataset. The `rownames` must be samples.

  - `Cell_ID_colname`: column name of the sample IDs for the reference and query dataset.
  - `Cell_type_colname`: column name of the cell types for the reference dataset.
  - `Cell_description_colname`: column name of the cell grouping in the query dataset. This should be some metadata that classifies the samples into diferent groups. It is optional.
  - `fig_width`: width of generated figures, in pixels.
  - `fig_height`: height of generated figures, in pixels.
  - `fig_res`: resolution of generated figures, in ppi.

Additionaly, there are some tool-specifi configurations. For example, you can set a cell type to be focused on the singleCellNet query violin plot.

You can leave these options as-is if you'd like to first make sure the workflow runs without error on your machine before using your own dataset and custom parameters.

## Additional configurations

The singularity containter `singularity.img` was generated from a Docker container named `gacrestani/ipsc_pipeline`. The code to generate the Docker container is in [`workflow/envs/docker`](workflow/envs/docker). The code used to generate the Singularity container from the Docker container was:
```sh
sudo singularity build singularity.img docker-daemon://gacrestani/ipsc_pipeline:latest
```
