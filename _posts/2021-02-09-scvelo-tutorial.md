---
layout: post
title: scVelo
subheading: using scVelo to analyze transcriptional dynamics in scRNA-seq data
author: Samuel Morabito
categories: pseudotime
banner: /assets/images/scVelo/scvelo_banner.png
tags: single-cell pseudotime velocity python scanpy
sidebar: []
usemathjax: true
---



## Introduction

In this tutorial, I will cover how to use the Python package [scVelo](https://scvelo.readthedocs.io/index.html) to perform RNA velocity analysis in single-cell RNA-seq data (scRNA-seq). scVelo was published in 2020 in [Nature Biotechnology](https://www.nature.com/articles/s41587-020-0591-3), making several improvements from the original [RNA velocity study](https://www.nature.com/articles/s41586-018-0414-6) and its accomanpying software [velocyto](http://velocyto.org/).

Briefly, RNA velocity analysis allows us to infer transcriptional dynamics that are not directly observed in a scRNA-seq experiment using a mathematical model of transcriptional kinetics. We can use RNA velocity to determine if a gene of interest is being induced or repressed in a give cell population of interest. Moreover, we can extrapolate this information to predict cell fate decision via pseudotime trajectories.

The majority of this tutorial is taken from the [scVelo](https://scvelo.readthedocs.io/index.html) documentation.

## Step -1: Convert data from Seurat to Python / anndata

For this tutorial, I am starting with a mouse brain dataset that contains cells from disease and control samples. I have already performed the primary data processing (filtering, normalization, clustering, batch alignment, dimensionality reduction) using Seurat in R. First I will go over the code that I used to convert my Seurat object into a format that is usable in scVelo (anndata format). It is possible to use the [SeuratDisk](https://mojaveazure.github.io/seurat-disk/reference/SeuratDisk-package.html) R package to convert between Seurat and anndata formats, but this software is in early development stages and I have had mixed results when using it.



```r

# assuming that you have some Seurat object called seurat_obj:

# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(out_data_dir, 'counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)
```


```python
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")

# load gene names:
with open("gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['SampleID'], frameon=False, save=True)

# save dataset as anndata format
adata.write('my_data.h5ad')

# reload dataset
adata = sc.read_h5ad('my_data.h5ad')
```

## Step 0: Constructing spliced and unspliced counts matrices

Rather than using the same UMI-based genes-by-counts matrix that we used in Seurat, we need to have a matrix for spliced and unspliced transcripts. We can construct this matrix using the [velocyto command line tool](https://velocyto.org/velocyto.py/tutorial/cli.html),  or using [Kallisto-Bustools](https://bustools.github.io/BUS_notebooks_R/velocity.html). Here I am using the velocyto command line tool, simply because I had a working script from before Kallisto supported RNA velocity, so I have personally never bothered to try Kallisto.

The velocyto command line tool has a function that works directly from the cellranger output directory, but it also can be used on any single-cell platform as long as you provide a .bam file. We also have to supply a reference .gtf annotation for your species (mm10 used here), and optionally you can provide a .gtf to mask repeat regions (recommended by velocyto).


```bash
repeats="/path/to/repeats/mm10_rmsk.gtf"
transcriptome="/path/to/annoation/file/gencode.vM25.annotation.gtf"
cellranger_output="/path/to/cellranger/output/"

velocyto run10x -m $repeats \
                $cellranger_output \
                $transcriptome
```

## Step 1: Load data

Now that we have our input data properly formatted, we can load it into python. Velocyto created a separate spliced and unspliced matrix for each sample, so we first have to merge the different samples into one object. Additionally, I am reformatting the cell barcodes to match my anndata object with the full genes-by-cells data.


```python
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
```


```python
adata = sc.read_h5ad('my_data.h5ad')
```


```python
# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('Sample1.loom', cache=True)
ldata2 = scv.read('Sample2.loom', cache=True)
ldata3 = scv.read('Sample3.loom', cache=True)
```

    Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    Variable names are not unique. To make them unique, call `.var_names_make_unique`.



```python
# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_9' for bc in barcodes]
ldata3.obs.index = barcodes

```


```python
# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
```


```python
# concatenate the three loom
ldata = ldata1.concatenate([ldata2, ldata3])
```


```python
# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)
```


```python
# plot umap to check
sc.pl.umap(adata, color='celltype', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')
```



![png](/assets/images/scVelo/output_14_0.png)



## Part 2: Computing RNA velocity using scVelo

Finally we can actually use scVelo to compute RNA velocity. scVelo allows the user to use the steady-state model from the original 2018 publication as well as their updated dynamical model from the 2020 publication.

First we inspect the proportion of spliced and unspliced transcripts in each of our cell clusters. Next we perform some pre-processing, and then we compute RNA velocity using the steady-state model (stochastic option).


```python
scv.pl.proportions(adata, groupby='celltype_full')
```



![png](/assets/images/scVelo/output_16_0.png)




```python
# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
```

    WARNING: Did not normalize X as it looks processed already. To enforce normalization, set `enforce=True`.
    Normalized count data: spliced, unspliced.
    WARNING: Did not modify X as it looks preprocessed already.
    computing neighbors
        finished (0:00:14) --> added
        'distances' and 'connectivities', weighted adjacency matrices (adata.obsp)
    computing moments based on connectivities
        finished (0:00:17) --> added
        'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)



```python
# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
```

    computing velocities
        finished (0:01:15) --> added
        'velocity', velocity vectors for each individual cell (adata.layers)
    computing velocity graph
        finished (0:15:00) --> added
        'velocity_graph', sparse matrix with cosine correlations (adata.uns)


## Part 2.1: Visualize velocity fields

We can get a broad visualiztion of RNA velocity across all genes and all cells by visualizing a vector field on top of the 2D dimensional reduction.


```python
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
```

    computing velocity embedding
        finished (0:00:06) --> added
        'velocity_umap', embedded velocity vectors (adata.obsm)
    saving figure to file ./figures/scvelo_embedding.pdf




![png](/assets/images/scVelo/output_20_1.png)




```python
scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltype', save='embedding_grid.pdf', title='', scale=0.25)
```

    saving figure to file ./figures/scvelo_embedding_grid.pdf




![png](/assets/images/scVelo/output_21_1.png)




```python
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['celltype', 'condition'], save='embedding_stream.pdf', title='')
```

    figure cannot be saved as pdf, using png instead.
    saving figure to file ./figures/scvelo_embedding_stream.pdf.png




![png](/assets/images/scVelo/output_22_1.png)



We can also visualize the dynamics of our favorite genes. Here we show the ratio of unspliced to spliced transcripts for Ptgds, as well as the velocity and expression values as feature plots.


```python
# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Ptgds'], color='celltype')
```



![png](/assets/images/scVelo/output_24_0.png)



## Part 3: Downstream analysis

Here we show how to identify highly dynamic genes, compute a measure of coherence among neighboring cells in terms of  velocity, and perform pseudotime inference. Using the pseudotime trajectory, we can identify predicted ancestors of individual cells, and we can orient the directionality of partition-based graph abstractions (PAGA).


```python
scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
```

    ranking velocity genes
        finished (0:00:23) --> added
        'rank_velocity_genes', sorted scores by group ids (adata.uns)
        'spearmans_score', spearmans correlation scores (adata.var)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AF1</th>
      <th>AF2</th>
      <th>AF3</th>
      <th>AF4</th>
      <th>AF5</th>
      <th>AF6</th>
      <th>CP</th>
      <th>EC1</th>
      <th>EC2</th>
      <th>MM1</th>
      <th>MM2</th>
      <th>MM3</th>
      <th>OA</th>
      <th>PE</th>
      <th>PF1</th>
      <th>PF2</th>
      <th>SMC</th>
      <th>UF</th>
      <th>UI</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Adam12</td>
      <td>Zfpm2</td>
      <td>Zfpm2</td>
      <td>Ptprd</td>
      <td>Pcdh7</td>
      <td>Nav1</td>
      <td>Sgip1</td>
      <td>Tmtc2</td>
      <td>St6galnac3</td>
      <td>Elmo1</td>
      <td>Rnf150</td>
      <td>Tmcc3</td>
      <td>Eya4</td>
      <td>Cobll1</td>
      <td>Nr3c2</td>
      <td>Slc1a3</td>
      <td>Dmd</td>
      <td>Cadm2</td>
      <td>Adam10</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Ghr</td>
      <td>Efemp2</td>
      <td>Mmp14</td>
      <td>Mtch1</td>
      <td>Rtl4</td>
      <td>Zfpm2</td>
      <td>Dlc1</td>
      <td>Ptprg</td>
      <td>Plcb1</td>
      <td>Mir142hg</td>
      <td>Ophn1</td>
      <td>Gm5086</td>
      <td>Soga3</td>
      <td>Atp13a5</td>
      <td>Fbxl7</td>
      <td>Frmpd4</td>
      <td>Ctnna3</td>
      <td>Dync1i1</td>
      <td>Trim37</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Dpyd</td>
      <td>Tmem208</td>
      <td>Efna5</td>
      <td>Cped1</td>
      <td>Adamtsl1</td>
      <td>Arl6ip4</td>
      <td>Grm7</td>
      <td>Mecom</td>
      <td>Cyyr1</td>
      <td>Rbm47</td>
      <td>Frmd4b</td>
      <td>Srgap2</td>
      <td>Lrrc4c</td>
      <td>Lnx2</td>
      <td>Slc47a2</td>
      <td>4930594M22Rik</td>
      <td>Cacnb2</td>
      <td>Nrp2</td>
      <td>Esam</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Tcf4</td>
      <td>Rpl29</td>
      <td>Lrrc7</td>
      <td>Ccbe1</td>
      <td>Slc24a3</td>
      <td>Rpl29</td>
      <td>Rgs6</td>
      <td>Dnm3</td>
      <td>Ptprm</td>
      <td>Gramd1b</td>
      <td>Plekhg5</td>
      <td>Slc7a8</td>
      <td>Frmd4a</td>
      <td>Apbb2</td>
      <td>Tbx15</td>
      <td>Acss3</td>
      <td>Slc38a11</td>
      <td>Ust</td>
      <td>Sgms1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Rbms3</td>
      <td>Mtch1</td>
      <td>Dclk1</td>
      <td>Cers6</td>
      <td>Slit3</td>
      <td>Efna5</td>
      <td>Trpc3</td>
      <td>Tjp1</td>
      <td>Lrch1</td>
      <td>Lyn</td>
      <td>Runx1</td>
      <td>Qk</td>
      <td>Ptprt</td>
      <td>Egflam</td>
      <td>Mmp16</td>
      <td>Kcnk2</td>
      <td>Sox6</td>
      <td>Sorbs1</td>
      <td>Emcn</td>
    </tr>
  </tbody>
</table>
</div>




```python
kwargs = dict(frameon=False, size=10, linewidth=1.5,
              add_outline='AF6, AF1')

scv.pl.scatter(adata, df['AF6'][:3], ylabel='AF6', frameon=False, color='celltype', size=10, linewidth=1.5)
scv.pl.scatter(adata, df['AF1'][:3], ylabel='AF1', frameon=False, color='celltype', size=10, linewidth=1.5)

```



![png](/assets/images/scVelo/output_27_0.png)





![png](/assets/images/scVelo/output_27_1.png)



* Speed: length of the velocity vector
* Coherence: how well a velocity vector correlates to its neighbors


```python
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])
```

    --> added 'velocity_length' (adata.obs)
    --> added 'velocity_confidence' (adata.obs)




![png](/assets/images/scVelo/output_29_1.png)




```python
scv.pl.velocity_graph(adata, threshold=.1, color='celltype')
```



![png](/assets/images/scVelo/output_30_0.png)




```python
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)
```



![png](/assets/images/scVelo/output_31_0.png)




```python
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')
```



![png](/assets/images/scVelo/output_32_0.png)




```python
# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
```


```python
scv.tl.paga(adata, groups='celltype')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
#df.style.background_gradient(cmap='Blues').format('{:.2g}')
```

    running PAGA using priors: ['velocity_pseudotime']
        finished (0:00:04) --> added
        'paga/connectivities', connectivities adjacency (adata.uns)
        'paga/connectivities_tree', connectivities subtree (adata.uns)
        'paga/transitions_confidence', velocity transitions (adata.uns)



```python
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)
```

    WARNING: Invalid color key. Using grey instead.




![png](/assets/images/scVelo/output_35_1.png)



## Part 4: Analyzing a specific cell population

In a scRNA-seq dataset comprised of multiple cell lineages, it may be more relevant to perfom RNA velocity analysis separately for each major cell population Here we are going to perform RNA velocity analysis using only our fibroblast clusters. Additionally, we are going to use the dynamical model (from 2020 paper) whereas earlier we used the steady-state model (2018 paper).


```python
cur_celltypes = ['AF1', 'AF2', 'AF3', 'AF4', 'AF5', 'AF6', 'PF1', 'PF2']
adata_subset = adata[adata.obs['celltype'].isin(cur_celltypes)]
```


```python
sc.pl.umap(adata_subset, color=['celltype', 'condition'], frameon=False, title=['', ''])
```



![png](/assets/images/scVelo/output_38_0.png)




```python
sc.pp.neighbors(adata_subset, n_neighbors=15, use_rep='X_pca')

# pre-process
scv.pp.filter_and_normalize(adata_subset)
scv.pp.moments(adata_subset)
```

    WARNING: Did not normalize X as it looks processed already. To enforce normalization, set `enforce=True`.
    WARNING: Did not normalize spliced as it looks processed already. To enforce normalization, set `enforce=True`.
    WARNING: Did not normalize unspliced as it looks processed already. To enforce normalization, set `enforce=True`.
    WARNING: Did not modify X as it looks preprocessed already.
    computing neighbors
        finished (0:00:03) --> added
        'distances' and 'connectivities', weighted adjacency matrices (adata.obsp)
    computing moments based on connectivities
        finished (0:00:11) --> added
        'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)
    computing velocities
        finished (0:00:31) --> added
        'velocity', velocity vectors for each individual cell (adata.layers)
    computing velocity graph
        finished (0:04:00) --> added
        'velocity_graph', sparse matrix with cosine correlations (adata.uns)


In this step, transcriptional dynamics are computed. This can take a lot longer than other steps, and may require the use of a high performance compute cluster such as HPC3 at UCI. This step took about 90 minutes to run on my laptop for ~15k cells.


```python
scv.tl.recover_dynamics(adata_subset)
```

    recovering dynamics
        finished (0:55:01) --> added
        'fit_pars', fitted parameters for splicing dynamics (adata.var)



```python
scv.tl.velocity(adata_subset, mode='dynamical')
scv.tl.velocity_graph(adata_subset)
```

    computing velocities
        finished (0:01:54) --> added
        'velocity', velocity vectors for each individual cell (adata.layers)
    computing velocity graph
        finished (0:01:37) --> added
        'velocity_graph', sparse matrix with cosine correlations (adata.uns)



```python
scv.pl.velocity_embedding_stream(adata_subset, basis='umap', color=['celltype', 'condition'], save='embedding_stream.pdf', title='')
```

    computing velocity embedding
        finished (0:00:08) --> added
        'velocity_umap', embedded velocity vectors (adata.obsm)
    figure cannot be saved as pdf, using png instead.
    saving figure to file ./figures/scvelo_embedding_stream.pdf.png




![png](/assets/images/scVelo/output_43_1.png)



Using the dynamical model, we can actually look at the transcritpion rate, splicing rate, and degradation rate.


```python
df = adata_subset.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata_subset, 'fit*', dropna=True).head()
```



![png](/assets/images/scVelo/output_45_0.png)






<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fit_alpha</th>
      <th>fit_beta</th>
      <th>fit_gamma</th>
      <th>fit_t_</th>
      <th>fit_scaling</th>
      <th>fit_std_u</th>
      <th>fit_std_s</th>
      <th>fit_likelihood</th>
      <th>fit_u0</th>
      <th>fit_s0</th>
      <th>fit_pval_steady</th>
      <th>fit_steady_u</th>
      <th>fit_steady_s</th>
      <th>fit_variance</th>
      <th>fit_alignment_scaling</th>
      <th>fit_r2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Sox17</th>
      <td>0.898695</td>
      <td>6.110614</td>
      <td>0.484805</td>
      <td>10.475681</td>
      <td>0.105126</td>
      <td>0.037980</td>
      <td>0.421499</td>
      <td>0.000018</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.492005</td>
      <td>0.104282</td>
      <td>1.640273</td>
      <td>1.071412</td>
      <td>1.816156</td>
      <td>0.661397</td>
    </tr>
    <tr>
      <th>Pcmtd1</th>
      <td>0.211919</td>
      <td>0.720424</td>
      <td>0.264111</td>
      <td>12.202008</td>
      <td>0.248776</td>
      <td>0.086398</td>
      <td>0.164096</td>
      <td>0.232746</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.336931</td>
      <td>0.270257</td>
      <td>0.446798</td>
      <td>0.749243</td>
      <td>2.372727</td>
      <td>0.108137</td>
    </tr>
    <tr>
      <th>Sgk3</th>
      <td>0.031516</td>
      <td>0.226350</td>
      <td>0.197472</td>
      <td>10.975034</td>
      <td>1.652986</td>
      <td>0.050927</td>
      <td>0.038999</td>
      <td>0.177014</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.159321</td>
      <td>0.164692</td>
      <td>0.126664</td>
      <td>1.650863</td>
      <td>3.383189</td>
      <td>0.087770</td>
    </tr>
    <tr>
      <th>Prex2</th>
      <td>0.209454</td>
      <td>0.229509</td>
      <td>0.231223</td>
      <td>17.965210</td>
      <td>1.644350</td>
      <td>0.334659</td>
      <td>0.288491</td>
      <td>0.303137</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.430445</td>
      <td>0.757209</td>
      <td>0.657950</td>
      <td>0.545063</td>
      <td>3.398314</td>
      <td>0.796931</td>
    </tr>
    <tr>
      <th>Sulf1</th>
      <td>0.184210</td>
      <td>0.221469</td>
      <td>0.183979</td>
      <td>14.692180</td>
      <td>1.693538</td>
      <td>0.300075</td>
      <td>0.217056</td>
      <td>0.243579</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.311456</td>
      <td>0.848147</td>
      <td>0.743920</td>
      <td>0.904572</td>
      <td>3.110569</td>
      <td>0.164060</td>
    </tr>
  </tbody>
</table>
</div>



Similar to pseudotime, 'latent time' is computed from the dynamical model.


```python
scv.tl.latent_time(adata_subset)
scv.pl.scatter(adata_subset, color='latent_time', color_map='gnuplot', size=80)
```

    computing latent time using root_cells as prior
        finished (0:00:48) --> added
        'latent_time', shared time (adata.obs)




![png](/assets/images/scVelo/output_47_1.png)




```python
top_genes = adata_subset.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata_subset, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100)

```



![png](/assets/images/scVelo/output_48_0.png)




```python
top_genes = adata_subset.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata_subset, color='celltype', basis=top_genes[:15], ncols=5, frameon=False)
```



![png](/assets/images/scVelo/output_49_0.png)




```python
var_names = ['Rps3', 'Col1a2', 'Tmeff2']
scv.pl.scatter(adata_subset, var_names, color='celltype', frameon=False)
scv.pl.scatter(adata_subset, x='latent_time', y=var_names, color='celltype', frameon=False)

```



![png](/assets/images/scVelo/output_50_0.png)





![png](/assets/images/scVelo/output_50_1.png)



# Conclusion and future directions:

In this tutorial we have learned how to perform RNA velocity analysis in scRNA-seq data using scVelo, as well as some downstream applications of RNA velocity such as identifying potential cell-state transitions and identifying genes that are being induced / repressed. RNA velocity is the starting point for a cell fate analysis program called [CellRank](https://cellrank.readthedocs.io/en/latest/auto_examples/index.html), which can provide further insights for the analysis of cell lineages.


```python

```
