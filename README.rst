

.. raw:: html

    <img src="logo_RASpy.png" height="300x">

=========================================================
Raspy : Single-Cell Reaction Activity Scores in Python
=========================================================


.. contents:: Table of Contents
   :depth: 2

***************
Install
***************

* Install Python (from release version 3.0 on)

* You can find the requirements in the `requirements.txt` file and install them through the following commands:
::
  
    pip install -r requirements.txt

* Using terminal, navigate to the favourite installation directory and run the following Git command:
::

    git clone https://github.com/compBtBs/RASpy.git


***************
Main Features
***************

* RAS computation
* RAS cluster analysis
* Differential Reaction Expression analysis
* Color a metabolic map (SVG Escher format) using RAS fold-change

***********************
Examples and Tutorials
***********************

To easily understand how to use RASPY, we invite you to try our notebook tutorial `Tutorial <https://github.com/CompBtBs/RASpy/blob/main/tutorial_RASpy.ipynb>`_

Other examples of use can be found in the `notebook_examples <https://github.com/CompBtBs/RASpy/tree/main/notebook_examples>`_ directory.

.. list-table:: Notebooks
   :widths: 25 15
   :header-rows: 1

   * - Name
     - Link
   * - How to prepare a count matrix for RAS computation
     - `Notebook1 <https://github.com/CompBtBs/RASpy/blob/main/notebook_examples/Pre-processing%20of%20the%20count%20matrix.ipynb>`_
   * - How to compute a RAS matrix
     - `Notebook2 <https://github.com/CompBtBs/RASpy/blob/main/notebook_examples/Ras%20computation.ipynb>`_
   * - How to perform RAS cluster analysis
     - `Notebook3 <https://github.com/CompBtBs/RASpy/blob/main/notebook_examples/Ras%20cluster%20analysis.ipynb>`_
   * - How to remove cell cycle effect from RAS matrix
     - `Notebook4 <https://github.com/CompBtBs/RASpy/blob/main/notebook_examples/Cell%20cycle%20removal%20on%20RAS%20matrix.ipynb>`_

**********************************************
Available datasets and metabolic models
**********************************************

Load a dataset
============================

You can load any count matrix file using any read function from scanpy tool. For example, you can load a h5ad file in the following way:

.. code-block:: python

    import scanpy as sc
    adata=sc.read_h5ad(name_of_count_matrix)

The dataset is saved in an anndata.AnnData object. We provided two example of public datasets in the datasets directory. 
* The  `first dataset  <https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-86618/downloads>`_ (E-GEOD-86618), is provided as TPM matrix and considers human lung epithelial cell types involved in the pathogenesis of Idiopathic pulmonary fibrosis (IPF). 
* The  `second dataset <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110949>`_ (GSE110949), is provided as Raw count matrix, and considers MDA-MB-231 cell line adapted to culture in media containing 0 mM or 2 mM metformin. 

Load a metabolic model
============================

You can load any metabolic model using any read function from cobrapy tool. For example, you can load a SBML model in the following way:

.. code-block:: python

   from cobra.io import read_sbml_model
   model=read_sbml_model(name_of_sbml_model)

As example of metabolic model, we have included the RECON3D model from http://bigg.ucsd.edu/models, converting the gene annotation both in ENSG and Gene symbol.


Compute RAS
============================

Once you have a metabolic model and a count matrix, you can use a RAS_computation object  in the following way:

.. code-block:: python

    from classRASpy import RAS_computation as rc
    import scanpy as sc
    #%% inizialize ras object
    ras_object=rc(adata,model)
    #%% Compute ras
    ras_adata=ras_object.compute()

The RAS dataset is saved in an anndata.AnnData object. See `Notebook2 <https://github.com/CompBtBs/RASpy/blob/main/notebook_examples/Ras%20computation.ipynb>`_ for more details.

**WARNING:** Make sure that gene annotation for count matrix and metabolic model must be the same.

RAS clustering
============================

Once you have computed the RAS dataset, you can perform a cluster analysis, using the Scanpy tool. For example,
you can use the following code to clusters the cells with the Leiden algorithm. 

.. code-block:: python

    import scanpy as sc
    sc.tl.pca(ras_adata, svd_solver='arpack'))
    sc.pp.neighbors(ras_adata)
    sc.tl.leiden(ras_adata)
    sc.tl.umap(ras_adata)
    sc.pl.umap(ras_adata, color=["leiden"])


See `Notebook3 <https://github.com/CompBtBs/RASpy/blob/main/notebook_examples/Ras%20cluster%20analysis.ipynb>`_ for more details.

Compute fold change between two groups
======================================================

Suppose that you want to characterize the metabolic differences between two groups of cells (e.g.  cancer vs normal cells). Starting from the ras_adata matrix
you can use the computer_diff method to obtain a list of reactions whose RASs results statistically different (up-regulated or down-regulated) as follow:

.. code-block:: python

    df=ras_object.compute_diff(ras_adata,name_feature)

where name_feature is the key of the observations grouping to consider. 

Color a metabolic map using RAS fold-change
========================================================

Once you have obtained the dataframe of statistically different reactions (up- or down- regulated) between two groups of cells, you can visualize it
on a metabolic map (in ESCHER svg format) using the colorMap method

.. code-block:: python

    from ras import RAS_map
    import numpy as np
    mappa=RAS_map()
    image=mappa.colorMap(mapNetwork,mapNetwork2,df_comparison)

where mapNetwork is the name of the SVG input metabolic map and mapNetwork2 is the name of the SVG output metabolic map. Up regulated reaction are coloured in red, whereas down-regulated reaction are coloured in blue.

*****************************
Team
*****************************

- Bruno Galuzzi <bruno.galuzzi@unimib.it> Implementation and conceptualization
- Davide Maspero <davide.maspero@unimib.it> Conceptualization
- Chiara Damiani <chiara.damiani@unimib.it> Conceptualization, Supervision
