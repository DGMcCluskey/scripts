import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import os
os.chdir("C:\\Users\\dan94\\OneDrive - University College London\\.UCL Senior Research Fellow\\scripts")
os.getcwd()

#%%
counts = csr_matrix(
    np.random.default_rng().poisson(1, size=(100, 2000)), dtype=np.float32
)
adata = ad.AnnData(counts)
adata

adata.X

#%%
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
print(adata.obs_names[:10])

#%%
ct = np.random.default_rng().choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
adata.obs
#%%
bdata = adata[adata.obs.cell_type == "B"]
bdata
#%%
adata.obsm["X_umap"] = np.random.default_rng().normal(0, 1, size=(adata.n_obs, 2))
adata.varm["gene_stuff"] = np.random.default_rng().normal(
    0, 1, size=(adata.n_vars, 5)
)
adata.obsm
adata
#%%
adata.uns["random"] = [1, 2, 3]
adata.uns
#%%
adata.layers["log_transformed"] = np.log1p(adata.X)
adata
#%%
(adata.X != adata.layers["log_transformed"]).nnz == 0
#%%
test = adata.to_df(layer="log_transformed")
#%%
os.chdir("C:\\Users\\dan94\\OneDrive - University College London\\.UCL Senior Research Fellow\\Python tutorial scRNA seq\\")
#%%
adata.write("my_results.h5ad", compression="gzip")
#%%
adata_new = ad.read_h5ad("my_results.h5ad")
adata_new
#%%
obs_meta = pd.DataFrame(
    {
        "time_yr": np.random.default_rng().choice([0, 2, 4, 8], adata.n_obs),
        "subject_id": np.random.default_rng().choice(
            ["subject 1", "subject 2", "subject 4", "subject 8"], adata.n_obs
        ),
        "instrument_type": np.random.default_rng().choice(
            ["type a", "type b"], adata.n_obs
        ),
        "site": np.random.default_rng().choice(["site x", "site y"], adata.n_obs),
    },
    index=adata.obs.index,  # these are the same IDs of observations as above!
)
#%%
adata = ad.AnnData(adata.X, obs=obs_meta, var=adata.var)
adata
#%%
adata_view = adata[:5, ["Gene_1", "Gene_3"]]
adata_view
#%%
adata
#%%
adata_subset = adata[:5, ["Gene_1", "Gene_3"]].copy()

