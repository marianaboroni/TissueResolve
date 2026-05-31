import anndata as ad
import numpy as np

from tissueresolve.protocol import detect as pd


def test_detect_protocol_from_anndata():
    adata = ad.AnnData(X=np.ones((5, 3)))
    adata.uns["platform"] = "10x Genomics"
    det = pd.detect_protocol_from_anndata(adata)
    assert "protocol" in det
    assert det["protocol"] != ""
