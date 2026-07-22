#!/usr/bin/env python3
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="star_ambient_fdr_h5ad.") as tmp:
        work = Path(tmp)
        library = work / "library"
        write_text(
            library / "matrix.mtx",
            """%%MatrixMarket matrix coordinate integer general
%
2 3 4
1 1 7
2 2 3
1 3 2
2 3 2
""",
        )
        write_text(library / "features.tsv", "guideA\tguideA\tCRISPR Guide Capture\nguideB\tguideB\tCRISPR Guide Capture\n")
        write_text(library / "barcodes.tsv", "cell1\ncell2\ncell3\n")
        write_text(
            library / "feature_per_cell.csv",
            "\n".join(
                [
                    "barcode,num_features,top_feature_index,total_deduped_umi",
                    "cell1,1,1,7",
                    "cell2,1,2,3",
                    "cell3,2,0,4",
                    "",
                ]
            ),
        )
        write_text(
            library / "pf_library_provenance.tsv",
            "key\tvalue\nlibrary_id\tlib1\nsample\tsample1\nfeature_type\tCRISPR Guide Capture\n",
        )

        counts = work / "counts.h5ad"
        obs = pd.DataFrame(index=pd.Index(["cell1-1", "cell2-1", "cell3-1"], name=None))
        var = pd.DataFrame(index=pd.Index(["gene1", "gene2"], name=None))
        ad.AnnData(X=sp.csr_matrix(np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float32)), obs=obs, var=var).write_h5ad(counts)

        crispr_analysis = work / "crispr_analysis"
        gmm_calls = crispr_analysis / "protospacer_calls_per_cell.csv"
        calls = crispr_analysis / "ambient_fdr" / "guide_fdr_calls_per_cell.csv"
        write_text(
            gmm_calls,
            "\n".join(
                [
                    "cell_barcode,num_features,feature_call,num_umis",
                    "cell1-1,1,guideA,7",
                    "cell2-1,1,guideB,3",
                    "cell3-1,0,None,0",
                    "",
                ]
            ),
        )
        write_text(
            calls,
            "\n".join(
                [
                    "cell_barcode,num_features,feature_call,num_umis,min_called_umi,max_called_umi,min_qvalue,num_features_at_default_fdr,call_status,default_fdr,caller",
                    "cell1-1,1,guideA,7,7,7,0.0001,1,singlet,0.01,ambient-fdr",
                    "cell2-1,2,guideA|guideB,9,3,6,0.002,2,multiplet,0.01,ambient-fdr",
                    "cell3-1,0,None,0,0,0,1,0,none,0.01,ambient-fdr",
                    "",
                ]
            ),
        )

        subprocess.run(
            [
                "python3",
                str(ROOT / "scripts" / "integrate_feature_library.py"),
                "--library-dir",
                str(library),
                "--feature-output-root",
                str(work / "feature_libraries"),
                "--counts-h5ad",
                str(counts),
                "--calls-csv",
                str(gmm_calls),
                "--set-generic-aliases",
            ],
            check=True,
            cwd=ROOT,
        )

        annotated = ad.read_h5ad(counts)
        expected = {
            "guide_fdr_num_umis",
            "guide_fdr_min_called_umi",
            "guide_fdr_max_called_umi",
            "guide_fdr_min_qvalue",
            "guide_fdr_call_status",
        }
        missing = expected.difference(annotated.obs.columns)
        assert not missing, missing
        assert annotated.obs.loc["cell2-1", "guide_fdr_num_umis"] == 9
        assert annotated.obs.loc["cell2-1", "guide_fdr_min_called_umi"] == 3
        assert annotated.obs.loc["cell3-1", "guide_fdr_call_status"] == "none"

        feature_h5ad = work / "feature_libraries" / "lib1" / "raw_feature_library.h5ad"
        feature_adata = ad.read_h5ad(feature_h5ad)
        assert feature_adata.obs.loc["cell1", "top_feature_name"] == "guideA"
        assert feature_adata.obs.loc["cell2", "top_feature_name"] == "guideB"
        assert feature_adata.obs.loc["cell3", "top_feature_name"] == ""

        atac = work / "atac"
        write_text(
            atac / "matrix.mtx",
            """%%MatrixMarket matrix coordinate integer general
%
1 3 3
1 1 5
1 2 6
1 3 7
""",
        )
        write_text(atac / "features.tsv", "chr1:10-20\tchr1:10-20\tPeaks\n")
        write_text(atac / "barcodes.tsv", "cell1-1\ncell2-1\ncell3-1\n")
        output_h5mu = work / "out.h5mu"
        subprocess.run(
            [
                "python3",
                str(ROOT / "scripts" / "build_multiome_mudata.py"),
                "--rna-h5ad",
                str(counts),
                "--atac-mex-dir",
                str(atac),
                "--output-h5mu",
                str(output_h5mu),
                "--all-barcodes-are-cells",
            ],
            check=True,
            cwd=ROOT,
        )

        mdata = md.read_h5mu(output_h5mu)
        assert "guide_fdr_num_umis" in mdata.obs.columns
        assert "guide_fdr_min_called_umi" in mdata.obs.columns
        assert mdata.obs.loc["cell1-1", "guide_fdr_num_umis"] == 7
        assert mdata.obs.loc["cell2-1", "guide_fdr_min_called_umi"] == 3


if __name__ == "__main__":
    main()
