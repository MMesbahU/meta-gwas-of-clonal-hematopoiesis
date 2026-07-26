import scanpy as sc
import os

manifest_path = "/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files/file_list.txt"

with open(manifest_path) as f:
    for i, line in enumerate(f, 1):
        parts = line.strip().split('\t')
        h5ad_file = parts[0]
        sample_name = os.path.basename(h5ad_file).replace('.h5ad', '')
        
        print(f"\n=== Line {i}: {sample_name} ===")
        
        if not os.path.exists(h5ad_file):
            print("  ERROR: file not found")
            continue
        
        # Read only obs (metadata) — much faster than loading full matrix
        adata = sc.read_h5ad(h5ad_file, backed='r')
        
        print(f"  N cells: {adata.n_obs:,}")
        print(f"  obs columns: {adata.obs.columns.tolist()}")
        
        # Print value counts for likely cell type columns
        for col in adata.obs.columns:
            if any(x in col.lower() for x in ['type', 'anno', 'cluster',
                                                'label', 'class', 'leiden']):
                n_unique = adata.obs[col].nunique()
                print(f"  [{col}] — {n_unique} unique values:")
                print(f"    {adata.obs[col].value_counts().head(5).to_dict()}")
        
        adata.file.close()  # important when using backed='r'
        
        
        
        
        
# === Line 1: BL_hashing ===
#   N cells: 9,278
#   obs columns: ['assignment', 'n_genes', 'n_counts', 'percent_mito', 'scale', 'leiden_labels', 'doublet_score', 'pred_dbl', 'anno']
#   [leiden_labels] — 12 unique values:
#     {'1': 1574, '2': 1486, '3': 1223, '4': 882, '5': 806}
#   [anno] — 10 unique values:
#     {'CD4+ naive T cells': 1574, 'CD14+ monocytes': 1486, 'Cytotoxic T cells': 1370, 'T helper cells': 1223, 'CD8+ naive T cells': 882}
# 
# === Line 2: BL_standard_design ===
#   N cells: 323,269
#   obs columns: ['n_genes', 'Channel', 'n_counts', 'percent_mito', 'scale', 'Group', 'leiden_labels', 'Donor', 'doublet_score', 'pred_dbl', 'anno']
#   [leiden_labels] — 20 unique values:
#     {'1': 41152, '2': 38627, '3': 36548, '4': 30154, '5': 29362}
#   [anno] — 15 unique values:
#     {'CD4+ naive T cells': 60215, 'Cytotoxic T cells': 48658, 'CD14+ monocytes': 38627, 'T helper cells': 36548, 'NK cells': 31498}
# 
# === Line 3: BM_pooling_and_control_filtered ===
#   N cells: 128,349
#   obs columns: ['n_genes', 'n_counts', 'percent_mito', 'Channel', 'assignment', 'anno']
#   [anno] — 20 unique values:
#     {'CD14+ monocytes': 26419, 'T helper cells': 13251, 'Cytotoxic T cells': 12869, 'Erythroid cells': 12660, 'CD4+ naive T cells': 10878}
# 
# === Line 4: BM_standard_design ===
#   N cells: 266,271
#   obs columns: ['n_genes', 'Channel', 'n_counts', 'percent_mito', 'scale', 'Group', 'leiden_labels', 'Donor', 'doublet_score', 'pred_dbl', 'leiden_labels_split', 'anno']
#   [leiden_labels] — 24 unique values:
#     {'1': 44009, '2': 24819, '3': 24797, '4': 23515, '5': 20877}
#   [leiden_labels_split] — 25 unique values:
#     {'1': 44009, '2': 24819, '3': 24797, '4': 23515, '5': 20877}
#   [anno] — 21 unique values:
#     {'CD4+ naive T cells': 44009, 'Cytotoxic T cells': 39450, 'CD14+ monocytes': 37570, 'T helper cells': 24797, 'Naive B cells': 23515}
# 
# === Line 5: CB_extra ===
#   N cells: 133,526
#   obs columns: ['n_genes', 'n_counts', 'Channel', 'Donor', 'percent_mito', 'scale', 'leiden_labels', 'doublet_score', 'pred_dbl', 'anno']
#   [leiden_labels] — 21 unique values:
#     {'1': 16041, '2': 13839, '3': 12996, '4': 12632, '5': 12238}
#   [anno] — 13 unique values:
#     {'T cells': 66718, 'CD8+ naive T cells': 17686, 'CD14+ monocytes': 16041, 'Naive B cells': 12632, 'CD4+ naive T cells': 8694}
# 
# === Line 6: CB_pooling_and_hashing ===
#   N cells: 20,989
#   obs columns: ['n_genes', 'n_counts', 'percent_mito', 'Channel', 'assignment', 'anno']
#   [anno] — 7 unique values:
#     {'T cells': 13992, 'Monocytes & DCs': 3008, 'B cells': 2394, 'NK cells': 981, 'HSCs': 248}
# 
# === Line 7: CB_standard_design ===
#   N cells: 239,544
#   obs columns: ['n_genes', 'Channel', 'n_counts', 'percent_mito', 'scale', 'Group', 'leiden_labels', 'Donor', 'doublet_score', 'pred_dbl', 'id2', 'leiden_labels_split', 'anno']
#   [leiden_labels] — 19 unique values:
#     {'1': 32701, '2': 32071, '3': 27183, '4': 25270, '5': 23964}
#   [leiden_labels_split] — 20 unique values:
#     {'1': 32701, '2': 32071, '3': 27183, '4': 25270, '5': 23964}
#   [anno] — 13 unique values:
#     {'CD4+ naive T cells': 116290, 'CD8+ naive T cells': 33643, 'Naive B cells': 32701, 'CD14+ monocytes': 23964, 'NK cells': 12066}
# 
# === Line 8: onek1k_sc_filtered ===
#   N cells: 1,174,789
#   obs columns: ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'donor_id', 'pool_number', 'predicted.celltype.l2', 'predicted.celltype.l2.score', 'age', 'tissue_ontology_term_id', 'assay_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'is_primary_data', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid', 'n_genes']
#   [predicted.celltype.l2] — 31 unique values:
#     {'CD4 TCM': 267649, 'CD4 Naive': 240144, 'NK': 160400, 'CD8 TEM': 153495, 'B naive': 60353}
#   [predicted.celltype.l2.score] — 947168 unique values:
#     {1.0: 67250, 1.0000000000000002: 51654, 0.9999999999999999: 38016, 0.9999999999999998: 27012, 1.0000000000000004: 14701}
#   [cell_type_ontology_term_id] — 29 unique values:
#     {'CL:0000904': 267649, 'CL:0000895': 240144, 'CL:0000623': 162130, 'CL:0000913': 153495, 'CL:0000788': 60353}
#   [suspension_type] — 1 unique values:
#     {'cell': 1174789}
#   [tissue_type] — 1 unique values:
#     {'tissue': 1174789}
#   [cell_type] — 29 unique values:
#     {'central memory CD4-positive, alpha-beta T cell': 267649, 'naive thymus-derived CD4-positive, alpha-beta T cell': 240144, 'natural killer cell': 162130, 'effector memory CD8-positive, alpha-beta T cell': 153495, 'naive B cell': 60353}