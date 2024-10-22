# TBL

This repo is for the GeoMx data from TBL + EPTB. The data is provided on gdrive as raw, normalized, and negative probe counts + low-res images. Note that the count matrices include PTB sample columns, but the scripts just ignore those. 

We do not have high-res images yet and may have to follow-up with NanoString people to get these. TBL16 and TBL42 were in an early batch (batch #2) while TBL62 and EPTB22.2 were in the last batch (batch #3).

## Processed count data:
1. geomx_integrated_q3_data_2024.csv: **genes** x **samples** (annotation | tissue | ROI) 
2. geomx_integrated_raw_data_2024.csv: **genes** x **samples** 
3. geomx_negprobes_2024.csv: **probes** x **samples** 
4. nuclei_2024.csv: **samples** x **(nuclei counts, surface area of spot, X coord, Y coord, aligned reads, aligned %)**

## There are four types of samples:
1. core - center of granuloma (macrophages)
2. mantle - periphery of granuloma (fibroblasts + lymphocytes)
3. giant - [core subtype, multinucleated macrophages](https://en.wikipedia.org/wiki/Langhans_giant_cell)
4. infiltrate - B cell cluster adjacent to granuloma, ignored for mantle vs core contrast 

Scripts are only for TBL but can be adapted for EPTB. 

### GeoMx_TBL_PCA.R
Script to generate PCA and dimensionality reduction plots with and without batch correction ([ComBat](https://bioconductor.org/packages/devel/bioc/vignettes/sva/inst/doc/sva.pdf)) 
### GeoMx_TBL_DGE.R
Script to run DGE for core vs mantle and identify enriched pathways 
### GeoMx_TBL_WGCNA.R
Script to run [WGCNA](https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0) on batch corrected counts data for 3k most variable genes
