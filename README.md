# Pan-cancer analysis

The code contained in this folder can be used to reproduce the pancancer analysis presented in Cabassi and Kirk (2020). The data used for the analysis are those of Hoadley et al. (2014), which can be downloaded from [Synapse](https://www.synapse.org/#!Synapse:syn2468297/wiki/64259). The links to each dataset can be found at the beginning of the R scripts where needed.

## Folder structure

0) **Preprocessing**. Each layer of the  dataset is preprocessed in the same way as in Hoadley et al. (2014). A clustering of the patient is found in each layer using the same clustering algorithm as Hoadley et al. (2014).
1) **COCA**. Integrative clustering is performed using Cluster-Of-Clusters Analysis (Cabassi and Kirk, 2020) and the clusters found in the previous step.
2) **Kernel**. For each layer, one kernel is built using Consensus Clustering (Monti et al. 2003) and the most appropriate clustering algorithm.
3) **KLIC**. Integrative clustering is performed using Kernel Learning Integrative Clustering (Cabassi and Kirk, 2020) and the kernels built in the previous step.

## References

- Cabassi, A., Kirk, P. D. W., 2020. Multiple kernel learning for integrative consensus clustering of genomic datasets. Bioinformatics, btaa593.

- Hoadley, K. A., Yau, C., Wolf, D. M., Cherniack, A. D., Tamborero, D., Ng, S., Leiserson, M. D., Niu, B., McLellan, M. D., Uzunangelov, V. and Zhang, J., 2014. Multiplatform analysis of 12 cancer types reveals molecular classification within and across tissues of origin. Cell, 158(4), pp.929-944.

- Monti, S., Tamayo, P., Mesirov, J. and Golub, T., 2003. Consensus clustering: a resampling-based method for class discovery and visualization of gene expression microarray data. Machine learning, 52(1-2), pp.91-118.

## Contact

Please feel free to contact me (Alessandra) should you need any help using this code. You can find my up-to-date email address in my GitHub profile.
