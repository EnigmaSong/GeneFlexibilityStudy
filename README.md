This repository contains code and scripts for the study “Flexibility in Gene Coexpression at Developmental and Evolutionary Timescales” (Fischer et al., 2025).

# Input Data

The required input is a matrix where rows represent samples and columns represent genes. 

# Features

* Random Projection Tests – test equality of covariance matrices across groups. (Test/rand_proj_cov.R)
* Correlation Network Construction – based on Cai & Liu (2016). (Test/cor_CL.R)
* Network Summary Plots – for exploratory visualization of network topology (Maugis et al., 2017). (nethist::violin_netsummary)
* Network Comparison Tests – two-sample network test (Shao et al., 2025). (Test/net_test.R)

# Reference

- Cai, T. T., & Liu, W. (2016). Large-Scale Multiple Testing of Correlations. Journal of the American Statistical Association, 111(513), 229–240. https://doi.org/10.1080/01621459.2014.999157
- Fischer, E. K., Song, Y., Zhou, W., & Hoke, K. L. (2025). Flexibility in gene coexpression at developmental and evolutionary timescales. Molecular Biology and Evolution, msaf194. https://doi.org/10.1093/molbev/msaf194
- Maugis, P. A. G., Olhede, S. C., & Wolfe, P. J. (2017). Topology reveals universal features for network comparison. arXiv preprint arXiv:1705.05677.
- Shao, M., Xia, D., Zhang, Y., Wu, Q., & Chen, S. (2025). Higher-Order Accurate Two-Sample Network Inference and Network Hashing. Journal of the American Statistical Association, 1–13. https://doi.org/10.1080/01621459.2025.2520459
