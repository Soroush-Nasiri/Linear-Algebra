# Image Compression using SVD

This project demonstrates low-rank image compression using **Singular Value Decomposition (SVD)** in Python. 

An image matrix $A \in \mathbb{R}^{m \times n}$ is loaded, converted to grayscale, and approximated using its top $r$ singular values:

$$A_r = U_r \Sigma_r V_r^T = \sum_{i=1}^r \sigma_i u_i v_i^T$$

### Compression Rationale
- **Original storage**: $m \times n$ values
- **Rank-$r$ storage**: $r(m + n + 1)$ values
- **Compression Ratio**: $\frac{r(m + n + 1)}{m \times n}$

---

## Results Summary

Tested on a $2961 \times 3000$ image:

| Rank ($r$) | Storage (% of original) | Energy Captured (%) | Rel Frobenius Error (%) | Quality |
| :---: | :---: | :---: | :---: | :--- |
| **$r = 10$** | 0.67% | 95.65% | 20.80% | Basic outlines visible |
| **$r = 50$** | 3.36% | 98.69% | 11.40% | Good visual fidelity |
| **$r = 100$** | 6.71% | 99.39% | 7.81% | Near-original quality with >93% storage reduction |

