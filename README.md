# Dominant Eigenvalue-Eigenvector Pair Estimation via Graph Infection

We present a novel method to estimate the dominant eigenvalue and eigenvector pair of any non-negative real matrix via graph infection.

The algorithm only requires a few lines of code:

```python
#!/usr/bin/env python3

for i in range(num_steps):
    x_new = x_old + A @ x_old * Delta_t
    I_new = np.sum(x_new)

    m = (math.log(I_new) - math.log(I_old)) / Delta_t
    dominant_eigenval = (math.exp(m * Delta_t) - 1) / Delta_t

    x_old = x_new
    I_old = I_new

# DONE
```

## Set up and Run

#### Python
```bash
pip install numpy networkx
python3 graph_infection_method.py
```

#### R Lang
```bash
Rscript graph_infection_method.r
```

## Publication

Research paper accepted by **[Proc. 16th International Conference on Graph Transformation (ICGT 2023), Leicester, UK](https://conf.researchr.org/details/icgt-2023/icgt-2023-papers/4/Dominant-Eigenvalue-Eigenvector-Pair-Estimation-via-Graph-Infection)**.

Extended abstract accepted by the Graph Signal Processing (GSP) Workshop 2023, Oxford, UK

Please see our arvix paper for more details: 
[Arvix https://arxiv.org/abs/2208.00982](https://arxiv.org/abs/2208.00982)

![This is an image](./fig_overview.png)

### Citation

```latex
@article{yangxia2022dominant,
  title={Dominant Eigenvalue-Eigenvector Pair Estimation via Graph Infection},
  author={Yang, Kaiyuan and Xia, Li and Tay, YC},
  journal={arXiv preprint arXiv:2208.00982},
  year={2022}
}
```
