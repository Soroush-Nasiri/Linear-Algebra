# PageRank and HITS Algorithms

This project contains C++ implementations of the PageRank and HITS (Hyperlink-Induced Topic Search) algorithms, applying them to a directed graph dataset (`wiki-Vote.txt`). 

## Project Overview

The main objective of this project is to analyze the importance of nodes within a network using linear algebra concepts. It reads the network structure, computes the PageRank scores, and calculates the Authority and Hub scores using the HITS algorithm iteratively.

To verify the correctness of the iterative HITS implementation, the project also exports the adjacency matrix and utilizes a Python script (`Help.py`) to compute the Singular Value Decomposition (SVD) of the matrix using `scipy.sparse.linalg.svds`. The principal singular vectors correspond to the Authority and Hub scores, which are then compared against the iteratively computed results.

## Requirements

*   A C++ compiler (e.g., `g++`)
*   Python 3
*   `numpy` and `scipy` Python libraries

## How to Run

1.  Compile the C++ code:
    ```bash
    g++ main.cpp
    ```
2.  Execute the compiled program:
    ```bash
    ./a.out
    ```

The program will output the PageRank, Authority, Hub, and Degree rankings, followed by a comparison with the SVD results.
