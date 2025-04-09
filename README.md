# CS175 Midterm 2 Review Notes

## Regulatory Motif Finding

**Goal:** Find short, recurring patterns (motifs) in DNA sequences.

**Inputs:**
* `t`: Number of sequences
* `n`: Length of sequences
* `L`: Length of motif

**Key Concepts:**
* **Profile Matrix:** `4 x L` matrix of base frequencies at each motif position.
* **Consensus String:** `L`-mer from the most frequent base at each position.
* **Consensus Score:** Measures motif conservation.
* **Motif Positions `S = (s1, ..., st)`:** Starting position of the motif in each sequence. Total `(n - L + 1)^t` possibilities.

**Problems:**

1.  **Motif Finding Problem:** Find `L`-mers (one from each sequence) maximizing the Consensus Score.
2.  **Median String Finding Problem:** Find an `L`-mer minimizing the total Hamming distance to the closest `L`-mer in each sequence.

**Algorithms:**

* **Brute-Force Motif Search:**
    * Checks all `(n - L + 1)^t` position sets.
    * Finds max score.
    * *Complexity:* `O((n - L + 1)^t * t * L)` - Impractical.
* **Brute-Force Median String Search:**
    * Checks all `4^L` possible `L`-mers.
    * Finds min total distance.
    * *Complexity:* `O(4^L * t * n * L)` - Feasible for small `L`.
    * *Post-processing:* Find motif `L`-mers closest to the median string.
* **Tree-Based Search (Median/Motif):**
    * Explore search space (`4^L` or `(n-L+1)^t`) using a tree.
    * **Branch and Bound:** Prune branches using optimistic bounds to speed up search.
* **Greedy Motif Search:**
    * Iteratively builds profile, adding best match from next sequence.
    * Faster, but potentially suboptimal.
    * *Complexity:* `~O(L * n^2 * t)`.

---

## Sequence Alignment

**Goal:** Quantify similarity by aligning sequences.

**Key Concepts:**
* **Alignment:** Insert gaps (`-`) to align characters in columns.
* **Operations:** Substitution (Mismatch), Insertion, Deletion.
* **Edit Distance:** Min number of operations to transform one string to another.
* **Similarity Score:** Score based on matches (reward) and mismatches/gaps (penalty). Maximize this score.

**Algorithms:**

1.  **Global Alignment (Needleman-Wunsch):**
    * Aligns *entire* sequences.
    * *Method:* Dynamic Programming table `V[i, j]`.
    * *Traceback:* From `V[n, m]` to `V[0, 0]` for alignment.
    * *Complexity:* `O(nm)` time, `O(nm)` space (can be `O(min(n, m))` space for score only).

2.  **Variations (Global):**
    * **LCS (Longest Common Subsequence):** Score: Match=1, Mismatch=-inf, Gap=0. *Complexity:* `O(nm)`.
    * **Hamming Distance:** Equal length strings, no gaps. Score: Match=0, Mismatch=1. *Complexity:* `O(n)`.

3.  **Semi-Global Alignment:**
    * No penalty for leading/trailing gaps.
    * Useful for finding overlaps or fitting short sequence in long one.

4.  **Local Alignment (Smith-Waterman):**
    * Finds best alignment between *substrings*.
    * *Method:* Dynamic Programming. Key difference: `V[i, j] = max(0, ...)` allows new alignments to start anywhere.
    * *Traceback:* Start from the *highest score* in the matrix, end at `0`.
    * *Complexity:* `O(nm)` time, `O(nm)` space.

---

## Gap Penalty Models

* **Linear:** `Cost = length * penalty`. Simple.
* **Affine:** `Cost = open_penalty + length * extension_penalty`. More realistic (penalizes starting gaps more).
* **Convex:** Extension penalty decreases with gap length. Complex.

---

## Similarity Score Matrices (Proteins)

Score amino acid substitutions based on observed frequencies.

* **PAM (Point Accepted Mutation):**
    * Based on evolutionary model (mutations per 100 residues).
    * *Low PAM (PAM30):* Close sequences.
    * *High PAM (PAM250):* Distant sequences.
* **BLOSUM (Blocks Substitution Matrix):**
    * Based on observed substitutions in conserved blocks.
    * BLOSUM`X`: Derived from blocks with >= `X`% identity.
    * *High BLOSUM (BLOSUM80):* Close sequences.
    * *Low BLOSUM (BLOSUM45):* Distant sequences.
    * *Common Default:* BLOSUM62.
    * Calculated using log-odds scores of observed vs. expected frequencies.

