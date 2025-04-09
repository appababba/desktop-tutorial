# CS175 Midterm 2 Review Notes

## Table of Contents
- [Regulatory Motif Finding](#regulatory-motif-finding)
- [Sequence Alignment](#sequence-alignment)
- [Gap Penalty Models](#gap-penalty-models)
- [Similarity Score Matrices (Proteins)](#similarity-score-matrices-proteins)
- [Bioinformatics Algorithms: Pseudocode Summary](#bioinformatics-algorithms-pseudocode-summary)
  - [Brute-force Motif Search](#brute-force-motif-search)
  - [Brute-force Median String Search](#brute-force-median-string-search)
  - [Branch-and-Bound Median String Search](#branch-and-bound-median-string-search)
  - [Greedy Motif Search](#greedy-motif-search)
  - [Needleman-Wunsch (Global Alignment)](#needleman-wunsch-global-alignment)
  - [Smith-Waterman (Local Alignment)](#smith-waterman-local-alignment)

    
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



# Bioinformatics Algorithms: Pseudocode Summary

This document provides pseudocode summaries for common motif finding and sequence alignment algorithms.

## Motif Finding Algorithms

### Brute-force Motif Search

* **Goal:** Find the set of starting positions (one per sequence) for L-mers that maximizes the consensus score.
* **Input:** A list of `t` DNA sequences (`DNA`), motif length `L`.
* **Output:** `best_motif`, a tuple `(s1, s2, ..., st)` of starting positions.

**Algorithm:**

1.  Initialize `best_score` to 0.
2.  Initialize `best_motif` (the set of starting positions) to null or `(1, 1, ..., 1)`.
3.  For each possible combination of starting positions `S = (s1, s2, ..., st)` (where each `si` goes from 1 to `n-L+1`, and `n` is the length of the sequences):
    * Calculate `current_score = Score(S, DNA)`. This score is based on the consensus of the L-mers starting at the positions specified in `S`.
    * If `current_score > best_score`:
        * Set `best_score = current_score`.
        * Set `best_motif = S`.
4.  Return `best_motif`.

---

### Brute-force Median String Search

* **Goal:** Find the L-mer (`word`) that minimizes the total Hamming distance to the closest L-mer in each of the input sequences.
* **Input:** A list of `t` DNA sequences (`DNA`), motif length `L`.
* **Output:** `best_word`, the L-mer that acts as the median string.

**Algorithm:**

1.  Initialize `best_distance` to infinity.
2.  Initialize `best_word` to the first possible L-mer (e.g., "AAA...A").
3.  For each possible L-mer `word` (from "AAA...A" to "TTT...T"):
    * Calculate `current_distance = TotalDistance(word, DNA)`.
        * *(Helper: To calculate `TotalDistance(word, DNA)`: For each sequence `i` in `DNA`, find the minimum Hamming distance `d(word, Seq_i)` between `word` and all possible L-mers within `Seq_i`. Sum these minimum distances across all sequences).*
    * If `current_distance < best_distance`:
        * Set `best_distance = current_distance`.
        * Set `best_word = word`.
4.  Return `best_word`.

---

### Branch-and-Bound Median String Search

* **Goal:** Find the median string more efficiently than brute-force by pruning the search space. Uses a conceptual tree of all possible L-mers.
* **Input:** A list of `t` DNA sequences (`DNA`), motif length `L`.
* **Output:** `best_word`, the L-mer that acts as the median string.

**Algorithm:**

1.  Initialize `best_distance` to infinity, `best_word` to null.
2.  Initialize tree traversal variables (e.g., current prefix array `s`, current level `i`). Start traversal (e.g., at the root or first level).
3.  While the traversal is not finished (e.g., using the logic `i > 0` from the reference slides):
    * **If the current node is not a leaf (`i < L`):**
        * `prefix` = string represented by the current node `s` up to level `i`.
        * Calculate `optimistic_distance = OptimisticDistance(prefix, DNA)` (sum of minimum distances of the `prefix` to each sequence, considering potential completions).
        * If `optimistic_distance > best_distance`:
            * Prune this branch: Move to the next node using a `Bypass(s, i, L, 4)` logic (find the next node at the same or higher level, skipping children).
        * Else (`optimistic_distance <= best_distance`):
            * Continue deeper: Move to the next node using `NextVertex(s, i, L, 4)` logic (go to the first child or the next sibling if no children exist).
    * **Else (current node is a leaf, `i == L`):**
        * `word` = L-mer string represented by the leaf node `s`.
        * Calculate `current_distance = TotalDistance(word, DNA)`.
        * If `current_distance < best_distance`:
            * Set `best_distance = current_distance`.
            * Set `best_word = word`.
        * Move to the next node using `NextVertex(s, i, L, 4)` logic (move to the next leaf sibling or backtrack if no more siblings).
4.  Return `best_word`.

*(Note: `NextVertex` and `Bypass` refer to specific functions for traversing the conceptual L-mer tree, often implemented using array manipulation to represent the current path/prefix).*

---

### Greedy Motif Search

* **Goal:** Find a plausible motif (set of starting positions) by iteratively building the motif profile, making locally optimal choices at each step.
* **Input:** A list of `t` DNA sequences (`DNA`), motif length `L`.
* **Output:** `best_motif`, a tuple `(s1, s2, ..., st)` of starting positions.

**Algorithm:**

1.  Initialize `best_motif` `S = (1, 1, ..., 1)`.
2.  **(Find best motif for the first 2 sequences):**
    * Initialize `best_score_so_far = -1`.
    * For `s1` from 1 to `n-L+1`:
        * For `s2` from 1 to `n-L+1`:
            * `current_S = (s1, s2)`
            * `current_score = Score(current_S, first 2 sequences)`
            * If `current_score > best_score_so_far`:
                * `best_score_so_far = current_score`
                * `best_motif[1] = s1`
                * `best_motif[2] = s2`
3.  **(Add motifs for remaining sequences):**
    * For `i` from 3 to `t`:
        * Initialize `best_score_for_i = -1` (or use the score based on `best_si_for_i = 1`).
        * Initialize `best_si_for_i = 1`.
        * For `si` from 1 to `n-L+1`:
            * Form `current_S = (best_motif[1], ..., best_motif[i-1], si)`.
            * Calculate `current_score = Score(current_S, first i sequences)`.
            * Calculate `best_known_score_for_i = Score((best_motif[1], ..., best_motif[i-1], best_si_for_i), first i sequences)` (Score using the best `si` found *so far* for this sequence `i`).
            * If `current_score > best_known_score_for_i`:
                * *(Update the best score for this iteration)* `best_known_score_for_i = current_score`
                * *(Update the best starting position for sequence i)* `best_si_for_i = si`
        * Set `best_motif[i] = best_si_for_i`.
4.  Return `best_motif`.

*(Note: The `Score` function typically involves building a profile matrix from the L-mers defined by the current motif `S` and calculating a consensus score based on that profile).*

---

## Sequence Alignment Algorithms

### Needleman-Wunsch (Global Alignment)

* **Goal:** Find the optimal alignment between two sequences across their entire lengths.
* **Input:** Sequence `S` (length `n`), Sequence `T` (length `m`), scoring scheme (match score, mismatch penalty, gap penalty).
* **Output:** The optimal global alignment and its score.

**Algorithm:**

1.  Create a grid `V` of size `(n+1) x (m+1)`.
2.  **Initialization:**
    * `V[0, 0] = 0`
    * For `i` from 1 to `n`: `V[i, 0] = V[i-1, 0] + gap_penalty` (penalty for gaps in T)
    * For `j` from 1 to `m`: `V[0, j] = V[0, j-1] + gap_penalty` (penalty for gaps in S)
3.  **Fill Grid:**
    * For `i` from 1 to `n`:
        * For `j` from 1 to `m`:
            * `match_mismatch_score = V[i-1, j-1] + score(S[i] vs T[j])` (Score depends on whether `S[i]` and `T[j]` match or mismatch)
            * `delete_score = V[i-1, j] + gap_penalty` (Align `S[i]` with a gap in `T`)
            * `insert_score = V[i, j-1] + gap_penalty` (Align `T[j]` with a gap in `S`)
            * `V[i, j] = max(match_mismatch_score, delete_score, insert_score)`
            * *(Optional: Store pointers indicating which choice(s) led to the max score for traceback)*
4.  **Final Score:** The optimal global alignment score is `V[n, m]`.
5.  **Traceback:**
    * Start from cell `V[n, m]`.
    * Follow the pointers (or recalculate the path based on scores) back to `V[0, 0]`.
    * Each step corresponds to an alignment operation (match/mismatch, insertion, or deletion), constructing the alignment string(s).

---

### Smith-Waterman (Local Alignment)

* **Goal:** Find the highest-scoring alignment between *subsequences* of two sequences.
* **Input:** Sequence `S` (length `n`), Sequence `T` (length `m`), scoring scheme (match score, mismatch penalty, gap penalty - usually positive for match, negative for mismatch/gap).
* **Output:** The optimal local alignment(s) and the maximum score.

**Algorithm:**

1.  Create a grid `V` of size `(n+1) x (m+1)`.
2.  **Initialization:**
    * For `i` from 0 to `n`: `V[i, 0] = 0`
    * For `j` from 0 to `m`: `V[0, j] = 0`
3.  **Fill Grid:**
    * Initialize `max_score = 0`.
    * Initialize `max_i = 0`, `max_j = 0` (to store the coordinates of the cell with the `max_score`).
    * For `i` from 1 to `n`:
        * For `j` from 1 to `m`:
            * `match_mismatch_score = V[i-1, j-1] + score(S[i] vs T[j])`
            * `delete_score = V[i-1, j] + gap_penalty`
            * `insert_score = V[i, j-1] + gap_penalty`
            * `V[i, j] = max(match_mismatch_score, delete_score, insert_score, 0)` *(Crucially, include 0 to allow alignments to start/end anywhere)*
            * If `V[i, j] > max_score`:
                * `max_score = V[i, j]`
                * `max_i = i`
                * `max_j = j`
            * *(Optional: Store pointers)*
4.  **Final Score:** The maximum local alignment score is `max_score`.
5.  **Traceback:**
    * Start from the cell `(max_i, max_j)` containing `max_score`.
    * Follow pointers (or recalculate path) back until a cell with a score of `0` is reached.
    * This path corresponds to the best local alignment. (Multiple paths might exist if there are ties for `max_score`).



