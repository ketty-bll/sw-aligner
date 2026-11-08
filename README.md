# Smith-Waterman Local Aligner
 
A command-line tool, written in Python, that finds the best **local alignment** between two DNA, RNA or protein sequences using the Smith-Waterman algorithm. It was a project for my master's degree in quantitative and computational biology.
 
## What it does
 
The Smith-Waterman algorithm compares two sequences and finds the most similar region shared by both, even if the rest of the sequences is very different. This is useful for finding conserved domains, motifs or related genes.
 
The program:
1. Builds a substitution matrix (match / mismatch scores, or BLOSUM62 for proteins)
2. Fills a scoring matrix and a path matrix using dynamic programming
3. Traces back from the highest score to recover the best alignment
4. Prints the alignment with the number of matches, mismatches and gaps
## Installation
 
You need Python 3.8 or newer.
 
```bash
git clone https://github.com/ketty-bll/sw-aligner.git
cd sw-aligner
pip install -r requirements.txt
```
 
## Usage
 
```bash
python smith_waterman.py -s1 SEQUENCE1 -s2 SEQUENCE2 [options]
```
 
| Option | Description | Default |
|---|---|---|
| `-s1`, `--sequence1` | First sequence (required) | |
| `-s2`, `--sequence2` | Second sequence (required) | |
| `-a`, `--alphabet` | `DNA`, `RNA`, `AA` (proteins, simple scores) or `BLOSUM62` (proteins, BLOSUM62 matrix downloaded from NCBI) | `DNA` |
| `-match` | Match score (positive integer) | `3` |
| `-mismatch` | Mismatch score (negative integer) | `-3` |
| `-gap` | Penalty for each gap position (negative integer) | `-2` |
| `-v`, `--verbosity` | Also print the substitution and scoring matrices | off |
 
Sequences are not case sensitive. The `-match` and `-mismatch` options are ignored with `BLOSUM62`, which downloads the matrix and therefore needs an internet connection.
 
## Example
 
```bash
python smith_waterman.py \
  -s1 TGCATCGTACGTAGCTAGCTAGGATCGTAGCTAGCATCGTACGTAGC \
  -s2 ACGTACGTGACGTTAGCTAGCTAGCATCGATCGTACGTAGCTAGCTAG
```
 
Output (`|` = match, `:` = mismatch, space = gap):
 
```
A-CGTACGTGACGTTAGCTA-G--C-TAGCATCG-ATCGTACGTAGC
| ||||||| | |:|||||| |  | |||| |:| ||||||||||||
ATCGTACGT-A-GCTAGCTAGGATCGTAGC-TAGCATCGTACGTAGC
Sequence length: 47
Matches: 36 (76.6%);
Mismatches: 2 (4.26%);
Gaps: 9 (19.15%)
Score 84
```
 
Protein example:
 
```bash
python smith_waterman.py -s1 HEAGAWGHEE -s2 PAWHEAE -a AA
```
 
```
AW-HEAE
|| || |
AWGHE-E
Sequence length: 7
Matches: 5 (71.43%);
Mismatches: 0 (0.0%);
Gaps: 2 (28.57%)
Score 11
```
 
If the two sequences have no positive-scoring region in common, the program says so instead of printing an alignment.
 
## Limitations
 
- Uses a **linear gap penalty** (every gap position costs the same). Affine gap penalties are not implemented.
- Reports only the single best local alignment.
- Meant for short sequences and learning purposes. For large-scale searches, use optimised tools such as BLAST.
## Author
 
Ketty Bello
