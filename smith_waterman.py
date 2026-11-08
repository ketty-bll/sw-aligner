#!/usr/bin/env python3

"""Smith-Waterman local sequence aligner.

Aligns two DNA, RNA or protein sequences with the Smith-Waterman algorithm
(local alignment, linear gap penalty) and prints the best local alignment.

Author: Ketty Bello

Example:
    python smith_waterman.py -s1 TGCATCGTACGTAGC -s2 ACGTACGTGACGTTAGC
"""

import argparse
import sys
import numpy as np
import pandas as pd

BLOSUM62_URL = "https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62"
ALPHABETS = {
    "DNA": "ACTG",
    "RNA": "ACUG",
    "AA": "ACDEFGHIKLMNPQRSTVWY",
}


# ---------------------------------------------------------------------------
# Command line arguments
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Smith-Waterman local aligner (DNA, RNA or protein)."
    )
    parser.add_argument("-s1", "--sequence1", required=True,
                        help="First sequence")
    parser.add_argument("-s2", "--sequence2", required=True,
                        help="Second sequence")
    parser.add_argument("-a", "--alphabet", default="DNA",
                        choices=["DNA", "RNA", "AA", "BLOSUM62"],
                        help="DNA, RNA, AA (proteins, simple match/mismatch "
                             "scores) or BLOSUM62 (proteins, BLOSUM62 matrix "
                             "downloaded from NCBI). Default: DNA")
    parser.add_argument("-v", "--verbosity", action="store_true",
                        help="also print the substitution and scoring matrices")
    parser.add_argument("-match", "--match_score", type=int, default=3,
                        help="positive integer, match score (default: 3)")
    parser.add_argument("-mismatch", "--mismatch_score", type=int, default=-3,
                        help="negative integer, mismatch score (default: -3)")
    parser.add_argument("-gap", "--gap_penalty", type=int, default=-2,
                        help="negative integer, penalty for each gap position "
                             "(default: -2)")
    args = parser.parse_args(argv)

    if args.match_score <= 0:
        parser.error("--match_score must be a positive integer")
    if args.mismatch_score >= 0:
        parser.error("--mismatch_score must be a negative integer")
    if args.gap_penalty >= 0:
        parser.error("--gap_penalty must be a negative integer")
    return args


# ---------------------------------------------------------------------------
# Substitution matrix
# ---------------------------------------------------------------------------
def create_substitution_matrix(alphabet, match_score, mismatch_score,
                               verbose=False):
    """Return a DataFrame with the score for every pair of letters."""
    if alphabet in ALPHABETS:
        letters = list(ALPHABETS[alphabet])
        values = np.full((len(letters), len(letters)), mismatch_score,
                         dtype=int)
        np.fill_diagonal(values, match_score)
        matrix = pd.DataFrame(values, index=letters, columns=letters)
        title = f"{alphabet} substitution matrix:"
    else:  # BLOSUM62
        try:
            import requests
            text = requests.get(BLOSUM62_URL, timeout=15).text
        except Exception as error:
            sys.exit(f"Could not download the BLOSUM62 matrix ({error}). "
                     "Check your internet connection or use -a AA instead.")
        lines = [l for l in text.splitlines() if l and not l.startswith("#")]
        letters = lines[0].split()
        values = [[int(x) for x in l.split()[1:]] for l in lines[1:]]
        matrix = pd.DataFrame(values, index=letters, columns=letters)
        title = "BLOSUM62 substitution matrix:"

    if verbose:
        print(title)
        print(matrix)
    return matrix


def check_sequence(sequence, matrix, name):
    """Make the sequence upper case and check every letter is allowed."""
    sequence = sequence.upper()
    allowed = set(matrix.index)
    bad = sorted(set(sequence) - allowed)
    if bad:
        sys.exit(f"{name} contains letters that are not valid for the chosen "
                 f"alphabet: {', '.join(bad)}")
    return list(sequence)


# ---------------------------------------------------------------------------
# Smith-Waterman algorithm
# ---------------------------------------------------------------------------
def calculate_score(score_M, r, c, S1, S2, substitution_M, gap_penalty):
    """Score of cell (r, c) and the direction it came from.

    Options: start over with 0 (ZERO), come from the diagonal (DIAG),
    from the cell above (UP) or from the cell on the left (LEFT).
    """
    pair_score = substitution_M.at[S1[r - 1], S2[c - 1]]
    diag_score = score_M[r - 1][c - 1] + pair_score
    up_score = score_M[r - 1][c] + gap_penalty
    left_score = score_M[r][c - 1] + gap_penalty

    scores = [0, diag_score, up_score, left_score]
    direction = ["ZERO", "DIAG", "UP", "LEFT"][int(np.argmax(scores))]
    return max(scores), direction


def create_matrices(S1, S2, substitution_M, gap_penalty):
    """Fill the scoring and path matrices.

    Returns score_M, path_M, max_score and max_pos (row, column of the
    highest score, where the best local alignment ends).
    """
    nrows, ncols = len(S1) + 1, len(S2) + 1
    score_M = np.zeros((nrows, ncols), dtype=int)
    path_M = np.full((nrows, ncols), "NULL", dtype=object)
    max_score, max_pos = -1000, None

    for r in range(1, nrows):
        for c in range(1, ncols):
            score, direction = calculate_score(score_M, r, c, S1, S2,
                                               substitution_M, gap_penalty)
            score_M[r][c] = score
            path_M[r][c] = direction
            if score >= max_score:
                max_score, max_pos = score, (r, c)

    return score_M, path_M, max_score, max_pos


def traceback(path_M, score_M, S1, S2, start_pos):
    """Walk back from the best cell until the score reaches 0."""
    aligned_S1, aligned_S2 = [], []
    r, c = start_pos

    while path_M[r][c] != "ZERO" and score_M[r][c] != 0:
        direction = path_M[r][c]
        if direction == "DIAG":
            aligned_S1.append(S1[r - 1])
            aligned_S2.append(S2[c - 1])
            r -= 1
            c -= 1
        elif direction == "UP":
            aligned_S1.append(S1[r - 1])
            aligned_S2.append("-")
            r -= 1
        elif direction == "LEFT":
            aligned_S1.append("-")
            aligned_S2.append(S2[c - 1])
            c -= 1
        else:
            raise RuntimeError(f"Unexpected value in path matrix at {r, c}")

    return "".join(reversed(aligned_S1)), "".join(reversed(aligned_S2))


# ---------------------------------------------------------------------------
# Printing
# ---------------------------------------------------------------------------
def print_matrix(score_M, S1, S2):
    row_labels = [" "] + S1
    col_labels = [" "] + S2
    print("Scoring matrix:")
    print(pd.DataFrame(score_M, index=row_labels, columns=col_labels))


def print_alignment(aligned_S1, aligned_S2):
    """Print the alignment (sequence 2 on top) and match/mismatch/gap stats."""
    middle = []
    n_matches = n_mismatches = n_gaps = 0
    for a, b in zip(aligned_S1, aligned_S2):
        if a == b:
            middle.append("|")
            n_matches += 1
        elif a != "-" and b != "-":
            middle.append(":")
            n_mismatches += 1
        else:
            middle.append(" ")
            n_gaps += 1

    length = len(aligned_S1)
    print(aligned_S2)
    print("".join(middle))
    print(aligned_S1)
    print(f"Sequence length: {length}")
    print(f"Matches: {n_matches} ({round(100 * n_matches / length, 2)}%);\n"
          f"Mismatches: {n_mismatches} ({round(100 * n_mismatches / length, 2)}%);\n"
          f"Gaps: {n_gaps} ({round(100 * n_gaps / length, 2)}%)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    substitution_M = create_substitution_matrix(
        args.alphabet, args.match_score, args.mismatch_score, args.verbosity)
    S1 = check_sequence(args.sequence1, substitution_M, "Sequence 1")
    S2 = check_sequence(args.sequence2, substitution_M, "Sequence 2")

    score_M, path_M, max_score, max_pos = create_matrices(
        S1, S2, substitution_M, args.gap_penalty)

    if args.verbosity:
        print_matrix(score_M, S1, S2)

    if max_score <= 0:
        print("No local alignment with a positive score was found.")
        return

    aligned_S1, aligned_S2 = traceback(path_M, score_M, S1, S2, max_pos)
    print_alignment(aligned_S1, aligned_S2)
    print(f"Score {max_score}")


if __name__ == "__main__":
    main()
