from typing import Dict, Tuple, List, Set

Direction = str  # 'D', 'U', 'L'

def load_substitution_matrix(path: str) -> Dict[Tuple[str, str], float]:
    """
   file format:

    # comment baij bolno
       A  C  G  T
    A  5 -1 -2 -1
    C -1  5 -3 -2
    G -2 -3  5 -2
    T -1 -2 -2  5
    """
    matrix: Dict[Tuple[str, str], float] = {}
    with open(path, encoding="utf-8") as f:
        lines = [line.strip() for line in f
                 if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError("Hooson file.")

    header = lines[0].split()
    for line in lines[1:]:
        parts = line.split()
        row_char = parts[0]
        scores = [float(x) for x in parts[1:]]
        if len(scores) != len(header):
            raise ValueError("Matrixiin mur baganiin tooo.")
        for col_char, score in zip(header, scores):
            matrix[(row_char, col_char)] = score
    return matrix


def create_simple_dna_matrix(match: float = 1.0,
                             mismatch: float = -1.0) -> Dict[Tuple[str, str], float]:
    """engiin DNA matrix (A,C,G,T)"""
    bases = "ACGT"
    matrix: Dict[Tuple[str, str], float] = {}
    for a in bases:
        for b in bases:
            matrix[(a, b)] = match if a == b else mismatch
    return matrix


def create_simple_protein_matrix(match: float = 1.0,
                                 mismatch: float = -1.0) -> Dict[Tuple[str, str], float]:
    """Engiin uurgiin matrix (20 amino acid)"""
    amino_acids = "ARNDCQEGHILKMFPSTWYV"
    matrix: Dict[Tuple[str, str], float] = {}
    for a in amino_acids:
        for b in amino_acids:
            matrix[(a, b)] = match if a == b else mismatch
    return matrix


def get_score(a: str, b: str,
              matrix: Dict[Tuple[str, str], float],
              default: float = -1.0) -> float:
    """Onoog unshij avah baihgui bol defualt."""
    if (a, b) in matrix:
        return matrix[(a, b)]
    if (b, a) in matrix:
        return matrix[(b, a)]
    return default


def build_dp(seq1: str,
             seq2: str,
             matrix: Dict[Tuple[str, str], float],
             gap_penalty: float = -2.0,
             alignment_type: str = "global"):
    """
    DP-n onooni bolon zamiin matrix
    alignment_type:
        - 'global' : Needleman–Wunsch
        - 'local'  : Smith–Waterman
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n, m = len(seq1), len(seq2)

    # score[i][j] – i urttai seq1, j urttai seq2 n hamgiin sain onoo
    score = [[0.0] * (m + 1) for _ in range(n + 1)]
    trace: List[List[List[Direction]]] = [[[] for _ in range(m + 1)] for _ in range(n + 1)]

    max_score = float("-inf")
    max_positions: List[Tuple[int, int]] = []

    if alignment_type == "global":
        for i in range(1, n + 1):
            score[i][0] = i * gap_penalty
            trace[i][0] = ["U"]  # deerees(gap seq2)
        for j in range(1, m + 1):
            score[0][j] = j * gap_penalty
            trace[0][j] = ["L"]  # zuunees(gap seq1)
        max_score = score[n][m]
        max_positions = [(n, m)]

    elif alignment_type == "local":
        max_score = 0.0
        max_positions = []
    else:
        raise ValueError("alignment_type global esvel local baih yostoi.")

    # DP
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = get_score(seq1[i - 1], seq2[j - 1], matrix)
            diag = score[i - 1][j - 1] + match_score
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty

            if alignment_type == "local":
                cell_score = max(0.0, diag, up, left)
            else:
                cell_score = max(diag, up, left)

            score[i][j] = cell_score
            trace[i][j] = []

            if alignment_type == "local" and cell_score == 0.0:
                continue

            if cell_score == diag:
                trace[i][j].append("D")
            if cell_score == up:
                trace[i][j].append("U")
            if cell_score == left:
                trace[i][j].append("L")

            if alignment_type == "local":
                if cell_score > max_score:
                    max_score = cell_score
                    max_positions = [(i, j)]
                elif cell_score == max_score and cell_score > 0:
                    max_positions.append((i, j))

    if alignment_type == "global":
        max_score = score[n][m]
        max_positions = [(n, m)]

    return score, trace, max_score, max_positions


def traceback_alignments(seq1: str,
                         seq2: str,
                         score,
                         trace,
                         start_positions: List[Tuple[int, int]],
                         alignment_type: str = "global",
                         max_alignments: int = 5):
    """
    DP matrix-ees alignments-g butsaana.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    results: List[Tuple[str, str, float]] = []
    seen: Set[Tuple[str, str]] = set()

    def backtrack(i: int, j: int, aln1: List[str], aln2: List[str]):
        if len(results) >= max_alignments:
            return

        if alignment_type == "global":
            if i == 0 and j == 0:
                a1 = "".join(reversed(aln1))
                a2 = "".join(reversed(aln2))
                key = (a1, a2)
                if key not in seen:
                    seen.add(key)
                    results.append((a1, a2, score[len(seq1)][len(seq2)]))
                return
        else:  # local
            if score[i][j] == 0:
                a1 = "".join(reversed(aln1))
                a2 = "".join(reversed(aln2))
                key = (a1, a2)
                if key not in seen:
                    seen.add(key)
                    results.append((a1, a2, score[i][j]))
                return

        for d in trace[i][j]:
            if d == "D":
                backtrack(i - 1, j - 1,
                          aln1 + [seq1[i - 1]],
                          aln2 + [seq2[j - 1]])
            elif d == "U":
                backtrack(i - 1, j,
                          aln1 + [seq1[i - 1]],
                          aln2 + ["-"])
            elif d == "L":
                backtrack(i, j - 1,
                          aln1 + ["-"],
                          aln2 + [seq2[j - 1]])

    for (i, j) in start_positions:
        backtrack(i, j, [], [])

    return results


def print_matrix(mat, seq1: str, seq2: str, title: str = "Score matrix"):
    """onooni esvel zamiig hevleh"""
    seq1 = " " + seq1.upper()
    seq2 = " " + seq2.upper()
    print(f"\n{title}")
    header = "    " + "  ".join(seq2)
    print(header)
    for i, row in enumerate(mat):
        line = seq1[i] + " "
        for val in row:
            if isinstance(val, list):         # traceback mtrix
                cell = "".join(val) if val else "."
            else:                             # score matrix
                cell = f"{val:4.0f}"
            line += f"{cell:>4}"
        print(line)


def main():
    seq_type = input("Protein or dna (protein/dna): ").strip().lower()
    if seq_type not in ("protein", "dna"):
        print("Error: 'protein' esvel 'dna' gej oruulna uu.")
        return

    seq1 = input("1st (seq1): ").strip().upper()
    seq2 = input("2nd (seq2): ").strip().upper()

    alignment_type = input("Global or dna (global/local): ").strip().lower()
    if alignment_type not in ("global", "local"):
        print("Error: 'global' esvel 'local' gej oruulna uu.")
        return

    if seq_type == "protein":
        print(" Score matrix:")
        print("   1) BLOSUM62")
        print("   2) PAM250 ")
        print("   3) simple  (match=1, mismatch=-1)")
        choice = input("Choose (1/2/3): ").strip()
        if choice in ("1", "2"):
            path = input("Path tohiruulah (ex: BLOSUM62.txt): ").strip()
            subst_matrix = load_substitution_matrix(path)
        else:
            subst_matrix = create_simple_protein_matrix()
    else:  # dna
        print("DNA score matrix:")
        print("   1) custom ")
        print("   2) simple (match=1, mismatch=-1)")
        choice = input("choose (1/2): ").strip()
        if choice == "1":
            path = input("path tohiruuulah (ex: DNA_MATRIX.txt): ").strip()
            subst_matrix = load_substitution_matrix(path)
        else:
            subst_matrix = create_simple_dna_matrix()

    try:
        gap_penalty = float(
            input("gap baival (ex: -2): ").strip() or "-2"
        )
    except ValueError:
        print("buruu utga oruullaa. defualtaar -2")
        gap_penalty = -2.0

    score, trace, max_score, max_positions = build_dp(
        seq1, seq2, subst_matrix,
        gap_penalty=gap_penalty,
        alignment_type=alignment_type
    )

    print_matrix(score, seq1, seq2, title="Score matrix")
    print_matrix(trace, seq1, seq2, title="Traceback matrix (D=diag, U=up, L=left)")

    alignments = traceback_alignments(
        seq1, seq2, score, trace,
        max_positions,
        alignment_type=alignment_type,
        max_alignments=5
    )
    if not alignments:
        print("Haritsuulalm oldsongui.")
        return

    if alignment_type == "global":
        print(f"result(global): {score[len(seq1)][len(seq2)]:.1f}")
    else:
        print(f"result(local): {max_score:.1f}")

    for idx, (a1, a2, sc) in enumerate(alignments, start=1):
        print(f"\n--- Haritsuulalm {idx} ---")
        print(a1)
        print(a2)


if __name__ == "__main__":
    main()
