import pyboolector
from pyboolector import Boolector
from utils import rotate_right, modular_addition

def get_diff(filename):
    diffs = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            # Ignore comments and blank lines.
            if not line or line.startswith("#"):
                continue
            a, b, c, d = line.split(',')
            diffs.append((int(a, 16), int(b, 16), int(c, 16), int(d, 16)))
    return diffs

# with considering key-schedule
def cham_quasidifferential_trails(diffs, nb_bits, key_words):
    """Build the extended-quasidifferential SMT model for CHAM."""
    btor = Boolector()
    btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)
    btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)

    nb_rounds = len(diffs) - 1

    # data : u[round_idx][u_idx]
    u = [[btor.Var(btor.BitVecSort(nb_bits), "u%d_%d" % (ui, i)) for ui in range(4)] for i in range(nb_rounds + 1)]
    # key : key masks depend on data masks only -> no more Variables
    rk_mask = [btor.Const(0, nb_bits)] * (2 * key_words)
    
    weight = btor.Const(0, nb_bits)
    for i in range(nb_rounds):
        # 1. data path
        if i % 2 == 1:  # odd
            alpha, beta = 1, 8
        else:  # even
            alpha, beta = 8, 1

        a = btor.Const(diffs[i][0], nb_bits)
        b = btor.Const(diffs[i][1], nb_bits)
        c = btor.Const(diffs[i + 1][3], nb_bits)

        b = btor.Rol(b, beta)
        c = btor.Ror(c, alpha)

        u_ = u[i][0]
        v_ = btor.Rol(u[i][1] ^ u[i + 1][0], beta)
        w_ = btor.Ror(u[i + 1][3], alpha)

        weight += modular_addition(
            btor, a, b, c, u_, v_, w_, nb_bits, i
        )

        # other branch
        btor.Assert(u[i][2] == u[i + 1][1])
        btor.Assert(u[i][3] == u[i + 1][2])

        # Accumulate the data-path mask entering RK[i mod 2m].
        rk_mask[i % (2 * key_words)] += v_

    # Enforce the transpose of CHAM's linear key expansion (paper, Eq. (5)).
    for i in range(key_words):
        v1 = rk_mask[i]
        v2 = rk_mask[(i + key_words) ^ 1]

        p = btor.Ror(v1, 8) ^ (v1 ^ v2) ^ btor.Ror(v1 ^ v2, 1) ^ btor.Ror(v1, 11)
        btor.Assert(p == 0)

    # end point
    for ui in range(4):
        btor.Assert(u[0][ui] == 0)
        btor.Assert(u[nb_rounds][ui] == 0)
        
    return btor, weight

def distinctness_constraint(btor, masks, solutions, word_size):
    distinctness_condition = btor.Const(1)
    for solution in solutions:
        condition = btor.Const(0)
        for i in range(len(masks)):  # round
            for j in range(len(masks[i])):  # word = 4
                condition |= (masks[i][j] != btor.Const(solution[i][j], word_size))
        distinctness_condition &= condition
    btor.Assume(distinctness_condition)

def solve_all(btor, weight, w, nb_rounds, word_size, key_words):
    # Get variables
    u = [[btor.Match_by_symbol("u%d_%d" % (ui, i)) for ui in range(4)] for i in range(nb_rounds + 1)]
    
    solutions = []
    while True:
        btor.Assume(weight == w)
        distinctness_constraint(btor, u, solutions, word_size)

        r = btor.Sat()
        if r != btor.SAT:
            return solutions

        solutions.append([
            (
                int(u[i][0].assignment, base=2),
                int(u[i][1].assignment, base=2),
                int(u[i][2].assignment, base=2),
                int(u[i][3].assignment, base=2)
            )
            for i in range(nb_rounds + 1)
        ])

def parity(x):
    return bin(x).count('1') % 2

def rotl(x, r, word_size):
   mask = (1 << word_size) - 1
   return ((x << r) & mask) | ((x & mask) >> (word_size - r))

def compute_sign_cham(differences, masks, word_size, key_words):
    """Compute the correlation sign of a CHAM quasidifferential trail."""

    def pseudoinverseM(t):
        t = t ^ ((t << 1) % 2 ** word_size)
        return t >> 1

    def complement(t):
        return (2 ** word_size - 1) ^ t

    s = 1
    for i in range(len(masks) - 1):
        if i % 2 == 1:  # odd
            alpha, beta = 1, 8
        else:  # even
            alpha, beta = 8, 1

        u = masks[i][0]
        v = rotl(masks[i][1] ^ masks[i + 1][0], beta, word_size)
        w = rotl(masks[i + 1][3], word_size - alpha, word_size)
        (u_, v_, _) = (u ^ w, v ^ w, u ^ v ^ w)

        (a, b_in, _, __) = differences[i]
        c_in             = differences[i + 1][3]

        b = rotl(b_in, beta, word_size)
        c = rotl(c_in, word_size - alpha, word_size)

        (a_, b_, c_) = (b ^ c, a ^ c, pseudoinverseM(a ^ b ^ c))
        p1 = parity(((complement(a_) & u_) ^ (c_ & v_)) & ((complement(b_) & v_) ^ (c_ & u_)))
        p2 = parity((u_ & v_) & (c_ ^ (a_ & b_ & complement(c_))))
        s *= (-1) ** (p1 ^ p2)

        # Walsh-character sign introduced by XORing CHAM's round constant.
        rc = i
        s *= (-1) ** (bin(u & rc).count("1"))
    return s
