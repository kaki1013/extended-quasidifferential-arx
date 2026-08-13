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
            a, b = line.split(',')
            diffs.append((int(a, 16), int(b, 16)))
    return diffs

def speck_1r(btor, diffs, nb_bits, data_mask, key_mask, weight, i, rk_mask, key_words):
    """Encode a SPECK round and accumulate its mask into RK[i mod 2m]."""
    u, v = data_mask

    a = btor.Const(diffs[i][0], nb_bits)
    b = btor.Const(diffs[i][1], nb_bits)
    c = btor.Const(diffs[i + 1][0], nb_bits)

    # data path
    if nb_bits == 16:
        u_ = btor.Ror(u[i], 7)
        v_ = btor.Ror(v[i + 1], 2) ^ v[i] 
        a  = btor.Ror(a, 7)
    else:
        u_ = rotate_right(btor, u[i], 8, nb_bits)
        v_ = rotate_right(btor, v[i + 1], 3, nb_bits) ^ v[i] 
        a  = rotate_right(btor, a, 8, nb_bits)
    w_ = u[i + 1] ^ v[i + 1]
    weight += modular_addition(
        btor, a, b, c, u_, v_, w_, nb_bits, i
    )

    # key path
    rk_mask[i % (2 * key_words)] += w_

    return weight

# with considering key-schedule
def speck_quasidifferential_trails(diffs, nb_bits, key_words):
    """Build the SMT model for SPECK with the CHAM-like linear key schedule."""
    btor = Boolector()
    btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)
    btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)

    nb_rounds = len(diffs) - 1

    # data
    u = [btor.Var(btor.BitVecSort(nb_bits), "u%d" % i) for i in range(nb_rounds + 1)]
    v = [btor.Var(btor.BitVecSort(nb_bits), "v%d" % i) for i in range(nb_rounds + 1)]
    # key (q : mask for the branch that key addition arises)
    p = [[btor.Var(btor.BitVecSort(nb_bits), "p%d_%d" % (mi, i)) for i in range(nb_rounds + 1)] for mi in range(key_words-1)]
    q = [btor.Var(btor.BitVecSort(nb_bits), "q%d" % i) for i in range(nb_rounds + 1)]

    # masks for cham-like key schedule
    rk_mask = [btor.Const(0, nb_bits)] * (2 * key_words)
    
    weight = btor.Const(0, nb_bits)
    data_mask = [u, v]
    key_mask = [p, q]

    # 0~ rounds
    for i in range(nb_rounds):
        # data & key
        weight = speck_1r(btor, diffs, nb_bits, data_mask, key_mask, weight, i, rk_mask, key_words)

    # Enforce the transpose of CHAM's linear key expansion (paper, Eq. (5)).
    for i in range(key_words):
        v1 = rk_mask[i]
        v2 = rk_mask[(i + key_words) ^ 1]

        # tmp = btor.Ror(v1, 8) ^ (v1 ^ v2) ^ btor.Ror(v1 ^ v2, 1) ^ btor.Ror(v1, 11)
        tmp = rotate_right(btor, v1, 8, nb_bits) ^ (v1 ^ v2) ^ rotate_right(btor, v1 ^ v2, 1, nb_bits) ^ rotate_right(btor, v1, 11, nb_bits)
        btor.Assert(tmp == 0)

    # key
    btor.Assert(q[0] == 0);             btor.Assert(q[nb_rounds] == 0)
    for mi in range(key_words-1):
        btor.Assert(p[mi][0] == 0);     btor.Assert(p[mi][nb_rounds] == 0)
    # cham-like key schedule doesn't need key mask
    for ri in range(nb_rounds+1):
        for mi in range(key_words-1):
            btor.Assert(p[mi][ri] == 0)
        btor.Assert(q[ri] == 0)

    # data
    btor.Assert(u[0] == 0);             btor.Assert(v[0] == 0)
    btor.Assert(u[nb_rounds] == 0);     btor.Assert(v[nb_rounds] == 0)

    return btor, weight

def distinctness_constraint(btor, masks, solutions, word_size, key_words):
    distinctness_condition = btor.Const(1)
    for solution in solutions:
        condition = btor.Const(0)
        for i in range(len(masks)):
            if i == 2: # [u,v,p,q][2] = p
                for mi in range(key_words-1):
                    for j in range(len(masks[i][mi])):
                        condition |= (masks[i][mi][j] != btor.Const(solution[j][i][mi], word_size))
                continue
            for j in range(len(masks[i])):  # round
                condition |= (masks[i][j] != btor.Const(solution[j][i], word_size))
        distinctness_condition &= condition
    btor.Assume(distinctness_condition)

def solve_all(btor, weight, w, nb_rounds, word_size, key_words):
    # key
    p = [[btor.Match_by_symbol("p%d_%d" % (mi, i)) for i in range(nb_rounds + 1)] for mi in range(key_words - 1)]
    q = [btor.Match_by_symbol("q%d" % i) for i in range(nb_rounds + 1)]

    # Get variables
    u = [btor.Match_by_symbol("u%d" % i) for i in range(nb_rounds + 1)]
    v = [btor.Match_by_symbol("v%d" % i) for i in range(nb_rounds + 1)]

    solutions = []
    while True:
        btor.Assume(weight == w)
        distinctness_constraint(btor, [u, v, p, q], solutions, word_size, key_words)

        r = btor.Sat()
        if r != btor.SAT:
            return solutions
        

        solutions.append([
            (
                int(u[i].assignment, base=2),
                int(v[i].assignment, base=2),
                tuple(int(p[mi][i].assignment, base=2) for mi in range(key_words - 1)),
                int(q[i].assignment, base=2)
            )
            for i in range(nb_rounds + 1)
        ])

def parity(x):
    return bin(x).count('1') % 2

def rotl(x, r, word_size):
   mask = (1 << word_size) - 1
   return ((x << r) & mask) | ((x & mask) >> (word_size - r))
    

def compute_sign_speck(differences, masks, word_size, key_words):
    """Compute the correlation sign of a SPECK quasidifferential trail."""

    def pseudoinverseM(t):
        t = t ^ ((t << 1) % 2 ** word_size)
        return t >> 1

    def complement(t):
        return (2 ** word_size - 1) ^ t

    s = 1

    # data path
    for i in range(len(masks) - 1):
        if word_size == 16:
            u = rotl(masks[i][0], word_size - 7, word_size)
            v = rotl(masks[i + 1][1], word_size - 2, word_size) ^ masks[i][1]
        else:
            u = rotl(masks[i][0], word_size - 8, word_size)
            v = rotl(masks[i + 1][1], word_size - 3, word_size) ^ masks[i][1]
        w = masks[i + 1][0] ^ masks[i + 1][1]
        (u_, v_, _) = (u ^ w, v ^ w, u ^ v ^ w)

        (l , b) = differences[i]
        c       = differences[i + 1][0]
        if word_size == 16:
            a = rotl(l, word_size - 7, word_size) 
        else:
            a = rotl(l, word_size - 8, word_size) 
        (a_, b_, c_) = (b ^ c, a ^ c, pseudoinverseM(a ^ b ^ c))
        p1 = parity(((complement(a_) & u_) ^ (c_ & v_)) & ((complement(b_) & v_) ^ (c_ & u_)))
        p2 = parity((u_ & v_) & (c_ ^ (a_ & b_ & complement(c_))))
        s *= (-1) ** (p1 ^ p2)

    # A linear key schedule only propagates masks and contributes no sign.
    return s
