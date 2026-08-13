
def is_power_of_two(x):
    return x & (x - 1) == 0

def next_power_of_two(x):
    """ Suppose x < 2^{32} """
    x -= 1
    x |= x >> 1
    x |= x >> 2
    x |= x >> 4
    x |= x >> 8
    x |= x >> 16
    x += 1
    return x

def shift_left(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return x << offset
    else:
        n = next_power_of_two(nb_bits)
        y = btor.Concat(btor.Const(0, n - nb_bits), x) << offset
        return y[nb_bits - 1:]

def shift_right(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return x >> offset
    else:
        n = next_power_of_two(nb_bits)
        y = btor.Concat(btor.Const(0, n - nb_bits), x) >> offset
        return y[nb_bits - 1:]

def rotate_right(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return btor.Ror(x, offset)
    else:
        return shift_left(btor, x, nb_bits - offset, nb_bits) | \
               shift_right(btor, x, offset, nb_bits)

def rotate_left(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return btor.Rol(x, offset)
    else:
        return shift_left(btor, x, offset, nb_bits) | \
               shift_right(btor, x, nb_bits - offset, nb_bits)

def M_pseudoinverse(btor, t, nb_bits):
    t = t ^ shift_left(btor, t, 1, nb_bits)
    return shift_right(btor, t, 1, nb_bits)

def M_transpose(btor, t, nb_bits):
    t = shift_right(btor, t, 1, nb_bits)
    i = 1
    while i < nb_bits:
        t = t ^ shift_right(btor, t, i, nb_bits)
        i *= 2
    return t

def hamming_weight_16(btor, x, nb_bits):
    x -= shift_right(btor, x, 1, nb_bits) & 0x5555
    x = (x & 0x3333) + (shift_right(btor, x, 2, nb_bits) & 0x3333)
    x = (x + shift_right(btor, x, 4, nb_bits)) & 0x0F0F
    x += shift_right(btor, x, 8, nb_bits)
    return x & 0x003F

def hamming_weight_32(btor, x, nb_bits):
    x -= shift_right(btor, x, 1, nb_bits) & 0x55555555
    x = (x & 0x33333333) + (shift_right(btor, x, 2, nb_bits) & 0x33333333)
    x = (x + shift_right(btor, x, 4, nb_bits)) & 0x0F0F0F0F
    x += shift_right(btor, x, 8, nb_bits)
    x += shift_right(btor, x, 16, nb_bits)
    return x & 0x0000003F

def hamming_weight(btor, x, nb_bits):
    if nb_bits == 16:
        return hamming_weight_16(btor, x, nb_bits)
    elif nb_bits == 24:
        return hamming_weight(btor, btor.Concat(btor.Const(0, 8), x), 32)[23:]
    elif nb_bits == 32:
        return hamming_weight_32(btor, x, nb_bits)
    elif nb_bits == 48:
        return btor.Concat(btor.Const(0, 16), hamming_weight(btor, x[31:], 32))\
             + btor.Concat(btor.Const(0, 32), hamming_weight(btor, x[:32], 16))
    elif nb_bits == 64:
        return btor.Concat(btor.Const(0, 32), hamming_weight(btor, x[31:], 32))\
             + btor.Concat(btor.Const(0, 32), hamming_weight(btor, x[:32], 32))

def modular_addition(btor, a, b, c, u, v, w, nb_bits, i):
    """Encode Theorem 5's nonzero-QDTM constraints and correlation weight.

    The return value is the negative-log2 magnitude contribution of this
    modular-addition transition.  Its sign is evaluated separately.
    """
    a_ = b ^ c
    b_ = a ^ c
    c_ = M_pseudoinverse(btor, a ^ b ^ c, nb_bits)

    u_ = u ^ w
    v_ = v ^ w
    w_ = M_transpose(btor, u ^ v ^ w, nb_bits)

    n = nb_bits - 1

    btor.Assert((u_ | v_) & ~(a_ | b_ | w_) == 0)
    btor.Assert((a_ & u_) ^ (b_ & v_) == (c_ & w_))
    btor.Assert(((a_[n] == 0) & (b_[n] == 0)) | (a_[n] & u_[n] == u_[n] ^ v_[n]))

    weight = hamming_weight(btor, w_ & ~a_ & ~b_, nb_bits) 
    extra = btor.Cond(
        (a_[n] | b_[n]) & ((u_[n] ^ v_[n]) == (a_[n] & u_[n])) & (u_[n] != v_[n]),
        btor.Const(1, nb_bits),
        btor.Const(0, nb_bits)
    )
    return weight - extra
