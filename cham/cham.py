import common
import time


def out(msg):
    print(msg)
    if save_file:
        f.write(msg + "\n")

        
# CHAM parameter (1): uncomment this block to analyze CHAM-64/128.
# block_size = 64
# num_rounds = 39
# weight = 64  # weight of the differential characteristic
# idx = 0  # distinguishes characteristics having the same weight
# key_size = 128

# CHAM parameter (2): default configuration (CHAM-128/128).
block_size = 128
num_rounds = 63
weight = 127  # weight of the differential characteristic
idx = 0  # distinguishes characteristics having the same weight
key_size = 128

# CHAM parameter (3): uncomment this block to analyze CHAM-128/256.
# block_size = 128
# num_rounds = 63
# weight = 127  # weight of trail
# idx = 0  # index (there would be other trails with same weight)
# key_size = 256

word_size = block_size // 4
key_words = key_size // word_size

# Validate the supported block-size/key-size combinations.
assert (block_size == 64 and key_size == 128) or \
       (block_size == 128 and key_size in [128, 256])

# Enumerate trails whose relative-correlation weight loss satisfies
# start_weight_loss <= w < max_weight_loss. A trail contributes sign * 2^-w
# after normalization by the all-zero-mask trail.
start_weight_loss = 0
max_weight_loss = 8

save_file = True

filename = f"./data/cham_{4*word_size}_{num_rounds}r_{weight}_{idx}.txt"
diffs = common.get_diff(filename)

output_file = f"./result_cham_{4*word_size}_{num_rounds}r_{weight}_{idx}_start_weight_loss_{start_weight_loss}_max_weight_{max_weight_loss}.txt"
if save_file:
    f = open(output_file, "w")

# quasidifferential trail search
format_length = word_size // 4

sols = []
btor, trail_weight  = common.cham_quasidifferential_trails(diffs, word_size, key_words)


# out("Characteristic:")
# for (a,b,c,d) in diffs:
#     out(f"    {a:0{format_length}x} {b:0{format_length}x} {c:0{format_length}x} {d:0{format_length}x}")


total_avg_proba = 0

# out("")
# out("Cumulative EDP:")
# out(f"{'trail_weight ':>8} {'EDP':>12} {'time':>10}")

for w in range(start_weight_loss, max_weight_loss):
    out("=" * 60)
    out(f"Weight loss = {w}")
    out("=" * 60)

    start = time.perf_counter()
    sols = []
    for solution in common.solve_all(btor, trail_weight , w, len(diffs)-1, word_size, key_words):
        s = common.compute_sign_cham(diffs, solution, word_size, key_words)
        sols.append((solution, s))
        total_avg_proba += s * 2**(-w)

    elapsed = time.perf_counter() - start

    # out("")
    # out("Trails:")
    # for solution, s in sols:
    #     out(f"Weight: {w} [{s:+}]")
    #     for ua,ub,uc,ud in solution:
    #         out(f"    {ua:0{format_length}x} {ub:0{format_length}x} {uc:0{format_length}x} {ud:0{format_length}x}")

    out("")
    out(f"Summary ({len(sols)} trails):")
    out(f"{'s':>3} {'weight':>8} {'value':>10}")

    for _, s in sols:
        out(f"{s:>3} {w:>8} {s*2**(-w):>10.6f}")

    out("")
    out("Cumulative EDP:")
    out(f"{'weight':>8} {'EDP':>12} {'weight time':>14}")
    out(f"{w:>8} {total_avg_proba:>12.6f} {elapsed:>13.2f}s")
    out("")

if save_file:
    f.close()
