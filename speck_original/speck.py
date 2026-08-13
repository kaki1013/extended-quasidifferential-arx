import common

# speck parameter
block_size = 32
num_rounds = 6
weight = 13  # weight of trail
idx = 0  # index (there would be other trails with same weight)
key_words = 4  # 32 : 4 / 48 : 3, 4 / 64 : 3, 4 / 96 : 2, 3 / 128 : 2, 3, 4 (among 2, 3, 4)

# Validate the supported block-size/key-word combinations.
assert block_size == 32 and key_words == 4 or \
       block_size == 48 and key_words in [3, 4] or \
       block_size == 64 and key_words in [3, 4] or \
       block_size == 96 and key_words in [2, 3] or \
       block_size == 128 and key_words in [2, 3, 4], \
    f"Invalid key_words value for block_size {block_size}."

# SMT parameter
max_weight_loss = 10

word_size = block_size // 2
filename = f"./data/speck_{2*word_size}_{num_rounds}r_{weight}_{idx}.txt"
diffs = common.get_diff(filename)

# quasidifferential trail search
sols = []
btor, weight = common.speck_quasidifferential_trails(diffs, word_size, key_words)
for w in range(max_weight_loss):
    print(f"w = {w} ...")
    for solution in common.solve_all(btor, weight, w, len(diffs) - 1, word_size, key_words):
        sols.append((w, solution))

format_length = word_size // 4
print("Characteristic:")
for (a, b) in diffs:
    print(" " * 4 + f"{a:0{format_length}x} {b:0{format_length}x}")

for (weight, solution) in sols:
    s = common.compute_sign_speck(diffs, solution, word_size, key_words)
    print("Weight: {} [{:+}]".format(weight, s))
    for (u, v, p_tuple, q) in solution:
        p_str = " ".join(f"{pi:0{format_length}x}" for pi in p_tuple)
        print(f"    {p_str} {q:0{format_length}x} {u:0{format_length}x} {v:0{format_length}x}")

# calculating EDP
print("Summary:")
print(f"{'s':>3} {'weight':>8} {'value':>15}")
total_avg_proba = 0
for (weight, solution) in sols:
    s = common.compute_sign_speck(diffs, solution, word_size, key_words)
    print(f"{s:>3} {weight:>8} {s * 2**(-weight):>15.9f}")
    total_avg_proba += s * 2**(-weight)
print("EDP:", total_avg_proba)

# Print the cumulative normalized correlation after each weight.
print()
last_weight = -1
total_avg_proba = 0
for (weight, solution) in sols:
    if last_weight >= 0 and last_weight != weight:
        print(last_weight, total_avg_proba)

    last_weight = weight

    s = common.compute_sign_speck(diffs, solution, word_size, key_words)
    total_avg_proba += s * 2**(-weight)
print(last_weight, total_avg_proba)
