import os
import re
import common
import time

# Characteristic directory relative to this script's working directory.
folder_path = './data'

def load_data(folder_path):
    # Return only SPECK characteristic files from the data directory.
    if os.path.exists(folder_path):
        files = [file for file in os.listdir(folder_path) if file.startswith('speck')]
        return files
    else:
        print("The 'data' directory does not exist.")
        return []

def get_data_info(filename):
    pattern = r"speck_(\d+)_(\d+)r_(\d+)_(\d+)\.txt"

    match = re.match(pattern, filename)
    if match:
        data_info = [int(x) for x in match.groups()]
        return data_info
    return [-1, -1, -1, -1]


def get_key_words_len(block_size):
    if block_size == 32:
        return [4]
    if block_size == 48:
        return [3, 4]
    if block_size == 64:
        return [3, 4]
    if block_size == 96:
        return [2, 3]
    if block_size == 128:
        return [2, 3, 4]


# SMT parameter
max_weight_loss = 10

save_file = True

if save_file:
    # Open the result.txt file in write mode (or append mode if you want to append results to the file)
    with open(f'result_max_weight_{max_weight_loss}.txt', 'w') as result_file:
        # Write the header to the file
        result_file.write(f"{'EDP':<10} {'block':<8} {'round':<8} {'weight':<8} {'idx':<8} {'key':<8} {'time(s)':<10}\n")
        
        # Loop over each file in the folder
        for filename in load_data(folder_path):
            block_size, num_rounds, trail_weight, idx = get_data_info(filename)

            word_size = block_size // 2
            diffs = common.get_diff(folder_path + '/' + filename)

            # Loop over each key length
            for key_words in get_key_words_len(block_size):
                start_time = time.perf_counter()

                # Quasidifferential trail search
                sols = []
                btor, weight = common.speck_quasidifferential_trails(diffs, word_size, key_words)
                for w in range(max_weight_loss):
                    for solution in common.solve_all(btor, weight, w, len(diffs) - 1, word_size, key_words):
                        sols.append((w, solution))

                total_avg_proba = 0
                for (weight, solution) in sols:
                    s = common.compute_sign_speck(diffs, solution, word_size, key_words)
                    total_avg_proba += s * 2**(-weight)

                elapsed = time.perf_counter() - start_time
                
                # Write the results to the file
                result_file.write(f"{total_avg_proba:<10.4f} {block_size:<8} {num_rounds:<8} {trail_weight:<8} {idx:<8} {key_words:<8} {elapsed:<10.3f}\n")
                result_file.flush()

else:
    print(f"{'EDP':<10} {'block':<8} {'round':<8} {'weight':<8} {'idx':<8} {'key':<8}")
    for filename in load_data(folder_path):
        block_size, num_rounds, trail_weight, idx = get_data_info(filename)

        word_size = block_size // 2
        diffs = common.get_diff(folder_path + '/' + filename)

        for key_words in get_key_words_len(block_size):
            start_time = time.perf_counter()

            # quasidifferential trail search
            sols = []
            btor, weight = common.speck_quasidifferential_trails(diffs, word_size, key_words)
            for w in range(max_weight_loss):
                for solution in common.solve_all(btor, weight, w, len(diffs) - 1, word_size, key_words):
                    sols.append((w, solution))

            total_avg_proba = 0
            for (weight, solution) in sols:
                s = common.compute_sign_speck(diffs, solution, word_size, key_words)
                total_avg_proba += s * 2**(-weight)

            elapsed = time.perf_counter() - start_time
            print(f"{total_avg_proba:<10.4f} {block_size:<8} {num_rounds:<8} {trail_weight:<8} {idx:<8} {key_words:<8} {elapsed:.3f}s")
