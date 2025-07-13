import pandas as pd
import numpy as np
import random


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

iterations = 0

def load_fpbase_data_as_map(filepath="fpbase.csv") -> dict:
    df = pd.read_csv(filepath)
    return {
        row["sequence"]: {
            "excitation_wavelength": row["excitation_wavelength"],
            "emission_wavelength": row["emission_wavelength"],
            "stokes_shift": row["stokes_shift"],
            "brightness": row["brightness"]
        }
        for _, row in df.iterrows()
    }

def mutate_sequence(seq: str) -> str:
    global iterations
    iterations += 1

    if iterations % 1000 == 0:
        print(f"mutate_sequence called {iterations} times, valid matches = {0}")

    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < 0.2:  # 20% substitution rate
            original = seq[i]
            seq[i] = random.choice([aa for aa in AMINO_ACIDS if aa != original])
    return ''.join(seq)

def generate_mutants(seq: str, num_mutants: int, sequence_map: dict[str, dict]) -> list[str]:
    mutants = set()
    while len(mutants) < num_mutants:
        mutant = mutate_until_valid(seq, sequence_map)
        if mutant and mutant != seq:
            mutants.add(mutant)
    return list(mutants)

def mutate_until_valid(base_seq: str, sequence_map: dict[str, dict]) -> str:
    max_attempts = 1000  # prevent infinite loops
    for _ in range(max_attempts):
        mutant = mutate_sequence(base_seq)
        if mutant in sequence_map:
            return mutant
    return None  # return None if no valid mutant found after max_attempts

def evaluate_gfp_variants(sequences: list[str], sequence_map: dict[str, dict]) -> dict[str, float]:
    brightness_weight = 0.5
    fitnesses: dict[str, float] = {}

    EXCITATION_MAX_THRESHOLD = 510
    EMISSION_MIN_THRESHOLD = 500

    for seq in sequences:
        if seq not in sequence_map:
            continue  # skip if unknown

        props = sequence_map[seq]

        fitness = props["brightness"] * brightness_weight

        if (props["excitation_wavelength"] > EXCITATION_MAX_THRESHOLD or
            props["emission_wavelength"] < EMISSION_MIN_THRESHOLD):
            fitness = 0  # penalty

        fitnesses[seq] = fitness

    return fitnesses

def get_top_keys(data: dict[str, float], x: int) -> list[str]:
    sorted_items = sorted(data.items(), key=lambda item: item[1], reverse=True)
    return [key for key, _ in sorted_items[:x]]

# Example usage
if __name__ == "__main__":
    sequence_map = load_fpbase_data_as_map()
    known_sequences = list(sequence_map.keys())

    num_mutants = 500
    num_generations = 5
    cutoff = 5
    new_mutant_sequences = []

    # First round: pick 5 random known sequences as parents
    base_parents = random.sample(known_sequences, 5)
    for base_sequence in base_parents:
        new_mutant_sequences += generate_mutants(base_sequence, num_mutants, sequence_map)

    # Filter to only valid mutants in the dataset
    valid_mutants = [s for s in new_mutant_sequences if s in sequence_map]

    print(f"Mutants generated: {len(new_mutant_sequences)}")
    print(f"Valid in map: {len(valid_mutants)}")

    # Evaluate and filter mutants based on local map
    fitness_scores = evaluate_gfp_variants(valid_mutants, sequence_map)
    top_keys = get_top_keys(fitness_scores, cutoff)

    # Evolution loop
    for _ in range(num_generations):
        new_mutant_sequences.clear()

        for sequence in top_keys:
            if sequence in sequence_map:
                new_mutant_sequences += generate_mutants(sequence, num_mutants, sequence_map)

        valid_mutants = [s for s in new_mutant_sequences if s in sequence_map]
        print(f"Generation produced {len(new_mutant_sequences)} mutants, {len(valid_mutants)} valid.")

        fitness_scores = evaluate_gfp_variants(valid_mutants, sequence_map)
        top_keys = get_top_keys(fitness_scores, cutoff)

    print(f"\nTop {len(top_keys)} Sequences:")
    for i, seq in enumerate(top_keys):
        print(f"{i+1}: {seq[:30]}...")  # Print first 30 characters of sequence