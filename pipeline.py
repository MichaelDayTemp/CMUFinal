import pandas as pd
import numpy as np
from Bio.SeqUtils import CodonAdaptationIndex
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
from Bio.Blast import NCBIWWW, NCBIXML
import requests
import random

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# -=-=-=-=-=Core Functions=-=-=-=-=-=-

# Generates 1 random mutation (sub, indel) in sequence
def mutate_sequence(seq: str) -> str:
    seq = list(seq)
    i = 0
    while i < len(seq):
        mutation_chance = random.random()
        if mutation_chance < 0.05:
            del seq[i]
            continue
        elif mutation_chance < 0.10:
            seq.insert(i, random.choice(AMINO_ACIDS))
            i += 1
        elif mutation_chance < 0.30:
            original = seq[i]
            new = random.choice([aa for aa in AMINO_ACIDS if aa != original])
            seq[i] = new
            i += 1
        else:
            i += 1
    return ''.join(seq)

#Generates array of size num_mutants of mutated versions of the sequence
def generate_mutants(seq: str, num_mutants: int) -> list[str]:
    mutants = set()
    while len(mutants) < num_mutants:
        mutant = mutate_sequence(seq)
        if mutant != seq:
            mutants.add(mutant)
    return list(mutants)

#returns a list of the top x best fitting sequences
def get_top_keys(data: dict[str, float], x: int) -> list[str]:
    sorted_items = sorted(data.items(), key=lambda item: item[1], reverse=True)
    return [key for key, _ in sorted_items[:x]]

#fetches fp properties from dpbase and returns them as a map
def get_fp_properties_from_name(fp_name):
    query = """
    {{
      fluorescentProtein(name: "{fp_name}") {{
        name
        brightness
        quantumYield
        extinctionCoefficient
        excitationMax
        emissionMax
      }}
    }}
    """
    response = requests.post('https://www.fpbase.org/graphql/', json={'query': query})
    return response.json()['data']['fluorescentProtein']

#returns map of fittnesses of a given amount of fp mutants 
def evaluate_gfp_variants(fp_variants):
    brightness_weight = 0.5
    fitnesses: dict[str, float] = {}

    # Constraints (example values — replace with real thresholds if you have them)
    EXCITATION_MAX_THRESHOLD = 510
    EMISSION_MIN_THRESHOLD = 500

    for fp in fp_variants:
        fp_properties = get_fp_properties_from_name(fp)
        if not fp_properties:
            continue  # skip if data is missing

        fitness = fp_properties['brightness'] * brightness_weight

        # Apply penalties
        if (fp_properties['excitationMax'] > EXCITATION_MAX_THRESHOLD or
            fp_properties['emissionMax'] < EMISSION_MIN_THRESHOLD):
            fitness -= 1000.0

        fitnesses[fp] = fitness

    return fitnesses

# ================== Main Evolution Loop ==================

if __name__ == "__main__":
    base_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

    num_mutants = 100
    num_generations = 50
    cutoff = 5  # top survivors
    mutant_gfp_names = []
    new_mutant_sequences = []

    # First round of mutation
    for _ in range(5):
        new_mutant_sequences += generate_mutants(base_sequence, num_mutants)

    for sequence in new_mutant_sequences:
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
        blast_record = NCBIXML.read(result_handle)
        top_hit = blast_record.alignments[0]
        mutant_gfp_names.append(top_hit.hit_def)

    df_results = evaluate_gfp_variants(mutant_gfp_names)
    top_keys = get_top_keys(df_results, cutoff)

    # Evolution loop
    for _ in range(num_generations):
        new_mutant_sequences.clear()
        mutant_gfp_names.clear()

        for sequence in top_keys:
            new_mutant_sequences += generate_mutants(sequence, num_mutants)

        for sequence in new_mutant_sequences:
            result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
            blast_record = NCBIXML.read(result_handle)
            top_hit = blast_record.alignments[0]
            mutant_gfp_names.append(top_hit.hit_def)

        df_results = evaluate_gfp_variants(mutant_gfp_names)
        top_keys = get_top_keys(df_results, cutoff)

    # Final results
    for i in range(cutoff):
        print(f"{i + 1}: {top_keys[i]}")
