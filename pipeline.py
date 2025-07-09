import pandas as pd
import numpy as np
from Bio.SeqUtils import CodonAdaptationIndex
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
from Bio.Blast import NCBIWWW, NCBIXML

#prints out the 5 best gfp sequences 
if __name__ == "__main__":

    #base sequence is for Aequorea victoria (first amino acid)
    baseSequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

    newMutantSequences = []
    numMutants = 100 #number of mutations for each sequence
    numGenerations = 50 #number of generations
    cutoff = 5 #amount of surviving mutant sequences
    for i in range(5):
        newMutantSequences = append(generate_mutants(baseSequence, numMutants)

    #gets actual protein name from sequence (if known)
    for sequence in newMutantSequences {
        # Run BLAST (this can take time and requires internet)
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
        blast_record = NCBIXML.read(result_handle)
        top_hit = blast_record.alignments[0]
        Mutant_gfp_names.append(top_hit.hit_def)

    #gets the corresponding fittness for each mutant
    df_results = evaluate_gfp_variants(Mutant_gfp_names)

    #array of the mutant names with the top 5 highest fittnesses
    top_keys = get_top_keys(df_results, cutoff)
    
    #repeats process with the keys for numGenerations amount of times
    for i in range(numGenerations):
        
        for sequence := topKeys:
            newMutantSequences.append(generate_mutants(sequence, numMutants)
    
        for sequence in newMutantSequences:
            # Run BLAST (this can take time and requires internet)
            result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
            blast_record = NCBIXML.read(result_handle)
            top_hit = blast_record.alignments[0]
            Mutant_gfp_names.append(top_hit.hit_def)
    
        df_results = evaluate_gfp_variants(Mutant_gfp_names)
        top_keys = get_top_keys(df_results, cutoff)
        Mutant_gfp_names.empty()
        newMutantSequences.empty()

    #prints out best resulting mutant names
    for i in range(cutoff):
        println(i + ": " top_keys[i])
        


#takes in a fp and returns the properties from fpbase as a map
def get_fp_properties_from_name(fp_name):
    import requests

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

    response = requests.post(
        'https://www.fpbase.org/graphql/',
        json={'query': query}
    )

    return response.json()['data']['fluorescentProtein']


#input: amino acid sequences
#output: map of coresponding fittnesses to those sequences
def evaluate_gfp_variants(fpVarients):
    brightnessWeight = 0.5

    fittnesses = dict[string, float] #key stirng, val float
    for fpVarient in fpVarients:
        fpProperties = get_fp_properties_from_name(fpVarient)
        fittness = fpProperties['brightness'] * brightnessWeight
        if fp['excitationMax'] > x || fp['emissionMax'] < y || fp['aggregation'] != z || fp['maturation'] != ab: #penalties/constraints
            fittness -= 1000.0
        fittnesses[fpVarient] = fittness
    


    return fittnesses

"""
Args:
    data: Dictionary with string keys and integer values.
    x: Number of top entries to return.

Returns:
    List of keys with the top x highest values.
"""
def get_top_keys(data: dict[str, int], x: int) -> list[str]:
    # Sort items by value in descending order
    sorted_items = sorted(data.items(), key=lambda item: item[1], reverse=True)

    # Extract just the keys of the top x items
    top_keys = [key for key, _ in sorted_items[:x]]
    
    return top_keys
import random

"""
Applies a realistic mix of mutations (substitution, insertion, deletion)
to a protein sequence.
"""
def mutate_sequence(seq: str) -> str:
    # Standard amino acids
    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
    seq = list(seq)
    i = 0
    while i < len(seq):
        mutation_chance = random.random()
        if mutation_chance < 0.05:
            # Deletion
            del seq[i]
            continue  # skip incrementing i to avoid skipping next position
        elif mutation_chance < 0.10:
            # Insertion
            seq.insert(i, random.choice(AMINO_ACIDS))
            i += 1  # move past inserted amino acid
        elif mutation_chance < 0.30:
            # Substitution
            original = seq[i]
            new = random.choice([aa for aa in AMINO_ACIDS if aa != original])
            seq[i] = new
            i += 1
        else:
            # No mutation
            i += 1

    return ''.join(seq)


"""
Generates a list of mutated amino acid sequences.

Args:
    seq: Original amino acid sequence.
    num_mutants: Number of mutants to generate.

Returns:
    List of mutated sequences (strings).
"""
def generate_mutants(seq: str, num_mutants: int) -> list[str]:
    mutants = set()
    while len(mutants) < num_mutants:
        mutant = mutate_sequence(seq)
        # Ensure the mutant is actually different
        if mutant != seq:
            mutants.add(mutant)
    return list(mutants)
