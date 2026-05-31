import torch
import torch.nn.functional as F
import esm
import numpy as np

# Load ESM-2
model_name = "esm2_t6_8M_UR50D"
model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
batch_converter = alphabet.get_batch_converter()
model.eval()

def get_embedding(seq):
    """Compute mean-pooled embedding for a sequence"""
    data = [("protein", seq)]
    _, _, batch_tokens = batch_converter(data)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[model.num_layers], return_contacts=False)
    token_representations = results["representations"][model.num_layers]
    embedding = token_representations[:, 1:-1].mean(1)
    return embedding

def evaluate_state(suffix_list):
    """
    Compute a score for a state (list of suffix sequences)
    using mean pairwise cosine similarity of embeddings
    """
    embeddings = [get_embedding(s) for s in suffix_list]
    embeddings = torch.cat(embeddings, dim=0)
    n = len(embeddings)
    if n <= 1:
        return 1.0  # single sequence, max similarity

    sim_sum = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            sim = F.cosine_similarity(
                embeddings[i].unsqueeze(0),
                embeddings[j].unsqueeze(0)
            ).item()
            sim_sum += sim
            count += 1
    return sim_sum / count  # mean pairwise similarity

# Example usage in beam search
beam_states = [
    ["PRTEINOSEQWENCE", "PRTEINOSEQYENCE", "AAAAABBB"],
    ["PRTEINOSEQWENCE", "MKVLYNLQSEK"]
]

for state in beam_states:
    score = evaluate_state(state)
    print(f"State {state} score: {score:.4f}")
    
    
