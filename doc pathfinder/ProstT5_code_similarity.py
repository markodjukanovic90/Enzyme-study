import torch
import torch.nn.functional as F
import esm
import numpy as np

# Example: a set of protein sequences
sequences = [
    "PRTEINOSEQWENCE",
    "PRTEINOSEQYENCE",
    "MKVLYNLQSEK",
    "GATCGATCGA"
]

# Load small ESM-2 model
model_name = "esm2_t6_8M_UR50D"
print("Loading ESM-2 model...")
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
    embedding = token_representations[:, 1:-1].mean(1)  # mean over residues
    return embedding

# Compute embeddings for all sequences
embeddings = []
for seq in sequences:
    emb = get_embedding(seq)
    embeddings.append(emb)

# Stack embeddings
embeddings = torch.cat(embeddings, dim=0)

# Compute cosine similarity matrix
num_seq = len(sequences)
similarity_matrix = np.zeros((num_seq, num_seq))
for i in range(num_seq):
    for j in range(num_seq):
        similarity_matrix[i, j] = F.cosine_similarity(
            embeddings[i].unsqueeze(0), embeddings[j].unsqueeze(0)
        ).item()

# Print similarity matrix
print("\nSequence similarity matrix (cosine similarity):")
print(similarity_matrix)
