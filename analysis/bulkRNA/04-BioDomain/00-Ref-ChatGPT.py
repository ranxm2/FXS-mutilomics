
from openai import OpenAI
import openai

client = OpenAI(
  api_key="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"
)

import pandas as pd
import openai
import time
from tqdm import tqdm

response = client.responses.create(
  model="gpt-4o-mini",
  input="Tell me a three sentence bedtime story about a unicorn."
)

# only show the response content
print(response.output[0].content[0].text)

# Load reference and target pathway data
ref = pd.read_csv("../../../data/ref/BioDomain_AD_Ref.csv")

# Only keep the pathway and Biodomain columns
ref = ref[['pathway', 'Biodomain']]

# count how many characters in the pathway and Biodomain columns
print(ref['pathway'].dropna().apply(lambda x: len(x.split())).sum())
print(ref['Biodomain'].dropna().apply(lambda x: len(x.split())).sum())

ref_cleaned = ref.dropna(subset=['pathway', 'Biodomain'])
ref_cleaned = ref_cleaned[ref_cleaned['Biodomain'].str.lower() != 'none']

# Create the context string
ref_context = "\n".join([
    f"{row['pathway']} -> {row['Biodomain']}"
    for _, row in ref_cleaned.iterrows()  # adjust number for more/less context
])

# Display result
print(ref_context)

# count how many characters in the context string
print(len(ref_context.split()))

new_pathway = "GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE"


# Construct prompt
prompt = f"""
You are a biomedical ontology expert. Below are known examples of pathway to biodomain mappings:

{ref_context}

Based on the above, assign the most appropriate biodomain to the following pathway:
Pathway: {new_pathway}
Biodomain:"""

# Send the request
response = client.chat.completions.create(
    model="gpt-4o-mini",
    messages=[
        {"role": "system", "content": "You are a biomedical ontology expert. Infer the biodomain from pathway names."},
        {"role": "user", "content": prompt}
    ],
    temperature=0
)

# excract the response
response_text = response.choices[0].message['content'].strip()
print(response_text)

# read the target data 

content = response.choices[0].message.content.strip()
biodomain = content.replace("Biodomain:", "").strip()


target = pd.read_csv("../03-Mutil_group_GSVA/results/03-Mix-effect/mixed_model_results_all.csv")

# do FDR correction for p_FXS
target['padj_FXS'] = target['p_FXS']


from statsmodels.stats.multitest import multipletests

# BH correction for p-values in column 'p_FXS'
_, padj, _, _ = multipletests(target['p_FXS'], method='fdr_bh')

# Assign corrected p-values to a new column
target['padj_FXS'] = padj

# only keep the rows with padj_FXS < 0.05
target = target[target['padj_FXS'] < 0.05]

pathway_list = target['pathway'].unique()
len(pathway_list)
pathway_list = pathway_list[:10]
results = []

# Loop through each pathway and query the model
for pathway_index in tqdm(range(len(pathway_list))):
    pathway = pathway_list[pathway_index]
    prompt = f"""
You are a biomedical ontology expert. Below are known examples of pathway to biodomain mappings:

{ref_context}

Based on the above, assign the most appropriate biodomain to the following pathway:
Pathway: {pathway}
Biodomain:"""

    try:
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": "You are a biomedical ontology expert. Infer the biodomain from pathway names."},
                {"role": "user", "content": prompt}
            ],
            temperature=0
        )

        content = response.choices[0].message.content.strip()
        biodomain = content.replace("Biodomain:", "").strip()

    except Exception as e:
        biodomain = f"ERROR: {e}"

    results.append({"pathway": pathway, "Biodomain": biodomain})

    time.sleep(5)  # rate limiting buffer

# Merge back with target
biodomain_df = pd.DataFrame(results)
# save the biodomain_df to csv
biodomain_df.to_csv("biodomain_results.csv", index=False)


target_with_biodomain = target.merge(biodomain_df, on="pathway", how="left")

# Save to CSV
target_with_biodomain.to_csv("target_with_biodomain.csv", index=False)

