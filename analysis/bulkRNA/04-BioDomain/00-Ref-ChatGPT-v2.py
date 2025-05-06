from openai import OpenAI
import openai

# read the API key from the OpenAI_APi.txt file
with open("../../../data/ref/OpenAI_API.txt", "r") as f:
    api_key = f.read()
    
client = OpenAI(
  api_key=api_key
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



# # 1. Create a new env named “openai” with Python 3.9 (or 3.10 if you prefer)
# conda create -n openai python=3.9 -y

# # 2. Activate it
# conda activate openai

# # 3. Install the OpenAI client (and any extras you like)
# pip install openai tqdm pandas

# # 4. (Optional) Verify
# python - <<'PYCODE'
# import openai, tqdm, pandas
# print("OpenAI lib version:", openai.__version__)
# print("tqdm version:", tqdm.__version__)
# print("pandas version:", pandas.__version__)
# PYCODE

import re
import time
import pandas as pd
import openai

from openai import OpenAI
from tqdm import tqdm
# from statsmodels.stats.multitest import multipletests

# (1) if you haven't already, instantiate a client
#    – this is optional if you set OPENAI_API_KEY in your env
# client = openai.OpenAI()

# --- 1) Load and clean AD reference, dropping any pathway that maps to multiple Biodomains
ref = pd.read_csv("./BioDomain_AD_Ref.csv")[['pathway', 'Biodomain']]
ref_cleaned = (
    ref
    .dropna(subset=['pathway', 'Biodomain'])
    .query("Biodomain.str.lower() != 'none'", engine='python')
)
dup_paths = ref_cleaned['pathway'][ref_cleaned['pathway'].duplicated(keep=False)].unique()
ref_cleaned = ref_cleaned[~ref_cleaned['pathway'].isin(dup_paths)]

# --- 2) Formatter to normalize pathway names
def format_pathway(p):
    p = re.sub(r'^GO[_]?BP[_]?', '', p)
    return p.replace('_', ' ').lower().strip()

ref_cleaned['pathway'] = ref_cleaned['pathway'].map(format_pathway)
ref_context = "\n".join(
    f"{row['pathway']} -> {row['Biodomain']}"
    for _, row in ref_cleaned.iterrows()
)

# --- 3) Single‐pathway test
new_pathway = "GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE"
new_pathway_fmt = format_pathway(new_pathway)

prompt = f"""
You are a biomedical ontology expert. Below are known pathway→Biodomain mappings:

{ref_context}

Based on these examples, assign the most appropriate Biodomain to the following pathway:
Pathway: {new_pathway_fmt}

Please output **only** the most appropriate Biodomain or 'unknown'.
"""

# >>> Here is the *new* v1 call:
response = client.chat.completions.create(
    model="gpt-4o-mini",
    messages=[
        {"role": "system", "content": "You are a biomedical ontology expert. Infer only the most appropriate Biodomain or 'unknown'."},
        {"role": "user",   "content": prompt}
    ],
    temperature=0
)
print("Prediction:", response.choices[0].message.content.strip())























# --- FDR‐correct your mixed‐model results
target = pd.read_csv("./mixed_model_results_all.csv")
_, padj, _, _ = multipletests(target['p_FXS'], method='fdr_bh')
target['padj_FXS'] = padj
sig_targets = target[target['padj_FXS'] < 0.05]
pathway_list = sig_targets['pathway'].unique()[:10]



# --- Loop and assign Biodomains
results = []
for pw in tqdm(pathway_list, desc="Assigning Biodomains"):
    pw_fmt = format_pathway(pw)
    prompt = f"""
You are a biomedical ontology expert. Below are known pathway→Biodomain mappings:

{ref_context}

Based on these examples, assign the most appropriate Biodomain to the following pathway:
Pathway: {pw_fmt}

Please output **only** the most appropriate Biodomain or 'unknown'.
"""
    try:
        resp = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": "You are a biomedical ontology expert. Infer only the most appropriate Biodomain or 'unknown'."},
                {"role": "user",   "content": prompt}
            ],
            temperature=0
        )
        bd = resp.choices[0].message.content.strip()
    except Exception as e:
        bd = f"ERROR: {e}"
    results.append({"pathway": pw, "Biodomain": bd})
    time.sleep(5)

biodomain_df = pd.DataFrame(results)
biodomain_df.to_csv("biodomain_results_0429.csv", index=False)
out = sig_targets.merge(biodomain_df, on="pathway", how="left")
out.to_csv("target_with_biodomain_0429.csv", index=False)
