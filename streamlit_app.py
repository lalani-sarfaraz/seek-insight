import streamlit as st
import pandas as pd


st.title("RNA-seq Results Interpreter")

uploaded_file = st.file_uploader(
    "Upload DESeq2 Results CSV",
    type=["csv"]
)

if uploaded_file:
    df = pd.read_csv(uploaded_file)

    st.write(df.head())


# Basic QC checks:

sig = df[df["padj"] < 0.05]

up = sig[sig["log2FoldChange"] > 0]
down = sig[sig["log2FoldChange"] < 0]

st.write(f"Total genes: {len(df)}")
st.write(f"Significant genes: {len(sig)}")
st.write(f"Upregulated: {len(up)}")
st.write(f"Downregulated: {len(down)}")




# Volcano plots:

import matplotlib.pyplot as plt
import numpy as np

df["neglog10_padj"] = -np.log10(df["padj"])

fig, ax = plt.subplots()

ax.scatter(
    df["log2FoldChange"],
    df["neglog10_padj"],
    s=5
)

ax.set_xlabel("log2 fold change")
ax.set_ylabel("-log10 adjusted p-value")

ax.set_title("Volcano plot")


st.pyplot(fig)




# LLM-based interpretation:

top_up = up.sort_values(
    "log2FoldChange",
    ascending=False
).head(10)

top_down = down.sort_values(
    "log2FoldChange"
).head(10)


st.write(f"""
         Top upregulated genes:""",
         top_up.head(10))

st.write(f"""         
         Top downregulated genes:""",
         top_down.head(10))




# Connect to OpenAI API:

from openai import OpenAI

api_key = st.secrets.get("OPENAI_API_KEY")

if not api_key:
    st.error("OPENAI_API_KEY not found.")
    st.stop()

client = OpenAI(api_key=api_key)

response = client.chat.completions.create(
    model="gpt-4.1-mini",
    messages=[
        {"role": "user", "content": prompt}
    ]
)

analysis = response.choices[0].message.content

st.write(analysis)