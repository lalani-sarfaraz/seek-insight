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

st.pyplot(fig)