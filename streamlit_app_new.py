import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import google.generativeai as genai

# -----------------------------
# PAGE CONFIG
# -----------------------------

st.set_page_config(
    page_title="RNA-seq Results Interpreter",
    layout="wide"
)

st.title("RNA-seq Results Interpreter")

st.markdown("""
Upload DESeq2 differential expression results and receive:
- QC summaries
- Volcano plots
- AI-generated biological interpretation
- Suggested next analyses
""")

# -----------------------------
# CONFIGURE GEMINI
# -----------------------------

api_key = st.secrets.get("GEMINI_API_KEY")

if not api_key:
    st.error("Gemini API key not found.")
    st.stop()

genai.configure(api_key=api_key)

model = genai.GenerativeModel("gemini-3.1-pro-preview")

# -----------------------------
# FILE UPLOAD
# -----------------------------

uploaded_file = st.file_uploader(
    "Upload DESeq2 Results CSV",
    type=["csv"]
)

if uploaded_file:

    df = pd.read_csv(uploaded_file)

    st.subheader("Preview")

    st.dataframe(df.head())

    required_cols = [
        "log2FoldChange",
        "padj"
    ]

    missing_cols = [
        col for col in required_cols
        if col not in df.columns
    ]

    if missing_cols:
        st.error(f"Missing required columns: {missing_cols}")
        st.stop()

    # -----------------------------
    # BASIC QC
    # -----------------------------

    sig = df[df["padj"] < 0.05]

    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]

    col1, col2, col3 = st.columns(3)

    col1.metric("Significant Genes", len(sig))
    col2.metric("Upregulated", len(up))
    col3.metric("Downregulated", len(down))

    # -----------------------------
    # VOLCANO PLOT
    # -----------------------------

    st.subheader("Volcano Plot")

    df["neglog10_padj"] = -np.log10(
        df["padj"].replace(0, 1e-300)
    )

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        df["log2FoldChange"],
        df["neglog10_padj"],
        s=8,
        alpha=0.6
    )

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10 adjusted p-value")

    st.pyplot(fig)

    # -----------------------------
    # TOP GENES
    # -----------------------------

    top_up = up.sort_values(
        "log2FoldChange",
        ascending=False
    ).head(10)

    top_down = down.sort_values(
        "log2FoldChange"
    ).head(10)

    with st.expander("Top Upregulated Genes"):
        st.dataframe(top_up)

    with st.expander("Top Downregulated Genes"):
        st.dataframe(top_down)

    # -----------------------------
    # AI INTERPRETATION
    # -----------------------------

    if st.button("Generate AI Interpretation"):

        prompt = f"""
        You are an expert bioinformatics assistant.

        Analyze the following RNA-seq differential expression summary.

        Significant genes: {len(sig)}
        Upregulated genes: {len(up)}
        Downregulated genes: {len(down)}

        Top upregulated genes:
        {top_up.to_string()}

        Top downregulated genes:
        {top_down.to_string()}

        Provide:
        1. Biological interpretation
        2. Potential pathways involved
        3. Potential technical concerns
        4. Suggested next analyses
        5. Concise manuscript-style summary
        """

        with st.spinner("Generating interpretation..."):

            try:

                response = model.generate_content(prompt)

                st.subheader("AI Interpretation")

                st.markdown(response.text)

            except Exception as e:
                st.error(f"Error generating interpretation: {e}")