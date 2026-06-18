import os
import io
import json
import textwrap
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import gseapy as gp
from groq import Groq


st.set_page_config(
    page_title="DESeq2 Interpreter",
    layout="wide",
    initial_sidebar_state="expanded",
)

DEFAULT_GENESET_CHOICES = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
]

DEFAULT_MODEL = "llama-3.3-70b-versatile"


def find_first_matching_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def standardize_deseq2_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    gene_col = find_first_matching_column(
        df,
        ["gene", "genes", "symbol", "gene_symbol", "external_gene_name", "Gene", "SYMBOL"]
    )
    log2fc_col = find_first_matching_column(
        df,
        ["log2FoldChange", "log2fc", "logFC", "lfc", "Log2FoldChange"]
    )
    padj_col = find_first_matching_column(
        df,
        ["padj", "adj_pval", "adj.P.Val", "FDR", "qvalue", "fdr", "Adjusted P-value"]
    )
    pval_col = find_first_matching_column(
        df,
        ["pvalue", "pval", "P.Value", "PValue"]
    )
    basemean_col = find_first_matching_column(
        df,
        ["baseMean", "mean", "AveExpr"]
    )

    rename_map = {}
    if gene_col:
        rename_map[gene_col] = "gene"
    if log2fc_col:
        rename_map[log2fc_col] = "log2FoldChange"
    if padj_col:
        rename_map[padj_col] = "padj"
    if pval_col:
        rename_map[pval_col] = "pvalue"
    if basemean_col:
        rename_map[basemean_col] = "baseMean"

    df = df.rename(columns=rename_map)

    required = ["gene", "log2FoldChange"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns after standardization: {missing}. "
            f"Need at least a gene column and a log2FoldChange column."
        )

    if "padj" not in df.columns and "pvalue" not in df.columns:
        raise ValueError("Need either a padj column or a pvalue column.")

    if "padj" not in df.columns:
        df["padj"] = df["pvalue"]

    df["gene"] = df["gene"].astype(str).str.strip()
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["gene", "log2FoldChange", "padj"])
    df = df[df["gene"] != ""]
    df = df.drop_duplicates(subset=["gene"], keep="first")

    tiny = 1e-300
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=tiny))

    return df


def classify_genes(
    df: pd.DataFrame,
    padj_thresh: float,
    lfc_thresh: float
) -> pd.DataFrame:
    df = df.copy()
    df["significant"] = (df["padj"] < padj_thresh) & (df["log2FoldChange"].abs() >= lfc_thresh)

    conditions = [
        (df["padj"] < padj_thresh) & (df["log2FoldChange"] >= lfc_thresh),
        (df["padj"] < padj_thresh) & (df["log2FoldChange"] <= -lfc_thresh),
    ]
    choices = ["Up", "Down"]
    df["direction"] = np.select(conditions, choices, default="NS")
    return df


def summarize_results(df: pd.DataFrame) -> Dict:
    up_df = df[df["direction"] == "Up"].copy()
    down_df = df[df["direction"] == "Down"].copy()
    sig_df = df[df["significant"]].copy()

    summary = {
        "n_total": int(df.shape[0]),
        "n_significant": int(sig_df.shape[0]),
        "n_up": int(up_df.shape[0]),
        "n_down": int(down_df.shape[0]),
        "top_up": up_df.sort_values(["padj", "log2FoldChange"], ascending=[True, False]).head(15),
        "top_down": down_df.sort_values(["padj", "log2FoldChange"], ascending=[True, True]).head(15),
    }
    return summary


def make_volcano_plot(df: pd.DataFrame, padj_thresh: float, lfc_thresh: float):
    color_map = {"Up": "#d62728", "Down": "#1f77b4", "NS": "#b0b0b0"}

    fig = px.scatter(
        df,
        x="log2FoldChange",
        y="neg_log10_padj",
        color="direction",
        color_discrete_map=color_map,
        hover_data=["gene", "padj", "log2FoldChange"],
        title="Volcano plot",
        opacity=0.75,
    )

    fig.add_vline(x=lfc_thresh, line_dash="dash", line_color="black")
    fig.add_vline(x=-lfc_thresh, line_dash="dash", line_color="black")
    fig.add_hline(y=-np.log10(max(padj_thresh, 1e-300)), line_dash="dash", line_color="black")

    return fig


@st.cache_data(show_spinner=False)
def run_enrichment(
    genes: List[str],
    libraries: List[str],
    organism: str = "human",
) -> pd.DataFrame:
    genes = [g.strip().upper() for g in genes if isinstance(g, str) and g.strip()]
    genes = list(dict.fromkeys(genes))

    if len(genes) < 5:
        return pd.DataFrame()

    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=libraries,
        organism=organism,
        outdir=None,
    )
    res = enr.results.copy()
    if res.empty:
        return res

    keep_cols = [
        c for c in [
            "Gene_set", "Term", "Overlap", "P-value",
            "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"
        ] if c in res.columns
    ]
    res = res[keep_cols].sort_values("Adjusted P-value", ascending=True)
    return res


def top_terms_for_llm(df: pd.DataFrame, n: int = 10) -> List[Dict]:
    if df is None or df.empty:
        return []

    out = []
    for _, row in df.head(n).iterrows():
        out.append({
            "gene_set": row.get("Gene_set", ""),
            "term": row.get("Term", ""),
            "adjusted_p": float(row.get("Adjusted P-value", np.nan)) if pd.notnull(row.get("Adjusted P-value", np.nan)) else None,
            "overlap": row.get("Overlap", ""),
            "combined_score": float(row.get("Combined Score", np.nan)) if pd.notnull(row.get("Combined Score", np.nan)) else None,
            "genes": row.get("Genes", ""),
        })
    return out


def build_analysis_payload(
    summary: Dict,
    top_up_table: pd.DataFrame,
    top_down_table: pd.DataFrame,
    up_enrich: pd.DataFrame,
    down_enrich: pd.DataFrame,
    params: Dict
) -> Dict:
    return {
        "thresholds": params,
        "summary": {
            "n_total": summary["n_total"],
            "n_significant": summary["n_significant"],
            "n_up": summary["n_up"],
            "n_down": summary["n_down"],
        },
        "top_up_genes": top_up_table[["gene", "log2FoldChange", "padj"]].to_dict(orient="records") if not top_up_table.empty else [],
        "top_down_genes": top_down_table[["gene", "log2FoldChange", "padj"]].to_dict(orient="records") if not top_down_table.empty else [],
        "up_enrichment_top_terms": top_terms_for_llm(up_enrich, n=10),
        "down_enrichment_top_terms": top_terms_for_llm(down_enrich, n=10),
    }


def build_system_prompt() -> str:
    return textwrap.dedent(
        """
        You are an expert bioinformatics assistant interpreting differential expression results.
        You must be careful not to overclaim causality.
        Use only the provided analysis payload.
        If information is insufficient, say so clearly.
        Return concise, scientifically literate answers.

        When asked for the default report, structure your answer with:
        1. Biological interpretation
        2. Potential pathways involved
        3. Potential technical concerns
        4. Suggested next analyses
        5. Concise manuscript-style summary
        """
    ).strip()


def build_default_user_prompt(payload: Dict) -> str:
    return (
        "Interpret the following DESeq2-derived analysis payload.\n\n"
        f"{json.dumps(payload, indent=2)}\n\n"
        "Please answer with the five requested sections."
    )


def call_groq(messages: List[Dict], model: str, temperature: float = 0.2) -> str:
    api_key = os.getenv("GROQ_API_KEY") or st.secrets.get("GROQ_API_KEY", None)
    if not api_key:
        raise ValueError("GROQ_API_KEY not found in environment or Streamlit secrets.")

    client = Groq(api_key=api_key)
    response = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=temperature,
        max_completion_tokens=1600,
    )
    return response.choices[0].message.content


def init_chat_state():
    if "chat_history" not in st.session_state:
        st.session_state.chat_history = []
    if "analysis_payload" not in st.session_state:
        st.session_state.analysis_payload = None
    if "analysis_report" not in st.session_state:
        st.session_state.analysis_report = None


def dataframe_to_csv_download(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


init_chat_state()

st.title("DESeq2 Results Interpreter")
st.caption("Upload a DESeq2 results CSV, explore summary statistics, enrichment, and LLM-assisted interpretation.")

with st.sidebar:
    st.header("Settings")
    padj_thresh = st.number_input("Adjusted p-value cutoff", min_value=0.0, max_value=1.0, value=0.05, step=0.01)
    lfc_thresh = st.number_input("Absolute log2FC cutoff", min_value=0.0, value=1.0, step=0.1)
    organism = st.selectbox("Organism for enrichment", ["human", "mouse", "yeast", "fly", "fish", "worm"], index=0)
    libraries = st.multiselect(
        "Enrichment libraries",
        options=DEFAULT_GENESET_CHOICES,
        default=["MSigDB_Hallmark_2020", "KEGG_2021_Human"]
    )
    llm_model = st.text_input("Groq model", value=DEFAULT_MODEL)

uploaded_file = st.file_uploader("Upload DESeq2 results CSV", type=["csv"])

if uploaded_file is not None:
    raw_df = pd.read_csv(uploaded_file)
    st.subheader("Input preview")
    st.dataframe(raw_df.head(20), use_container_width=True)

    try:
        df = standardize_deseq2_df(raw_df)
        df = classify_genes(df, padj_thresh=padj_thresh, lfc_thresh=lfc_thresh)
        summary = summarize_results(df)

        st.subheader("Summary")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total genes", summary["n_total"])
        c2.metric("Significant", summary["n_significant"])
        c3.metric("Upregulated", summary["n_up"])
        c4.metric("Downregulated", summary["n_down"])

        st.subheader("Volcano plot")
        volcano_fig = make_volcano_plot(df, padj_thresh=padj_thresh, lfc_thresh=lfc_thresh)
        st.plotly_chart(volcano_fig, use_container_width=True)

        up_genes = df.loc[df["direction"] == "Up", "gene"].tolist()
        down_genes = df.loc[df["direction"] == "Down", "gene"].tolist()

        st.subheader("Top genes")
        left, right = st.columns(2)
        with left:
            st.markdown("**Top upregulated**")
            st.dataframe(summary["top_up"][["gene", "log2FoldChange", "padj"]], use_container_width=True)
        with right:
            st.markdown("**Top downregulated**")
            st.dataframe(summary["top_down"][["gene", "log2FoldChange", "padj"]], use_container_width=True)

        if libraries:
            st.subheader("Pathway enrichment")
            with st.spinner("Running enrichment for upregulated and downregulated genes..."):
                up_enrich = run_enrichment(up_genes, libraries=libraries, organism=organism)
                down_enrich = run_enrichment(down_genes, libraries=libraries, organism=organism)

            e1, e2 = st.columns(2)
            with e1:
                st.markdown("**Upregulated gene enrichment**")
                if up_enrich.empty:
                    st.info("Not enough upregulated genes or no enrichment results returned.")
                else:
                    st.dataframe(up_enrich.head(20), use_container_width=True)
                    st.download_button(
                        "Download up enrichment CSV",
                        dataframe_to_csv_download(up_enrich),
                        file_name="upregulated_enrichment.csv",
                        mime="text/csv"
                    )

            with e2:
                st.markdown("**Downregulated gene enrichment**")
                if down_enrich.empty:
                    st.info("Not enough downregulated genes or no enrichment results returned.")
                else:
                    st.dataframe(down_enrich.head(20), use_container_width=True)
                    st.download_button(
                        "Download down enrichment CSV",
                        dataframe_to_csv_download(down_enrich),
                        file_name="downregulated_enrichment.csv",
                        mime="text/csv"
                    )
        else:
            up_enrich = pd.DataFrame()
            down_enrich = pd.DataFrame()

        params = {
            "padj_threshold": padj_thresh,
            "abs_log2fc_threshold": lfc_thresh,
            "organism": organism,
            "enrichment_libraries": libraries,
        }

        analysis_payload = build_analysis_payload(
            summary=summary,
            top_up_table=summary["top_up"],
            top_down_table=summary["top_down"],
            up_enrich=up_enrich,
            down_enrich=down_enrich,
            params=params,
        )
        st.session_state.analysis_payload = analysis_payload

        st.subheader("LLM interpretation")
        if st.button("Generate report", type="primary"):
            with st.spinner("Generating interpretation..."):
                system_prompt = build_system_prompt()
                user_prompt = build_default_user_prompt(analysis_payload)
                report = call_groq(
                    messages=[
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": user_prompt},
                    ],
                    model=llm_model,
                    temperature=0.2,
                )
                st.session_state.analysis_report = report
                st.session_state.chat_history = []

        if st.session_state.analysis_report:
            st.markdown(st.session_state.analysis_report)

        st.subheader("Ask your own question")
        user_question = st.text_area(
            "Example: Are the results more suggestive of inflammatory signaling or cell-cycle suppression?",
            height=100
        )

        if st.button("Ask question"):
            if not st.session_state.analysis_payload:
                st.warning("Please upload and process a dataset first.")
            elif not user_question.strip():
                st.warning("Enter a question first.")
            else:
                with st.spinner("Generating answer..."):
                    base_context = (
                        "Use only the analysis payload below to answer questions about this dataset.\n\n"
                        f"{json.dumps(st.session_state.analysis_payload, indent=2)}"
                    )

                    messages = [{"role": "system", "content": build_system_prompt()}]
                    messages.append({"role": "user", "content": base_context})

                    for turn in st.session_state.chat_history:
                        messages.append({"role": "user", "content": turn["user"]})
                        messages.append({"role": "assistant", "content": turn["assistant"]})

                    messages.append({"role": "user", "content": user_question})

                    answer = call_groq(
                        messages=messages,
                        model=llm_model,
                        temperature=0.2,
                    )
                    st.session_state.chat_history.append(
                        {"user": user_question, "assistant": answer}
                    )

        if st.session_state.chat_history:
            for i, turn in enumerate(st.session_state.chat_history, start=1):
                with st.expander(f"Q{i}: {turn['user']}", expanded=(i == len(st.session_state.chat_history))):
                    st.markdown(turn["assistant"])

        st.subheader("Processed data")
        st.dataframe(df.head(50), use_container_width=True)
        st.download_button(
            "Download processed DE table",
            dataframe_to_csv_download(df),
            file_name="processed_deseq2_results.csv",
            mime="text/csv"
        )

    except Exception as e:
        st.error(f"Error while processing file: {e}")
else:
    st.info("Upload a DESeq2 CSV file to begin.")