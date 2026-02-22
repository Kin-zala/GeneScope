import streamlit as st
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(page_title="GeneScope", layout="wide")

st.title("GeneScope 🧬")
st.write("Interactive Gene Expression Visualization & Exploratory Analysis Tool")

st.markdown("""
### Expected File Format
- First column: Gene names  
- Remaining columns: Numeric expression values  
- Rows: Genes  
- Columns: Samples  
""")

# --------------------------------------------------
# DATA INPUT
# --------------------------------------------------

df = None

if st.button("Load Demo Dataset"):
    genes = [f"Gene_{i}" for i in range(2000)]
    samples = [f"Sample_{i}" for i in range(50)]
    demo_data = np.random.poisson(lam=100, size=(2000, 50))
    df = pd.DataFrame(demo_data, index=genes, columns=samples)

uploaded_file = st.file_uploader("Upload your CSV file", type=["csv"])

if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file)
        df = df.set_index(df.columns[0])
        df = df.apply(pd.to_numeric, errors="coerce")
        df = df.dropna(how="all")
    except Exception as e:
        st.error(f"Error reading CSV file: {e}")

# --------------------------------------------------
# CONTINUE IF DATA EXISTS
# --------------------------------------------------

if df is not None:

    st.subheader("Data Preview")
    st.dataframe(df.head())

    st.sidebar.header("Analysis Settings")

    analysis_mode = st.sidebar.radio(
        "Select Analysis Mode:",
        ("Standard Visualization", "PCA for Large Datasets")
    )

    log_transform = st.sidebar.checkbox("Apply Log2(x + 1) Transformation")

    if log_transform:
        try:
            df = np.log2(df + 1)
        except Exception as e:
            st.error(f"Error applying log2 transformation: {e}")

    # ==================================================
    # STANDARD VISUALIZATION MODE
    # ==================================================

    if analysis_mode == "Standard Visualization":

        top_n = st.sidebar.slider("Number of top variable genes", 5, 50, 20)

        try:
            df["Variance"] = df.var(axis=1)
            top_genes = df.sort_values("Variance", ascending=False).head(top_n)
            heatmap_df = top_genes.drop(columns=["Variance"])
        except Exception as e:
            st.error(f"Error computing variance or selecting top genes: {e}")
            heatmap_df = pd.DataFrame()

        # ------------------ BOX PLOT ------------------

        if not heatmap_df.empty:
            st.subheader("Interactive Boxplot")
            box_df = heatmap_df.T
            colors = px.colors.qualitative.Safe  # colorblind friendly

            try:
                fig_box = go.Figure()
                for i, gene in enumerate(box_df.columns):
                    fig_box.add_trace(go.Box(
                        y=box_df[gene],
                        name=gene,
                        boxpoints="outliers",
                        marker_color=colors[i % len(colors)]
                    ))

                fig_box.update_layout(
                    xaxis_title="Genes",
                    yaxis_title="Expression Level",
                    height=600
                )

                st.plotly_chart(fig_box, use_container_width=True)

                img_box = fig_box.to_image(format="png", width=1000, height=600, scale=2)
                st.download_button(
                    "Download Boxplot as PNG",
                    data=img_box,
                    file_name="genescope_boxplot.png",
                    mime="image/png"
                )
            except Exception as e:
                st.error(f"Error generating boxplot: {e}")

        # ------------------ HEATMAP ------------------

        if not heatmap_df.empty:
            st.subheader("Interactive Heatmap")
            try:
                fig_heat = go.Figure(data=go.Heatmap(
                    z=heatmap_df.values,
                    x=heatmap_df.columns,
                    y=heatmap_df.index,
                    colorscale="Purples",
                    text=heatmap_df.round(1).values,
                    texttemplate="%{text}",
                    colorbar=dict(
                        title="Expression (log2)" if log_transform else "Expression"
                    )
                ))

                fig_heat.update_layout(
                    xaxis_title="Samples",
                    yaxis_title="Genes",
                    height=700
                )

                st.plotly_chart(fig_heat, use_container_width=True)

                img_heat = fig_heat.to_image(format="png", width=1000, height=800, scale=2)
                st.download_button(
                    "Download Heatmap as PNG",
                    data=img_heat,
                    file_name="genescope_heatmap.png",
                    mime="image/png"
                )
            except Exception as e:
                st.error(f"Error generating heatmap: {e}")

    # ==================================================
    # PCA MODE
    # ==================================================

    elif analysis_mode == "PCA for Large Datasets":

        st.subheader("PCA + KMeans Clustering")

        df_transposed = df.T

        try:
            # Variance filtering if extremely large
            if df_transposed.shape[1] > 5000:
                variances = df.var(axis=1)
                top_genes = variances.sort_values(ascending=False).head(5000).index
                df_transposed = df.loc[top_genes].T
                st.info("Using top 5000 most variable genes for PCA to improve performance.")

            # Scale and PCA
            scaler = StandardScaler()
            df_scaled = scaler.fit_transform(df_transposed)

            pca = PCA(n_components=2, svd_solver="randomized", random_state=42)
            components = pca.fit_transform(df_scaled)

            pca_df = pd.DataFrame(
                components,
                columns=["PC1", "PC2"],
                index=df_transposed.index
            )

            # KMeans clustering
            st.sidebar.subheader("Clustering Settings")
            n_clusters = st.sidebar.slider("Number of clusters", 2, 6, 3)

            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            pca_df["Cluster"] = kmeans.fit_predict(pca_df).astype(str)

            # Interactive PCA plot
            fig_pca = px.scatter(
                pca_df,
                x="PC1",
                y="PC2",
                color="Cluster",
                hover_name=pca_df.index,
                title="PCA + KMeans Clustering",
                color_discrete_sequence=px.colors.qualitative.Safe
            )

            fig_pca.update_traces(marker=dict(size=8, opacity=0.7))
            fig_pca.update_layout(
                xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance)",
                yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance)"
            )

            st.plotly_chart(fig_pca, use_container_width=True)

            img_pca = fig_pca.to_image(format="png", width=900, height=700, scale=2)
            st.download_button(
                "Download PCA Plot as PNG",
                data=img_pca,
                file_name="genescope_pca_plot.png",
                mime="image/png"
            )

            # CSV download
            csv = pca_df.to_csv().encode("utf-8")

            st.subheader("Cluster Assignment Table")
            st.dataframe(pca_df)

            st.markdown(
                f"**Explained Variance:** "
                f"PC1 = {pca.explained_variance_ratio_[0]*100:.2f}% | "
                f"PC2 = {pca.explained_variance_ratio_[1]*100:.2f}%"
            )
            st.download_button(
                "Download PCA + Cluster Results (CSV)",
                data=csv,
                file_name="genescope_pca_results.csv",
                mime="text/csv"
            )

        except Exception as e:
            st.error(f"PCA or clustering failed: {e}")