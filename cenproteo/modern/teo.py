import pandas as pd
import networkx as nx
from .base import ModernBase
import numpy as np


class TEO(ModernBase):
    """
    Calculates Topology-Expression-Ontology (TEO) score using PPI network (with GO similarity),
    and gene expression data.
    Inherits from ModernBase for common functionalities.
    Reference: Zhang W, Xu J, Li Y, et al. IEEE/ACM transactions on computational biology and bioinformatics, 2016, 15(1): 109-116.
    """

    def __init__(self, ppi_file, gene_expression_file):
        """
        Initializes the TEO algorithm class.

        Args:
            ppi_file (str): Path to the PPI network file which also includes GO similarity info (CSV format).
                            Expected columns like: 'Protein A', 'Protein B', 'GO similarity under BP term',
                            'GO similarity under MF term', 'GO similarity under CC term'. Column order might be used as fallback.
            gene_expression_file (str): Path to the gene expression file (CSV format).
                                        Expected first column as gene ID, rest as expression values.
        """
        self.ppi_file_path = ppi_file  # Store path for reference if needed

        try:
            # Use pandas for robust CSV reading, handle potential bad lines if necessary
            self.df_ppi = pd.read_csv(ppi_file, on_bad_lines="warn")
            print(
                f"PPI and GO similarity data loaded from {ppi_file}. Shape: {self.df_ppi.shape}"
            )
            if self.df_ppi.empty:
                raise ValueError(f"PPI file {ppi_file} loaded an empty DataFrame.")
        except FileNotFoundError:
            print(f"Error: PPI file not found at {ppi_file}")
            raise
        except Exception as e:
            print(f"Error loading data from {ppi_file}: {e}")
            raise

        # Initialize base class - pass ppi_file=None to avoid reloading graph in base
        # We override _load_graph below to handle graph creation from self.df_ppi
        super().__init__(ppi_file=None, gene_expression_file=gene_expression_file)

        # Manually call the overridden _load_graph
        self.graph = self._load_graph()
        if self.graph is None or len(self.graph.nodes) == 0:
            raise ValueError("Graph could not be loaded or is empty.")

        # Load GO similarity data from the preloaded DataFrame
        self.go_similarity_data = self._load_go_similarity_from_df(self.df_ppi)

        # Check if data loading was successful
        if self.gene_expression_data is None:
            print(
                "Warning: Gene expression data not loaded. PCC calculation will return 0."
            )
        if self.go_similarity_data is None:
            raise ValueError(
                "GO similarity data is required for TEO calculation and could not be loaded from the PPI file."
            )

    def _load_graph(self):
        """Overrides ModernBase._load_graph to use the preloaded self.df_ppi."""
        if self.df_ppi is None or self.df_ppi.empty:
            print("Error: PPI DataFrame not loaded or empty.")
            return None
        try:
            # Define expected column names for edges
            source_col = "Protein A"
            target_col = "Protein B"
            # Check if standard columns exist, otherwise fallback to first two
            if (
                source_col not in self.df_ppi.columns
                or target_col not in self.df_ppi.columns
            ):
                if len(self.df_ppi.columns) >= 2:
                    source_col = self.df_ppi.columns[0]
                    target_col = self.df_ppi.columns[1]
                    print(
                        f"Warning: Columns 'Protein A' or 'Protein B' not found. Using columns '{source_col}' and '{target_col}' for graph edges."
                    )
                else:
                    raise ValueError(
                        "PPI DataFrame has less than 2 columns, cannot determine edges."
                    )

            # Drop rows where source or target protein is NaN/missing
            df_clean = self.df_ppi.dropna(subset=[source_col, target_col])
            if len(df_clean) < len(self.df_ppi):
                print(
                    f"Warning: Dropped {len(self.df_ppi) - len(df_clean)} rows with missing protein IDs in edge columns."
                )

            graph = nx.from_pandas_edgelist(
                df_clean, source=source_col, target=target_col
            )
            # Remove self-loops if necessary (usually not desired in PPI)
            graph.remove_edges_from(nx.selfloop_edges(graph))
            print(
                f"Graph created successfully from preloaded data: {len(graph.nodes())} nodes, {len(graph.edges())} edges (after removing self-loops)."
            )
            return graph
        except Exception as e:
            print(f"Error creating graph from preloaded DataFrame: {e}")
            return None

    def _load_go_similarity_from_df(self, df_go):
        """Loads GO similarity data from the provided DataFrame into a nested dictionary."""
        go_dict = {}
        protein_a_col, protein_b_col = None, None
        bp_col, mf_col, cc_col = None, None, None

        try:
            # --- Robust Column Name Detection ---
            cols = df_go.columns.str.lower()  # Case-insensitive matching

            # Protein columns (try specific names first, then positional)
            if "protein a" in cols and "protein b" in cols:
                protein_a_col = df_go.columns[cols.get_loc("protein a")]
                protein_b_col = df_go.columns[cols.get_loc("protein b")]
            elif len(df_go.columns) >= 2:
                protein_a_col = df_go.columns[0]
                protein_b_col = df_go.columns[1]
                print(
                    f"Warning: Using positional columns '{protein_a_col}', '{protein_b_col}' for proteins."
                )
            else:
                raise ValueError("Cannot determine protein columns A and B.")

            # GO term columns (try specific names, then positional)
            term_map = {}
            expected_terms = {
                "bp": "go similarity under bp term",
                "mf": "go similarity under mf term",
                "cc": "go similarity under cc term",
            }
            found_terms_by_name = {}

            for term_key, expected_name in expected_terms.items():
                if expected_name in cols:
                    found_terms_by_name[term_key] = df_go.columns[
                        cols.get_loc(expected_name)
                    ]

            if len(found_terms_by_name) == 3:  # Found all by name
                bp_col = found_terms_by_name["bp"]
                mf_col = found_terms_by_name["mf"]
                cc_col = found_terms_by_name["cc"]
                print("Found GO term columns by name.")
            elif (
                len(df_go.columns) >= 5
            ):  # Fallback to positions 2, 3, 4 if names failed
                bp_col = df_go.columns[2]
                mf_col = df_go.columns[3]
                cc_col = df_go.columns[4]
                print(
                    f"Warning: Using positional columns '{bp_col}', '{mf_col}', '{cc_col}' for GO terms BP, MF, CC respectively."
                )
            else:
                raise ValueError(
                    f"Insufficient columns in PPI file for GO similarity. Expected at least 5, found {len(df_go.columns)}."
                )

            # --- Data Extraction ---
            required_data_cols = [protein_a_col, protein_b_col, bp_col, mf_col, cc_col]

            for _, row in df_go.iterrows():
                u_raw, v_raw = row[protein_a_col], row[protein_b_col]

                # Convert protein IDs to string and handle potential NaN/None
                u = str(u_raw) if pd.notna(u_raw) else None
                v = str(v_raw) if pd.notna(v_raw) else None
                if u is None or v is None:
                    continue  # Skip rows with missing protein IDs

                # Use .get(col, 0.0) for safer access in case column exists but value is missing
                bp_sim_raw = row.get(bp_col, 0.0)
                mf_sim_raw = row.get(mf_col, 0.0)
                cc_sim_raw = row.get(cc_col, 0.0)  # Use 'CC' key

                # Ensure values are numeric, default to 0.0 if conversion fails or NaN
                bp_sim = pd.to_numeric(bp_sim_raw, errors="coerce")
                mf_sim = pd.to_numeric(mf_sim_raw, errors="coerce")
                cc_sim = pd.to_numeric(cc_sim_raw, errors="coerce")

                bp_sim = 0.0 if pd.isna(bp_sim) else bp_sim
                mf_sim = 0.0 if pd.isna(mf_sim) else mf_sim
                cc_sim = 0.0 if pd.isna(cc_sim) else cc_sim

                # Store similarity for both (u, v) and (v, u) pairs for easy access
                # Ensure nested dict exists
                go_dict.setdefault((u, v), {})["BP"] = bp_sim
                go_dict.setdefault((u, v), {})["MF"] = mf_sim
                go_dict.setdefault((u, v), {})["CC"] = cc_sim  # Use 'CC' as the key

                go_dict.setdefault((v, u), {})["BP"] = bp_sim
                go_dict.setdefault((v, u), {})["MF"] = mf_sim
                go_dict.setdefault((v, u), {})["CC"] = cc_sim  # Use 'CC' as the key

            print(
                f"GO similarity data extracted for {len(go_dict)//2} unique pairs (symmetric)."
            )  # Approx unique pairs
            if not go_dict:
                print("Warning: GO similarity dictionary is empty after processing.")
                return None
            return go_dict

        except Exception as e:
            print(f"Error processing GO similarity from DataFrame: {e}")
            import traceback

            traceback.print_exc()
            return None

    def _calculate_ecc(self, u, v):
        """
        Calculates Edge Clustering Coefficient (ECC) for an edge (u, v).
        Formula: ECC = (N_uv)^3 / min(Degree(u)-1, Degree(v)-1)
        where N_uv is the number of triangles involving edge (u,v).
        """
        if u not in self.graph or v not in self.graph:
            # This check might be redundant if called only within calculate() loop
            # print(f"Warning: Node {u} or {v} not found in graph during ECC calculation.")
            return 0.0
        # Check if the edge (u, v) actually exists in the graph (might have been removed)
        if not self.graph.has_edge(u, v):
            return 0.0

        try:
            # Common neighbors directly gives the number of triangles for the edge
            # Need to handle potential generator exhaustion if iterating neighbors
            # common_neighbors = list(nx.common_neighbors(self.graph, u, v)) # Less efficient for just count
            # More efficient way to get count:
            neighbors_u = set(self.graph.neighbors(u))
            neighbors_v = set(self.graph.neighbors(v))
            # Common neighbors are nodes forming triangles with edge (u,v)
            # Exclude u and v themselves if they were somehow included (shouldn't be)
            num_common_neighbors = len((neighbors_u & neighbors_v) - {u, v})

            degree_u = self.graph.degree(u)
            degree_v = self.graph.degree(v)

            # Denominator: min(Degree(u)-1, Degree(v)-1)
            # Handle cases where degree is 1 (denominator would be 0)
            min_degree_minus_1 = min(degree_u - 1, degree_v - 1)

            if min_degree_minus_1 <= 0:
                # If min degree is 1, ECC is typically considered 0.
                # Also handles division by zero.
                return 0.0
            else:
                # Formula: (N_uv)^3 / min(du - 1, dv - 1)
                ecc_value = (num_common_neighbors**3) / min_degree_minus_1
                # Handle potential floating point inaccuracies if needed, though usually fine
                return ecc_value
        except Exception as e:
            # print(f"Error calculating ECC for edge ({u}, {v}): {e}")
            return 0.0

    def _calculate_pcc(self, u, v):
        """
        Calculates the SIGNED Pearson Correlation Coefficient (PCC) between two proteins
        based on their gene expression profiles. Returns 0 if data is missing,
        insufficient, or standard deviation is zero.
        Follows formula (2) and usage in formula (4) of the reference paper.
        """
        if self.gene_expression_data is None:
            return 0.0
        # Check if proteins exist in the expression data index
        if (
            u not in self.gene_expression_data.index
            or v not in self.gene_expression_data.index
        ):
            # print(f"Warning: Protein {u} or {v} not found in gene expression data.")
            return 0.0

        try:
            # Select expression vectors for u and v
            expr_u = self.gene_expression_data.loc[u]
            expr_v = self.gene_expression_data.loc[v]

            # Ensure data is numeric, converting non-numeric to NaN
            data_u = pd.to_numeric(expr_u, errors="coerce")
            data_v = pd.to_numeric(expr_v, errors="coerce")

            # Create a DataFrame from the two series to easily drop rows with NaN in *either* series
            combined = pd.DataFrame({"u": data_u, "v": data_v}).dropna()

            # Check if sufficient data remains after dropping NaNs
            if len(combined) < 2:
                # Need at least 2 data points for correlation
                return 0.0

            # Extract the cleaned series
            clean_u = combined["u"]
            clean_v = combined["v"]

            # Check for zero standard deviation in the cleaned data
            std_u = clean_u.std()
            std_v = clean_v.std()
            if std_u == 0 or std_v == 0 or pd.isna(std_u) or pd.isna(std_v):
                # Correlation is undefined or meaningless if variance is zero
                return 0.0

            # Calculate SIGNED PCC using pandas .corr()
            # The result of series.corr(other_series) is a single float value
            pcc = clean_u.corr(clean_v, method="pearson")

            # Handle potential NaN result from .corr() (though std checks should prevent this)
            if pd.isna(pcc):
                # print(f"Warning: PCC calculation resulted in NaN for ({u}, {v}) despite checks.")
                return 0.0

            # Return the signed PCC value as per the paper's formula usage
            return pcc

        except KeyError:
            # This might happen if index lookup fails unexpectedly
            # print(f"Warning: Key error during PCC calculation for ({u}, {v}).")
            return 0.0
        except Exception as e:
            # Catch other potential calculation errors
            # print(f"Warning: Could not calculate PCC for ({u}, {v}): {e}")
            return 0.0

    def _get_go_similarity(self, u, v, go_term):
        """
        Retrieves the GO similarity score for a given pair (u, v) and term ('BP', 'MF', 'CC').
        Accesses the pre-loaded symmetric dictionary.
        """
        # Access using the (u, v) tuple directly, as the dict is symmetric
        # Default to an empty dict if the pair (u,v) isn't found, then default to 0.0 if the term isn't found
        return self.go_similarity_data.get((u, v), {}).get(go_term, 0.0)

    def calculate(self, go_term="BP"):
        """
        Calculates the TEO score for all proteins based on the specified GO term.

        Args:
            go_term (str): The Gene Ontology term to use ('BP', 'MF', or 'CC'). Defaults to 'BP'.

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by TEO score in descending order.
                              Columns: ['Protein', 'Score']. Returns None on failure.
        """
        if not self.graph:
            print("Error: Graph not loaded. Cannot calculate TEO.")
            return None
        # Standardize valid GO terms to 'BP', 'MF', 'CC'
        valid_go_terms = ["BP", "MF", "CC"]
        if go_term not in valid_go_terms:
            print(
                f"Error: Invalid GO term '{go_term}'. Must be one of {valid_go_terms}."
            )
            return None
        if self.go_similarity_data is None:
            print("Error: GO similarity data not loaded. Cannot calculate TEO.")
            return None

        print(f"Calculating TEO scores using GO term: {go_term}...")
        teo_scores = {node: 0.0 for node in self.graph.nodes()}

        processed_edges = 0
        total_edges = len(self.graph.edges())
        # Cache PCC values to avoid redundant calculations for each edge direction
        pcc_cache = {}

        # Iterate through edges for efficiency, calculating contributions to both nodes
        for u, v in self.graph.edges():
            # Calculate components for the edge (u, v) ONCE
            ecc_uv = self._calculate_ecc(u, v)
            go_sim_uv = self._get_go_similarity(
                u, v, go_term
            )  # Uses the specified go_term

            # Calculate or retrieve SIGNED PCC from cache
            pair_key = tuple(sorted((u, v)))  # Use sorted tuple for cache key
            if pair_key in pcc_cache:
                pcc_uv = pcc_cache[pair_key]
            else:
                pcc_uv = self._calculate_pcc(u, v)
                # Only cache valid numerical PCC values
                if pd.notna(pcc_uv):
                    pcc_cache[pair_key] = pcc_uv
                else:
                    pcc_uv = 0.0  # Default to 0 if calculation failed

            # Calculate the edge contribution term: Ecc * (GO + PCC)
            # Ensure all components are numeric before calculating contribution
            if not all(
                isinstance(val, (int, float)) for val in [ecc_uv, go_sim_uv, pcc_uv]
            ):
                # print(f"Warning: Non-numeric value encountered for edge ({u}, {v}). ECC={ecc_uv}, GO={go_sim_uv}, PCC={pcc_uv}. Skipping contribution.")
                edge_contribution = 0.0
            elif pd.isna(ecc_uv) or pd.isna(go_sim_uv) or pd.isna(pcc_uv):
                # print(f"Warning: NaN value encountered for edge ({u}, {v}). ECC={ecc_uv}, GO={go_sim_uv}, PCC={pcc_uv}. Skipping contribution.")
                edge_contribution = 0.0
            else:
                edge_contribution = ecc_uv * (go_sim_uv + pcc_uv)

            # Add contribution to both nodes u and v (since TEO sums over neighbors)
            # Lock might be needed if using parallel processing, but not here.
            if u in teo_scores:
                teo_scores[u] += edge_contribution
            if v in teo_scores:
                teo_scores[v] += edge_contribution

            processed_edges += 1
            if processed_edges % 5000 == 0 or processed_edges == total_edges:
                print(f"Processed {processed_edges}/{total_edges} edges...")

        print("TEO calculation finished.")

        # Convert scores to DataFrame
        df_scores = pd.DataFrame(list(teo_scores.items()), columns=["Protein", "Score"])

        # Handle potential NaN or Inf values introduced during summation (less likely but possible)
        df_scores.replace([np.inf, -np.inf], np.nan, inplace=True)
        # Optionally fill NaN scores with 0 or drop them
        # df_scores.fillna(0, inplace=True) # Option 1: Fill NaN with 0
        df_scores.dropna(
            subset=["Score"], inplace=True
        )  # Option 2: Drop proteins with NaN scores

        # Sort proteins by score in descending order
        sorted_df = df_scores.sort_values(by="Score", ascending=False).reset_index(
            drop=True
        )

        self.results = sorted_df  # Store results internally
        print(f"Generated TEO scores for {len(sorted_df)} proteins.")
        return sorted_df
