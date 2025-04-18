import pandas as pd
import numpy as np
import networkx as nx
from .base import ModernBase
from collections import defaultdict


class TGSO(ModernBase):
    """
    Calculates essentiality scores using the TGSO iteration model.
    Combines PPI network topology, gene expression (PCC, CEN),
    subcellular localization (CLN), and gene orthology (O_score).
    Inherits from ModernBase for common functionalities.
    Reference: Li S, Zhang Z, Li X, et al. BMC bioinformatics, 2021, 22: 1-25.
    """

    def __init__(
        self,
        ppi_file,
        gene_expression_file,
        subcellular_localization_file,
        gene_orthology_file,
        alpha=0.3,
        max_iter=100,
        tol=1e-6,
    ):
        """
        Initializes the TGSO algorithm class.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format, 1st 2 cols are proteins).
            gene_expression_file (str): Path to the gene expression file (CSV format, 1st col=protein ID).
            subcellular_localization_file (str): Path to the subcellular localization file (CSV format).
                                                 Expected columns like 'protein', 'GO_term' (cols 0 and 2 used).
            gene_orthology_file (str): Path to the gene orthology file (CSV format).
                                       Expected columns: 'Protein' (as index or col 0), 'O_score' (col 2 or named).
            alpha (float): Weight parameter for initial score in iteration. Defaults to 0.3.
            max_iter (int): Maximum number of iterations for score calculation. Defaults to 100.
            tol (float): Tolerance for convergence check. Defaults to 1e-6.
        """
        # Initialize base class with PPI and gene expression
        super().__init__(
            ppi_file,
            gene_expression_file,
            subcellular_localization_file=subcellular_localization_file,
            gene_orthology_file=gene_orthology_file,
        )

        # Load additional data using helper methods
        self.localization_data, self.loc_counts, self.total_loc_proteins = (
            self._load_localization_data(subcellular_localization_file)
        )
        self.orthology_data = self._load_orthology_data(
            gene_orthology_file
        )  # Contains O_score

        # Store algorithm parameters
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol

        # Pre-calculate S_scores needed for CLN
        self.sub_scores, self.s_scores = self._precalculate_localization_scores()

        if self.gene_expression_data is None:
            print(
                "Warning: Gene expression data not loaded. PCC/CEN calculation will be affected."
            )
        if self.localization_data is None:
            print(
                "Warning: Subcellular localization data not loaded. CLN calculation will be affected."
            )
        if self.orthology_data is None:
            raise ValueError(
                "Gene orthology data (O_score) is required for TGSO calculation and could not be loaded."
            )

    def _load_localization_data(self, file_path):
        """Loads and processes subcellular localization data."""
        try:
            df_loc = self._load_generic_csv("subcellular_localization", file_path)
            if df_loc is None:
                return None, None, 0

            # Process into a dictionary: protein -> set of GO terms (localizations)
            loc_dict = defaultdict(set)
            loc_counts = defaultdict(int)  # Count proteins per location
            processed_proteins = set()  # Track unique proteins with localization

            # Assuming columns are 'protein' (col 0) and 'GO_term' (col 2) based on example
            protein_col_idx = 0
            go_col_idx = 2

            for _, row in df_loc.iterrows():
                # Ensure protein and GO term are strings and handle potential NaN
                protein = str(row.iloc[protein_col_idx])
                go_term = str(row.iloc[go_col_idx])
                if protein == "nan" or go_term == "nan":
                    continue

                loc_dict[protein].add(go_term)
                # Only count protein for a location once for sub_score calculation
                if protein not in processed_proteins:
                    loc_counts[go_term] += 1  # Increment count for this location
                    processed_proteins.add(
                        protein
                    )  # Mark protein as counted for loc distribution

            total_loc_proteins = sum(
                loc_counts.values()
            )  # Sum of counts across all locations
            print(
                f"Subcellular localization data loaded for {len(loc_dict)} proteins across {len(loc_counts)} locations."
            )
            return dict(loc_dict), dict(loc_counts), total_loc_proteins
        except Exception as e:
            print(f"Error processing subcellular localization file {file_path}: {e}")
            return None, None, 0

    def _load_orthology_data(self, file_path):
        """Loads gene orthology data (O_score)."""
        try:
            # Try loading with index first
            df_ortho = self._load_generic_csv("gene_orthology", file_path, index_col=0)

            if df_ortho is None:  # Try loading without index if first fails
                df_ortho = self._load_generic_csv(
                    "gene_orthology", file_path, index_col=None
                )
                if df_ortho is None:
                    return None
                # Assume Protein is col 0, O_score is col 2 if index load failed
                if df_ortho.shape[1] < 3:
                    raise ValueError(
                        "Insufficient columns for Protein and O_score in orthology file."
                    )
                df_ortho = df_ortho.set_index(df_ortho.columns[0])
                df_ortho["O_score"] = pd.to_numeric(
                    df_ortho.iloc[:, 1], errors="coerce"
                )  # Use second non-index column (index 1)
                df_ortho = df_ortho[["O_score"]]  # Keep only O_score

            elif "O_score" not in df_ortho.columns:
                # If index worked but no 'O_score' column, assume it's the second column (index 1)
                if df_ortho.shape[1] < 2:
                    raise ValueError(
                        "Insufficient columns for O_score in orthology file."
                    )
                print(
                    "Warning: 'O_score' column not found, assuming it's the second column."
                )
                df_ortho["O_score"] = pd.to_numeric(
                    df_ortho.iloc[:, 1], errors="coerce"
                )
                df_ortho = df_ortho[["O_score"]]  # Keep only O_score
            else:
                # Ensure O_score is numeric
                df_ortho["O_score"] = pd.to_numeric(
                    df_ortho["O_score"], errors="coerce"
                )
                df_ortho = df_ortho[["O_score"]]  # Keep only O_score

            # Convert to dictionary: protein -> O_score value
            ortho_dict = df_ortho["O_score"].to_dict()
            # Handle potential NaNs introduced by coercion
            ortho_dict = {k: (v if pd.notna(v) else 0.0) for k, v in ortho_dict.items()}
            print(
                f"Gene orthology data (O_score) loaded for {len(ortho_dict)} proteins."
            )
            return ortho_dict
        except Exception as e:
            print(f"Error processing gene orthology file {file_path}: {e}")
            return None

    def _precalculate_localization_scores(self):
        """Precalculates sub_score(i) and S_score(u) based on Eqs 4 & 5."""
        if self.localization_data is None or self.total_loc_proteins == 0:
            print(
                "Warning: Cannot precalculate localization scores (sub_score, S_score)."
            )
            return {}, {}

        sub_scores = {}  # Store sub_score for each location i
        print("Calculating sub_scores (Eq 4)...")
        for loc, count in self.loc_counts.items():
            sub_scores[loc] = count / self.total_loc_proteins
        print(f"Calculated sub_scores for {len(sub_scores)} locations.")

        s_scores = {}  # Store S_score for each protein u
        print("Calculating S_scores (Eq 5)...")
        # Iterate through proteins present in the graph AND localization data
        proteins_to_calc = set(self.graph.nodes()) & set(self.localization_data.keys())
        for protein_u in proteins_to_calc:
            score = 0.0
            protein_locations = self.localization_data.get(protein_u, set())
            for loc_i in protein_locations:
                score += sub_scores.get(loc_i, 0.0)  # Add sub_score if location exists
            s_scores[protein_u] = score
        print(f"Calculated S_scores for {len(s_scores)} proteins.")
        return sub_scores, s_scores

    def _calculate_pcc(self, u, v):
        """Calculates the absolute Pearson Correlation Coefficient (PCC). Returns 0.0 on error."""
        # Ensure u and v are strings (graph nodes are guaranteed strings now)
        u_str, v_str = str(u), str(v)

        if (
            self.gene_expression_data is None
            or u_str not in self.gene_expression_data.index
            or v_str not in self.gene_expression_data.index
        ):
            return 0.0
        try:
            expr_u = self.gene_expression_data.loc[u_str]
            expr_v = self.gene_expression_data.loc[v_str]

            # Handle potential Series or DataFrame rows, exclude metadata if needed
            # Assuming expression values are numeric or coercible to numeric
            data_u = pd.to_numeric(expr_u, errors="coerce").values
            data_v = pd.to_numeric(expr_v, errors="coerce").values

            # Simple NaN handling: ignore pairs with NaNs
            valid_mask = ~np.isnan(data_u) & ~np.isnan(data_v)
            data_u_valid = data_u[valid_mask]
            data_v_valid = data_v[valid_mask]

            n = len(data_u_valid)
            if n <= 1:  # Need more than 1 point for correlation
                return 0.0

            std_u = np.std(data_u_valid)
            std_v = np.std(data_v_valid)

            if std_u == 0 or std_v == 0:
                return 0.0  # No correlation if one series is constant

            # Calculate PCC using numpy for efficiency
            pcc = np.corrcoef(data_u_valid, data_v_valid)[0, 1]

            # Handle potential NaN result from corrcoef if input was problematic
            return abs(pcc) if pd.notna(pcc) else 0.0

        except Exception as e:
            # print(f"Debug: Error calculating PCC for ({u_str}, {v_str}): {e}") # Optional debug
            return 0.0

    def _calculate_adn(self, u, v):
        """Calculates Aggregation Degree Network (ADN) score (Eq 1)."""
        # Graph ensures u,v exist if they are neighbors
        try:
            neighbors_u = set(self.graph.neighbors(u))
            neighbors_v = set(self.graph.neighbors(v))
            common_neighbors = neighbors_u & neighbors_v
            num_common_neighbors = len(common_neighbors)

            # Use graph.degree for actual neighbors
            degree_u = self.graph.degree(u)
            degree_v = self.graph.degree(v)
            min_degree = min(degree_u, degree_v)

            if min_degree == 0:
                # If min degree is 0, they cannot be neighbors unless graph is empty/single node
                # The +1 term still applies even if common_neighbors is 0
                return 0.0  # Or perhaps 1/0 is undefined? Let's return 0 for non-neighbors/isolated
                # Revisit: Paper formula is (|NG intersect NG| + 1) / min(|NG|, |NG|)
                # If min_degree is 0, denominator is 0. Let's return 0.
            else:
                return (num_common_neighbors + 1) / min_degree
        except (
            nx.NetworkXError
        ):  # Handle cases where u or v might not be in graph (shouldn't happen if called on neighbors)
            return 0.0

    def _calculate_cen(self, u, v, pcc_cache):
        """Calculates Co-Expression Network (CEN) score (Eq 3)."""
        # CEN calculation depends on PCC and common neighbors
        # Get PCC(u, v) - use cache
        pcc_uv = pcc_cache.get((u, v))
        if pcc_uv is None:
            pcc_uv = self._calculate_pcc(u, v)
            pcc_cache[(u, v)] = pcc_uv
            pcc_cache[(v, u)] = pcc_uv

        cen_score = pcc_uv

        try:
            neighbors_u = set(self.graph.neighbors(u))
            neighbors_v = set(self.graph.neighbors(v))
            common_neighbors = neighbors_u & neighbors_v

            for w in common_neighbors:
                # Get PCC(u, w)
                pcc_uw = pcc_cache.get((u, w))
                if pcc_uw is None:
                    pcc_uw = self._calculate_pcc(u, w)
                    pcc_cache[(u, w)] = pcc_uw
                    pcc_cache[(w, u)] = pcc_uw
                # Get PCC(v, w)
                pcc_vw = pcc_cache.get((v, w))
                if pcc_vw is None:
                    pcc_vw = self._calculate_pcc(v, w)
                    pcc_cache[(v, w)] = pcc_vw
                    pcc_cache[(w, v)] = pcc_vw

                cen_score += pcc_uw * pcc_vw
        except nx.NetworkXError:
            pass  # If u or v somehow not in graph, common_neighbors will be empty

        return cen_score

    def _calculate_cln(self, u, v):
        """Calculates Co-Localization Network (CLN) score (Eq 6 - CORRECTED)."""
        if self.localization_data is None:
            return 0.0

        loc_u = self.localization_data.get(u, set())
        loc_v = self.localization_data.get(v, set())

        if not loc_u and not loc_v:
            return 0.0

        common_loc = loc_u & loc_v
        union_loc = loc_u | loc_v

        if not union_loc:  # Should not happen if loc_u or loc_v is non-empty
            return 0.0

        # Overlap fraction part of Eq 6
        overlap_fraction = len(common_loc) / len(union_loc)

        # S_score part of Eq 6 (using precalculated S_scores)
        s_score_u = self.s_scores.get(u, 0.0)
        s_score_v = self.s_scores.get(v, 0.0)

        cln_score = overlap_fraction * (s_score_u + s_score_v) / 2.0
        return cln_score

    def _calculate_lsg(self, node, cache):
        """Calculates Local Summation Global (LSG) score (Eq 7) for a single node."""
        if node in cache["lsg"]:
            return cache["lsg"][node]

        lsg_score = 0.0
        try:
            for neighbor in self.graph.neighbors(node):
                # Ensure neighbor is also a string
                neighbor_str = str(neighbor)

                # Get ADN(node, neighbor)
                adn_key = tuple(
                    sorted((node, neighbor_str))
                )  # Use sorted tuple for cache key
                adn_val = cache["adn"].get(adn_key)
                if adn_val is None:
                    adn_val = self._calculate_adn(node, neighbor_str)
                    cache["adn"][adn_key] = adn_val

                # Get CLN(node, neighbor)
                cln_key = tuple(sorted((node, neighbor_str)))
                cln_val = cache["cln"].get(cln_key)
                if cln_val is None:
                    cln_val = self._calculate_cln(node, neighbor_str)
                    cache["cln"][cln_key] = cln_val

                # Get CEN(node, neighbor)
                cen_key = tuple(sorted((node, neighbor_str)))
                cen_val = cache["cen"].get(cen_key)
                if cen_val is None:
                    # Pass the PCC cache to CEN calculation
                    cen_val = self._calculate_cen(node, neighbor_str, cache["pcc"])
                    cache["cen"][cen_key] = cen_val

                lsg_score += adn_val * (cln_val + cen_val)
        except nx.NetworkXError:
            pass  # Node might not be in graph if data mismatch occurred

        cache["lsg"][node] = lsg_score
        return lsg_score

    def _calculate_initial_vector_p0(self):
        """Initializes the scores P0 for the iteration based on O_score (Eq 11 - CORRECTED)."""
        initial_scores = {}
        if self.orthology_data is None:
            raise ValueError("Orthology data required for P0 initialization.")

        # Ensure all nodes in the graph have an initial score, default to 0 if no O_score found
        for node in self.graph.nodes():
            initial_scores[node] = self.orthology_data.get(node, 0.0)

        # Optional: Normalize P0 so it sums to 1 (like a probability distribution)
        # total_o_score = sum(initial_scores.values())
        # if total_o_score > 0:
        #     initial_scores = {k: v / total_o_score for k, v in initial_scores.items()}
        # The paper doesn't explicitly state normalization for P0, but it's common in PageRank.
        # Let's stick to the direct O_score based on Eq 11.

        print(
            f"Initialized P0 vector for {len(initial_scores)} nodes based on O_score."
        )
        return initial_scores

    def calculate(self):
        """
        Calculates the TGSO scores iteratively (CORRECTED ITERATION).

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by TGSO score in descending order.
                              Columns: ['Protein', 'Score']
        """
        if not self.graph or len(self.graph.nodes()) == 0:
            print("Error: Graph not loaded or empty.")
            return None
        if self.orthology_data is None:  # P0 depends on this
            print("Error: Orthology data not loaded.")
            return None

        print("Starting TGSO calculation...")
        nodes = list(self.graph.nodes())  # Use nodes present in the actual graph
        num_nodes = len(nodes)

        # --- Calculate LSG Scores ---
        print("Calculating LSG scores...")
        # Initialize cache for this run
        cache = {"pcc": {}, "adn": {}, "cen": {}, "cln": {}, "lsg": {}}
        lsg_scores = {}
        nodes_with_neighbors = 0
        for i, node in enumerate(nodes):
            if self.graph.degree(node) > 0:  # Only calculate LSG if node has neighbors
                lsg_scores[node] = self._calculate_lsg(node, cache)
                nodes_with_neighbors += 1
            else:
                lsg_scores[node] = 0.0  # LSG is 0 for isolated nodes

            if (i + 1) % 500 == 0:
                print(f"  Calculated LSG for {i+1}/{num_nodes} nodes...")
        print(f"LSG scores calculated ({nodes_with_neighbors} nodes had neighbors).")

        total_lsg = sum(lsg_scores.values())
        if total_lsg == 0:
            print(
                "Warning: Total LSG score is zero. PCIN matrix will be zero. Iteration might not change scores significantly."
            )
            # Avoid division by zero later, but recognize the implication
            total_lsg_norm = 1.0
        else:
            total_lsg_norm = total_lsg
        print(f"Total LSG: {total_lsg:.4f}")

        # --- Initialize Scores P0 ---
        P0 = self._calculate_initial_vector_p0()

        # --- Iteration ---
        print(
            f"Starting iteration (alpha={self.alpha}, max_iter={self.max_iter}, tol={self.tol})..."
        )
        # Initialize current score vector P. Ensure all graph nodes have a starting score.
        P = {node: P0.get(node, 0.0) for node in nodes}
        iter_time = 0

        for iteration in range(self.max_iter):
            iter_time += 1
            P_new = {}
            max_diff = 0.0
            sum_p_new = 0.0  # For potential normalization check

            # --- Calculate Σ [PCIN(i, j) * P(j)] ---
            # This can be slow. Pre-calculating PCIN matrix is possible but memory intensive.
            # Let's calculate contributions on the fly.

            # Precompute P(j) values for this iteration
            P_j_values = {node_j: P.get(node_j, 0.0) for node_j in nodes}

            for i, protein_i in enumerate(nodes):
                lsg_i = lsg_scores.get(protein_i, 0.0)
                sum_pcin_pj = 0.0

                # Sum contributions from all nodes j
                for protein_j in nodes:
                    lsg_j = lsg_scores.get(protein_j, 0.0)
                    p_j_old = P_j_values.get(protein_j, 0.0)

                    # Calculate PCIN(i, j) based on Eq 8
                    if protein_i == protein_j:
                        pcin_ij = lsg_i / total_lsg_norm
                    else:
                        # Check if i and j are neighbors - paper doesn't restrict PCIN to neighbors
                        # Eq 8 seems global based on LSG scores
                        pcin_ij = min(lsg_i, lsg_j) / total_lsg_norm

                    sum_pcin_pj += pcin_ij * p_j_old

                # Apply Eq 12: Pt+1(i) = (1 - α) * Σ[PCIN(i,j)Pt(j)] + α * P0(i)
                term1 = (1 - self.alpha) * sum_pcin_pj
                term2 = self.alpha * P0.get(protein_i, 0.0)
                p_i_new = term1 + term2

                P_new[protein_i] = p_i_new
                sum_p_new += p_i_new  # Track sum

                # Track difference for convergence check
                diff = abs(p_i_new - P.get(protein_i, 0.0))
                if diff > max_diff:
                    max_diff = diff

            # --- Convergence Check ---
            print(
                f"Iteration {iter_time}: Max score change = {max_diff:.2e}, Sum(P_new) = {sum_p_new:.4f}"
            )
            if max_diff < self.tol:
                print(f"Converged after {iter_time} iterations.")
                P = P_new
                break

            P = P_new  # Update scores for the next iteration

        else:  # Loop finished without break
            print(
                f"Warning: Did not converge within {self.max_iter} iterations. Using scores from last iteration."
            )

        # --- Finalize ---
        print("TGSO calculation finished.")
        # Handle potential NaN/Inf results before sorting
        final_scores = {
            k: (v if pd.notna(v) and np.isfinite(v) else 0.0) for k, v in P.items()
        }
        df_scores = pd.DataFrame(
            list(final_scores.items()), columns=["Protein", "Score"]
        )

        sorted_df = df_scores.sort_values(by="Score", ascending=False).reset_index(
            drop=True
        )

        self.results = sorted_df
        return sorted_df
