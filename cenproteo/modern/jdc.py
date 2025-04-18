import pandas as pd
from .base import ModernBase
import networkx as nx


class JDC(ModernBase):
    """
    Calculates Joint Degree Centrality (JDC) using PPI network and gene expression data.
    Inherits from ModernBase for common functionalities.
    Reference: Zhong, J., Tang, C., Peng, W. et al. BMC Bioinformatics 22, 248 (2021).
    """

    def __init__(self, ppi_file, gene_expression_file):
        """
        Initializes the JDC algorithm class.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format).
            gene_expression_file (str): Path to the gene expression file (CSV format).
                                        Expected first column as gene ID, rest as expression values,
                                        with optional 'mean' and 'std' columns at the end.
        """
        # Pass gene_expression_file to the base class initializer
        super().__init__(ppi_file, gene_expression_file=gene_expression_file)
        if self.gene_expression_data is None:
            raise ValueError("Gene expression data is required for JDC calculation.")

    def _calculate_ecc(self, u, v):
        """Calculates Edge Clustering Coefficient (ECC) for an edge (u, v)."""
        # Use nx.common_neighbors for potentially better performance/readability
        common_neighbors = list(nx.common_neighbors(self.graph, u, v))
        z_uv = len(common_neighbors)
        try:
            # Get degrees, handle potential KeyError if node not in graph (shouldn't happen for edges)
            degree_u = self.graph.degree[u]
            degree_v = self.graph.degree[v]
        except KeyError as e:
            # Log a warning if a node is unexpectedly missing
            print(
                f"Warning: Node {e} not found in graph during ECC calculation for edge ({u}, {v})."
            )
            return 0.0

        # Calculate the denominator for ECC
        min_degree_minus_1 = min(degree_u - 1, degree_v - 1)

        # Avoid division by zero
        if min_degree_minus_1 <= 0:
            return 0.0
        else:
            # Calculate ECC according to the reference paper's formula
            return z_uv / min_degree_minus_1

    def _calculate_jaccard(self, u, v, active_profile):
        """Calculates Jaccard similarity based on active time points."""
        # Retrieve the sets of active time point indices for proteins u and v
        profile_u = active_profile.get(u, set())
        profile_v = active_profile.get(v, set())

        # Check if profiles are valid sets (basic type check)
        if not isinstance(profile_u, set) or not isinstance(profile_v, set):
            print(
                f"Warning: Invalid active profile type for Jaccard calculation ({u}, {v}). Expected set."
            )
            return 0.0

        # Calculate intersection and union sizes
        intersection_size = len(profile_u & profile_v)
        union_size = len(profile_u | profile_v)

        # Avoid division by zero if the union is empty
        if union_size == 0:
            return 0.0
        else:
            # Calculate Jaccard index
            return intersection_size / union_size

    def calculate(self):
        """
        Calculates the JDC score for all proteins.

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by JDC score in descending order.
                              Columns: ['Protein', 'Score']
        """
        # Ensure graph and expression data are loaded
        if not self.graph or self.gene_expression_data is None:
            print("Error: Graph and Gene Expression data must be loaded.")
            return None

        print("Calculating JDC scores...")
        # Initialize scores dictionary for all nodes in the graph
        jdc_scores = {node: 0.0 for node in self.graph.nodes()}
        # Dictionary to store the set of active time point indices for each protein
        active_profile = {}

        # --- Pre-calculate active profiles based on gene expression volatility ---
        print("Calculating gene expression volatility and active profiles...")
        # Find proteins present in both the graph and the expression data
        common_proteins = list(
            set(self.graph.nodes()) & set(self.gene_expression_data.index)
        )
        # Handle case where there is no overlap
        if not common_proteins:
            print(
                "Warning: No overlap between proteins in PPI and gene expression data."
            )
            # Return DataFrame with zero scores for all graph nodes
            self.results = pd.DataFrame(
                list(jdc_scores.items()), columns=["Protein", "Score"]
            )
            return self.results.sort_values(by="Score", ascending=False).reset_index(
                drop=True
            )

        # Identify actual expression columns (exclude 'mean', 'std' if they exist by name)
        expr_cols = self.gene_expression_data.columns
        data_cols = [col for col in expr_cols if col.lower() not in ["mean", "std"]]

        # Check if data columns were correctly identified
        if len(data_cols) == len(expr_cols):
            print(
                "Info: 'mean' and 'std' columns not found or not excluded by name. Assuming all columns are expression data."
            )
        elif len(data_cols) == 0:
            raise ValueError(
                "No expression data columns found after excluding 'mean' and 'std'. Check column names."
            )
        else:
            print(f"Info: Using columns {data_cols} for expression data.")

        # Filter the DataFrame to include only common proteins and actual data columns
        expr_data_filtered = self.gene_expression_data.loc[common_proteins, data_cols]

        # Calculate mean and std dev based *only* on the actual expression data columns
        means = expr_data_filtered.mean(axis=1)
        stds = expr_data_filtered.std(axis=1)
        # Handle potential NaN std dev (e.g., if only one time point) by replacing with 0
        stds = stds.fillna(0)

        # Calculate volatility and threshold G(u)
        # If std=0, volatility becomes 1
        volatility = 1 / (1 + stds**2)
        # If std=0, the second term becomes 0, threshold = mean
        thresholds = means + 2 * stds * volatility
        # Convert thresholds to dictionary for faster scalar lookup
        thresholds_dict = thresholds.to_dict()

        # Determine active time points using only expression data columns
        num_time_points = len(data_cols)  # Number of actual data columns
        # Get numpy array for potentially faster row/column access
        expr_data_values = expr_data_filtered.values

        # Iterate through common proteins to determine their active time points
        for i, protein in enumerate(
            common_proteins
        ):  # Use enumerate for index access to numpy array
            threshold_g = thresholds_dict[
                protein
            ]  # Get the specific threshold for this protein
            active_points = set()
            protein_expression = expr_data_values[
                i, :
            ]  # Get the expression row for this protein
            # Iterate through each time point
            for t in range(num_time_points):
                # Compare expression value at time t with the threshold
                if protein_expression[t] > threshold_g:
                    active_points.add(t)  # Add the index 't' of the active time point
            active_profile[protein] = active_points  # Store the set of active indices
        print("Active profiles calculated.")

        # --- Calculate JDC for each protein ---
        print("Calculating ECC and Jaccard for edges...")
        processed_edges = 0
        total_edges = len(self.graph.edges())
        # Optional: Cache ECC values to avoid recalculation for symmetric edges
        ecc_cache = {}

        # Iterate through all edges in the graph
        for u, v in self.graph.edges():
            # Process only if both proteins have expression data (are in active_profile)
            if u in active_profile and v in active_profile:
                try:
                    # Calculate or retrieve ECC from cache
                    # Check both (u, v) and (v, u) in cache for symmetry
                    edge_key = tuple(
                        sorted((u, v))
                    )  # Use sorted tuple as consistent key
                    if edge_key in ecc_cache:
                        ecc_uv = ecc_cache[edge_key]
                    else:
                        ecc_uv = self._calculate_ecc(u, v)
                        ecc_cache[edge_key] = ecc_uv

                    # Calculate Jaccard similarity using the active profiles
                    jaccard_uv = self._calculate_jaccard(u, v, active_profile)

                    # Calculate JDC contribution for the edge
                    jc_uv = jaccard_uv * ecc_uv

                    # Add the contribution to the JDC scores of both nodes
                    jdc_scores[u] += jc_uv
                    jdc_scores[v] += jc_uv
                except Exception as e:
                    # Catch and report potential errors during calculation for a specific edge
                    print(
                        f"Warning: Error calculating JDC contribution for edge ({u}, {v}): {e}"
                    )

            # Update progress counter
            processed_edges += 1
            # Report progress periodically
            if (
                processed_edges % 5000 == 0 or processed_edges == total_edges
            ):  # Report every 5000 edges and at the end
                print(f"Processed {processed_edges}/{total_edges} edges...")

        print("JDC calculation finished.")
        # Convert the scores dictionary to a DataFrame
        df_scores = pd.DataFrame(list(jdc_scores.items()), columns=["Protein", "Score"])

        # Handle potential NaN or Inf values (though unlikely for JDC) before sorting
        df_scores.replace([float("inf"), -float("inf")], float("nan"), inplace=True)
        df_scores.dropna(subset=["Score"], inplace=True)

        # Sort the DataFrame by score in descending order
        sorted_df = df_scores.sort_values(by="Score", ascending=False).reset_index(
            drop=True
        )

        self.results = sorted_df  # Store the final sorted results internally
        return sorted_df
