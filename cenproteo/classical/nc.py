import networkx as nx
import pandas as pd
from .base import ClassicalBase


class NC(ClassicalBase):
    """
    Calculates Neighborhood Centrality (NC) for nodes in a PPI network.
    Neighborhood Centrality, NC(u), is defined as:
        NC(u) = ∑_{v ∈ N_u} ECC(u, v)
              = ∑_{v ∈ N_u} (z_{u, v} / min(d_u - 1, d_v - 1))
        where,
        - N_u represents the neighborhood of node u,
        - ECC(u, v) is the Edge Clustering Coefficient between nodes u and v,
        - z_{u, v} is the number of common neighbors between nodes u and v,
        - d_u and d_v are the degrees of nodes u and v, respectively.

        reference:
        Wang J, Li M, Wang H, Pan Y.
        Identification of essential proteins based on edge clustering coefficient.
        IEEE/ACM Trans Comput Biol Bioinform. 2011;9(4):1070-80
    Inherits from ClassicalBase for common functionalities.
    """

    def __init__(self, ppi_file):
        """
        Initializes the DC algorithm class.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format).
        """
        super().__init__(ppi_file)  # Initialize the base class

    def calculate(self) -> pd.DataFrame:
        """
        Calculates the Neighborhood Centrality for all proteins.

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by NC score in descending order.
                              Columns: ['Protein', 'Score']
        """
        if not self.graph:
            print("Error: Graph not loaded. Cannot calculate NC.")
            return None

        # Calculate Neighborhood Centrality using networkx
        nc_scores = {}
        for u in self.graph.nodes():
            d_u = self.graph.degree(u)
            neighbors_u = list(self.graph.neighbors(u))
            NC_u = 0
            # Iterate through neighbors of u
            # and calculate common neighbors with each neighbor v
            for v in neighbors_u:
                common_neighbors = list(nx.common_neighbors(self.graph, u, v))
                z_uv = len(common_neighbors)
                d_v = self.graph.degree(v)
                # Calculate ECC(u, v)
                # and add to NC(u) if min(d_u - 1, d_v - 1) > 0
                if min(d_u - 1, d_v - 1) > 0:
                    ECC_uv = z_uv / min(d_u - 1, d_v - 1)
                    NC_u += ECC_uv
            nc_scores[u] = NC_u

        # Convert to DataFrame
        df_scores = pd.DataFrame(list(nc_scores.items()), columns=["Protein", "Score"])

        # Sort by score
        sorted_df = df_scores.sort_values(by="Score", ascending=False).reset_index(
            drop=True
        )

        self.results = sorted_df  # Store results internally
        return sorted_df
