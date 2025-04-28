import networkx as nx
import pandas as pd
from .base import ClassicalBase


class IC(ClassicalBase):
    """
    Calculates Information Centrality (IC) for nodes in a PPI network.
    Inherits from ClassicalBase for common functionalities.
    """

    def __init__(self, ppi_file):
        """
        Initializes the IC algorithm class.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format).
        """
        super().__init__(ppi_file)  # Initialize the base class

    def calculate(self) -> pd.DataFrame:
        """
        Calculates the Information Centrality for all proteins.

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by IC score in descending order.
                              Columns: ['Protein', 'Score']
        """
        if not self.graph:
            print("Error: Graph not loaded. Cannot calculate IC.")
            return None

        # Calculate Information Centrality using networkx
        # Normalize by N-1 for classical definition if graph is not empty
        n_nodes = len(self.graph)
        if n_nodes > 1:
            largest_cc = max(
                nx.connected_components(self.graph), key=len
            )  # compute the max connected components
            subgraph = self.graph.subgraph(
                largest_cc
            )  # construct the max connected subgraph and compute its information_centrality
            ic_scores = nx.current_flow_betweenness_centrality(subgraph)
        elif n_nodes == 1:  # Handle single node graph
            node = list(self.graph.nodes())[0]
            ic_scores = {node: 0.0}
        else:  # Handle empty graph
            ic_scores = {}

        # Convert to DataFrame
        df_scores = pd.DataFrame(list(ic_scores.items()), columns=["Protein", "Score"])

        # Sort by score
        sorted_df = df_scores.sort_values(by="Score", ascending=False).reset_index(
            drop=True
        )

        self.results = sorted_df  # Store results internally
        return sorted_df
