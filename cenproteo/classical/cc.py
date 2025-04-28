import networkx as nx
import pandas as pd
from .base import ClassicalBase


class CC(ClassicalBase):
    """
    Calculates Closeness Centrality (CC) for nodes in a PPI network.
    Inherits from ClassicalBase for common functionalities.
    """

    def __init__(self, ppi_file):
        """
        Initializes the CC algorithm class.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format).
        """
        super().__init__(ppi_file)  # Initialize the base class

    def calculate(self) -> pd.DataFrame:
        """
        Calculates the Closeness Centrality for all proteins.

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by CC score in descending order.
                              Columns: ['Protein', 'Score']
        """
        if not self.graph:
            print("Error: Graph not loaded. Cannot calculate CC.")
            return None

        # Calculate Closeness Centrality using networkx
        # Normalize by N-1 for classical definition if graph is not empty
        n_nodes = len(self.graph)
        if n_nodes > 1:
            cc_scores = nx.closeness_centrality(self.graph)
        elif n_nodes == 1:  # Handle single node graph
            node = list(self.graph.nodes())[0]
            cc_scores = {node: 0.0}
        else:  # Handle empty graph
            cc_scores = {}

        # Convert to DataFrame
        df_scores = pd.DataFrame(list(cc_scores.items()), columns=["Protein", "Score"])

        # Sort by score
        sorted_df = df_scores.sort_values(by="Score", ascending=False).reset_index(
            drop=True
        )

        self.results = sorted_df  # Store results internally
        return sorted_df
