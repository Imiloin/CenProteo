import networkx as nx
import pandas as pd
from .base import ClassicalBase


class EC(ClassicalBase):
    """
    Calculates Eigenvector Centrality (EC) for nodes in a PPI network.
    Inherits from ClassicalBase for common functionalities.
    """

    def __init__(self, ppi_file):
        """
        Initializes the EC algorithm class.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format).
        """
        super().__init__(ppi_file)  # Initialize the base class

    def calculate(self) -> pd.DataFrame:
        """
        Calculates the Eigenvector Centrality for all proteins.

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by EC score in descending order.
                              Columns: ['Protein', 'Score']
        """
        if not self.graph:
            print("Error: Graph not loaded. Cannot calculate EC.")
            return None

        # Calculate Eigenvector Centrality using networkx
        # Normalize by N-1 for classical definition if graph is not empty
        n_nodes = len(self.graph)
        if n_nodes > 1:
            ec_scores = nx.eigenvector_centrality(self.graph)
        elif n_nodes == 1:  # Handle single node graph
            node = list(self.graph.nodes())[0]
            ec_scores = {node: 0.0}
        else:  # Handle empty graph
            ec_scores = {}

        # Convert to DataFrame
        df_scores = pd.DataFrame(list(ec_scores.items()), columns=["Protein", "Score"])

        # Sort by score
        sorted_df = df_scores.sort_values(by="Score", ascending=False).reset_index(
            drop=True
        )

        self.results = sorted_df  # Store results internally
        return sorted_df
