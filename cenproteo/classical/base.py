import pandas as pd
import networkx as nx
from abc import ABC, abstractmethod


class ClassicalBase(ABC):
    """
    Base class for classical centrality algorithms.
    Handles common tasks like graph loading, result export, and comparison.
    """

    def __init__(self, ppi_file):
        """
        Initializes the base class for classical algorithms.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format).
                            Expected columns: 'Protein A', 'Protein B'.
        """
        self.ppi_file = ppi_file
        self.graph = self._load_graph()
        self.results = None  # To store the calculated scores

    def _load_graph(self):
        """Loads the PPI network from a CSV file into a NetworkX graph."""
        try:
            df = pd.read_csv(self.ppi_file)
            # Assuming the first two columns are the interacting proteins
            graph = nx.from_pandas_edgelist(
                df, source=df.columns[0], target=df.columns[1]
            )
            print(
                f"Graph loaded successfully: {len(graph.nodes())} nodes, {len(graph.edges())} edges."
            )
            return graph
        except FileNotFoundError:
            print(f"Error: PPI file not found at {self.ppi_file}")
            raise
        except Exception as e:
            print(f"Error loading graph from {self.ppi_file}: {e}")
            raise

    @abstractmethod
    def calculate(self):
        """
        Abstract method to calculate the specific centrality score.
        Must be implemented by subclasses.

        Returns:
            pandas.DataFrame: DataFrame with proteins sorted by score.
                              Columns: ['Protein', 'Score']
        """
        pass

    def export_results_to_csv(self, output_file):
        """
        Exports the calculated centrality scores to a CSV file.

        Args:
            output_file (str): Path to save the results CSV file.
        """
        if self.results is None:
            print("Error: No results to export. Please run calculate() first.")
            return
        try:
            self.results.to_csv(output_file, index=False)
            print(f"Results successfully exported to {output_file}")
        except Exception as e:
            print(f"Error exporting results to {output_file}: {e}")

    def _load_gold_standard(self, gold_standard_file):
        """Loads the gold standard essential proteins from a CSV file."""
        try:
            df_gold = pd.read_csv(gold_standard_file)
            # Assume essential proteins can be in any of the first few columns
            # Combine all potential columns into a single set for efficient lookup
            essential_proteins = set()
            for col in df_gold.columns:
                # Check if column likely contains protein names (simple check)
                if df_gold[col].dtype == "object":
                    essential_proteins.update(
                        df_gold[col].dropna().astype(str).unique()
                    )
            if not essential_proteins:
                print(
                    f"Warning: No essential proteins found in {gold_standard_file}. Check file format."
                )
            return essential_proteins
        except FileNotFoundError:
            print(f"Error: Gold standard file not found at {gold_standard_file}")
            raise
        except Exception as e:
            print(f"Error loading gold standard file {gold_standard_file}: {e}")
            raise

    def compare_top_n(self, n, gold_standard_file) -> int:
        """
        Compares the top N predicted proteins with the gold standard essential proteins.

        Args:
            n (int): The number of top proteins to consider.
            gold_standard_file (str): Path to the gold standard essential proteins CSV file.

        Returns:
            int: The number of correctly identified essential proteins in the top N.
                 Returns -1 if results or gold standard are not available.
        """
        if self.results is None:
            print("Error: No results to compare. Please run calculate() first.")
            return -1

        try:
            essential_proteins = self._load_gold_standard(gold_standard_file)
            if not essential_proteins:
                return 0  # No gold standard proteins to compare against

            # Get the top N predicted proteins
            top_n_predicted = self.results["Protein"].head(n).tolist()

            # Count the intersection
            correct_count = len(set(top_n_predicted) & essential_proteins)
            print(f"Comparison: Top {n} predicted proteins vs Gold Standard.")
            print(f"Number of correctly predicted essential proteins: {correct_count}")
            return correct_count

        except Exception as e:
            print(f"Error during comparison: {e}")
            return -1
