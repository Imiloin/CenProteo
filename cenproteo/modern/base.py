import pandas as pd
import networkx as nx
from abc import ABC, abstractmethod


class ModernBase(ABC):
    """
    Base class for modern centrality algorithms using multiple data sources.
    Handles common tasks like data loading, result export, and comparison.
    """

    def __init__(self, ppi_file, gene_expression_file=None, **kwargs):
        """
        Initializes the base class for modern algorithms.

        Args:
            ppi_file (str): Path to the PPI network file (CSV format).
            gene_expression_file (str, optional): Path to the gene expression file. Defaults to None.
            **kwargs: Additional file paths or parameters needed by specific algorithms.
        """
        self.ppi_file = ppi_file
        self.gene_expression_file = gene_expression_file
        self.graph = self._load_graph()
        self.gene_expression_data = (
            self._load_gene_expression() if gene_expression_file else None
        )
        self.results = None  # To store the calculated scores
        # Store other potential data paths if provided
        self.other_data_files = kwargs

    def _load_graph(self):
        """Loads the PPI network from a CSV file into a NetworkX graph."""
        try:
            df = pd.read_csv(self.ppi_file)
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

    def _load_gene_expression(self):
        """Loads gene expression data from a CSV file."""
        if not self.gene_expression_file:
            return None
        try:
            # Assuming the first column is the gene ID/protein name
            df_expr = pd.read_csv(self.gene_expression_file, index_col=0)
            print(f"Gene expression data loaded for {len(df_expr)} genes.")
            return df_expr
        except FileNotFoundError:
            print(
                f"Error: Gene expression file not found at {self.gene_expression_file}"
            )
            # Depending on the algorithm, this might be critical or optional
            # For now, we'll allow it to proceed, subclasses can raise errors if needed
            return None
        except Exception as e:
            print(
                f"Error loading gene expression data from {self.gene_expression_file}: {e}"
            )
            return None

    # Method to load other data types can be added here or in subclasses as needed
    def _load_generic_csv(self, file_key, file_path, index_col=None):
        """Generic function to load data from a CSV file."""
        if not file_path:
            print(f"Warning: No file path provided for {file_key}.")
            return None
        try:
            df = pd.read_csv(file_path, index_col=index_col)
            print(f"Data loaded successfully from {file_path} ({file_key}).")
            return df
        except FileNotFoundError:
            print(f"Error: File not found at {file_path} ({file_key}).")
            return None  # Allow proceeding, subclasses handle criticality
        except Exception as e:
            print(f"Error loading data from {file_path} ({file_key}): {e}")
            return None

    @abstractmethod
    def calculate(self, **kwargs):
        """
        Abstract method to calculate the specific centrality score.
        Must be implemented by subclasses.

        Args:
            **kwargs: Algorithm-specific parameters.

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
        # This method can be identical to the one in ClassicalBase
        # Or potentially customized if modern algorithms need different handling
        try:
            df_gold = pd.read_csv(gold_standard_file)
            essential_proteins = set()
            for col in df_gold.columns:
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
        # This method can also be identical to the one in ClassicalBase
        if self.results is None:
            print("Error: No results to compare. Please run calculate() first.")
            return -1

        try:
            essential_proteins = self._load_gold_standard(gold_standard_file)
            if not essential_proteins:
                return 0

            top_n_predicted = self.results["Protein"].head(n).tolist()
            correct_count = len(set(top_n_predicted) & essential_proteins)
            print(f"Comparison: Top {n} predicted proteins vs Gold Standard.")
            print(f"Number of correctly predicted essential proteins: {correct_count}")
            return correct_count

        except Exception as e:
            print(f"Error during comparison: {e}")
            return -1
