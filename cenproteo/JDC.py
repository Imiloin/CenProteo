import pandas as pd
import networkx as nx

"""
The JDC method is based on the PPI network data and gene expression data. 
The JDC method offers a dynamic threshold method to binarize gene expression data. 
After that, it combines the degree centrality and Jaccard similarity index to calculate the JDC score for each protein in the PPI network. 

reference:
Zhong J, Tang C, Peng W, et al. A novel essential protein identification method based on PPI networks and gene expression data[J]. BMC bioinformatics, 2021, 22(1): 248.
"""


class JDC:
    def __init__(self, ppi_file, gene_expression_file):
        # Load PPI network data
        self.ppi_file = ppi_file
        try:
            data1 = pd.read_csv(ppi_file)
        except FileNotFoundError:
            raise Exception(f"File {ppi_file} not found")

        G = nx.Graph()
        edges = data1.apply(
            lambda row: (row["Protein A"], row["Protein B"]), axis=1
        ).tolist()
        G.add_edges_from(edges)
        self.G = G

        # Load gene expression data
        self.gene_expression_file = gene_expression_file
        try:
            data2 = pd.read_csv(gene_expression_file, index_col=0)
        except FileNotFoundError:
            raise Exception(f"File {gene_expression_file} not found")

        self.gene_expression_dict = self.binariztion_gene_expression(data2)
        self.ecc = None
        self.jaccard = None

        self.sorted_score = self.calculate_jdc()

    def _get_essential_protein(self, df):
        essential_pro = []
        for _, row in df.iterrows():
            pro = row.iloc[1]
            essential_pro.append(pro)
        return essential_pro

    def edge_clustering_coefficient(self):
        if self.ecc is not None:
            return self.ecc

        ecc = {}
        for u, v in self.G.edges():
            common_neighbors = list(nx.common_neighbors(self.G, u, v))
            actual_triangles = len(common_neighbors)
            possible_triangles = min(self.G.degree[u], self.G.degree[v]) - 1
            if possible_triangles > 0:
                ecc[(u, v)] = actual_triangles / possible_triangles
            else:
                ecc[(u, v)] = 0
        self.ecc = ecc
        return ecc

    def binariztion_gene_expression(self, df):
        gene_expression_dict = {}
        for index, row in df.iterrows():
            data = row[:-2].tolist()
            mean = row.iloc[-2]
            std = row.iloc[-1]
            volatility = 1 / (1 + std)
            threshold = mean + 2 * std + volatility
            binary_data = [1 if value > threshold else 0 for value in data]
            gene_expression_dict[index] = {
                "binary_data": binary_data,
                "mean": mean,
                "std": std,
            }
        return gene_expression_dict

    def Jaccard_similarity_index(self):
        if self.jaccard is not None:
            return self.jaccard

        jaccard = {}
        for u, v in self.G.edges():
            try:
                gene_expression_u = set(self.gene_expression_dict[u]["binary_data"])
                gene_expression_v = set(self.gene_expression_dict[v]["binary_data"])

                intersection = len(gene_expression_u.intersection(gene_expression_v))
                union = len(gene_expression_u.union(gene_expression_v))

                jaccard_similarity_index = intersection / union if union != 0 else 0
                jaccard[(u, v)] = jaccard_similarity_index
            except KeyError:
                jaccard[(u, v)] = 0
        self.jaccard = jaccard
        return jaccard

    def calculate_jdc(self):
        self.edge_clustering_coefficient()
        self.Jaccard_similarity_index()

        jdc_dict = {}
        for node in self.G.nodes():
            jdc_value = 0
            for neighbor in self.G.neighbors(node):
                try:
                    jaccard_value = self.jaccard[(node, neighbor)]
                    ecc_value = self.ecc[(node, neighbor)]
                    jdc_value += jaccard_value * ecc_value
                except KeyError:
                    continue

            jdc_dict[node] = jdc_value
        sorted_jdc = sorted(jdc_dict.items(), key=lambda x: x[1], reverse=True)
        self.sorted_jdc = sorted_jdc
        return sorted_jdc

    def export_result_to_csv(self, save_path):
        """
        Export the result data to a CSV file.

        Args:
        result (list of tuples): The result data to be exported, typically the output of a centrality calculation.
        file_name (str): Name of the file to save the results to.
        """
        result_df = pd.DataFrame(
            self.sorted_jdc, columns=["Protein", "JDC Centrality Score"]
        )
        result_df.to_csv(save_path, index=False)

    def first_n_comparison(self, n, real_essential_protein_file):
        df_essential = pd.read_csv(real_essential_protein_file)
        self.essential_protein_list = self._get_essential_protein(df_essential)
        count = 0

        for protein_tuple in self.sorted_jdc[:n]:
            protein_name, score = protein_tuple
            if protein_name in self.essential_protein_list:
                count = count + 1
        print(
            f"There're {count} essential proteins in the top {n} predicted by algorism."
        )
        return count
