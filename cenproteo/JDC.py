import pandas as pd
import networkx as nx
import argparse

"""
The JDC method is based on the PPI network data and gene expression data. 
The JDC method offers a dynamic threshold method to binarize gene expression data. 
After that, it combines the degree centrality and Jaccard similarity index to calculate the JDC score for each protein in the PPI network. 

reference:
Zhong J, Tang C, Peng W, et al. A novel essential protein identification method based on PPI networks and gene expression data[J]. BMC bioinformatics, 2021, 22(1): 248.
"""

class JDC:
    def __init__(self, method, ppi_file, gene_expression_file, n, saved_file_path = None, **kargs):
        # Load PPI network data
        self.ppi_file = ppi_file
        try:
            data1 = pd.read_csv(ppi_file)
        except FileNotFoundError:
            raise Exception(f"File {ppi_file} not found")
        
        G = nx.Graph()
        edges = data1.apply(lambda row: (row['Protein A'], row['Protein B']), axis=1).tolist()
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
        self.n = n
        self.saved_file_path = saved_file_path
        self.method = method

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
                'binary_data': binary_data,
                'mean': mean,
                'std': std
            }
        return gene_expression_dict
    
    def Jaccard_similarity_index(self):
        if self.jaccard is not None:
            return self.jaccard
        
        jaccard = {}
        for u, v in self.G.edges():
            try:
                gene_expression_u = set(self.gene_expression_dict[u]['binary_data'])
                gene_expression_v = set(self.gene_expression_dict[v]['binary_data'])

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
        sorted_jdc = sorted(jdc_dict.items(), key=lambda x: x[1], reverse=True)[:self.n]
        self.sorted_jdc = sorted_jdc
        return sorted_jdc
    
    def export_result_to_csv(self):
        """
        Export the result data to a CSV file.
        
        Args:
        result (list of tuples): The result data to be exported, typically the output of a centrality calculation.
        file_name (str): Name of the file to save the results to.
        """
        if self.saved_file_path == None:
            raise ValueError("saved_file_path can't be None")
        result_df = pd.DataFrame(self.sorted_jdc, columns=['Protein', 'JDC Centrality Score'])
        result_df.to_csv(self.saved_file_path, index=False)

    def start(self):
        if self.method == "calculate_jdc":
            self.calculate_jdc()
        if self.method == "export_result_to_csv":
            self.export_result_to_csv()
        


# Example usage
#jdc_instance = JDC('/Users/xiaoyao/Documents/ppi/data/DIP_data_with_combined_scores.csv', '/Users/xiaoyao/Documents/ppi/filtered_GE_matrix.csv', 100)
#top_nodes = jdc_instance.calculate_jdc

def parse_command_line_args():
    parser = argparse.ArgumentParser(description='JDC Algorithm')
    parser.add_argument("--method", required=True, type=str, help="please select the method")
    parser.add_argument('-ppi', '--ppi_file', required=True, type=str,help='protein protein network file path for algorithm JDC')
    parser.add_argument('-ge', '--gene_expression_file',  required=True, type=str, help='gene exprrssion file path for algorithm JDC')
    parser.add_argument('-n', '--n', type=int, required=True, help='the number of essential proteins output')
    parser.add_argument('-s', '--saved_file_path', type=str, default=None, help='the file path to save the result')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_command_line_args()
    jdc = JDC(**vars(args))
    jdc.start()
