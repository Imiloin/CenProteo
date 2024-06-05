import networkx as nx
import pandas as pd


class classical_algorithms:
    def __init__(self, ppi_file):
        """
        Initialize the class with a file path to a protein-protein interaction data CSV.
        This method reads the CSV file and constructs a graph using NetworkX.

        Args:
        ppi_file (str): Path to the CSV file containing protein interaction data.
        essential_protein_file(str): Path to the CSV file containing essential protein data
        """
        self.ppi_file = ppi_file
        df = pd.read_csv(ppi_file)
        G = nx.Graph()
        edges = df.apply(
            lambda row: (row["Protein A"], row["Protein B"]), axis=1
        ).tolist()
        G.add_edges_from(edges)
        self.G = G

    # find essential protein by computing degree centrality
    def DC(self):
        DC = nx.degree_centrality(self.G)
        sorted_DC = sorted(DC.items(), key=lambda x: x[1], reverse=True)
        return sorted_DC

    # find essential protein by computing betweeness centrality
    def BC(self):
        BC = nx.betweenness_centrality(self.G)
        sorted_BC = sorted(BC.items(), key=lambda x: x[1], reverse=True)
        return sorted_BC

    # find essential protein by computing neighbor centrality
    def NC(self):
        """
        Calculate the neighborhood centrality for each protein in the graph.
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
        """
        NC = {}
        for u in self.G.nodes():
            d_u = self.G.degree(u)
            neighbors_u = list(self.G.neighbors(u))
            NC_u = 0
            for v in neighbors_u:
                common_neighbors = list(nx.common_neighbors(self.G, u, v))
                z_uv = len(common_neighbors)
                d_v = self.G.degree(v)
                if min(d_u - 1, d_v - 1) > 0:
                    ECC_uv = z_uv / min(d_u - 1, d_v - 1)
                    NC_u += ECC_uv
            NC[u] = NC_u
        sorted_NC = sorted(NC.items(), key=lambda x: x[1], reverse=True)
        return sorted_NC

    # find essential protein by computing closeness centrality
    def cCC(self):
        CC = nx.closeness_centrality(self.G)
        sorted_CC = sorted(CC.items(), key=lambda x: x[1], reverse=True)
        return sorted_CC

    # find essential protein by computing eigenvector centrality
    def EC(self):
        EC = nx.eigenvector_centrality(self.G)
        sorted_EC = sorted(EC.items(), key=lambda x: x[1], reverse=True)
        return sorted_EC

    # find essential protein by computing information centrality(current flow closeness centrality)
    def IC(self):
        largest_cc = max(
            nx.connected_components(self.G), key=len
        )  # compute the max connected components
        subgraph = self.G.subgraph(
            largest_cc
        )  # construct the max connected subgraph and compute its information_centrality
        IC = nx.current_flow_betweenness_centrality(subgraph)
        sorted_IC = sorted(IC.items(), key=lambda x: x[1], reverse=True)
        return sorted_IC

    # find essential protein by computing subgraph centrality
    def SC(self):
        SC = nx.subgraph_centrality(self.G)
        sorted_SC = sorted(SC.items(), key=lambda x: x[1], reverse=True)
        return sorted_SC

    def export_result_to_csv(self, sorted_result, file_name):
        """
        Export the result data to a CSV file.

        Args:
        result (list of tuples): The result data to be exported, typically the output of a centrality calculation.
        file_name (str): Name of the file to save the results to.
        """
        result_df = pd.DataFrame(sorted_result, columns=["Protein", "Centrality Score"])
        result_df.to_csv(file_name, index=False)

    def _get_essential_protein(self, df):
        essential_pro = []
        for _, row in df.iterrows():
            pro = row.iloc[1]
            essential_pro.append(pro)
        return essential_pro

    def first_n_comparison(self, n, result, real_essential_protein_file):
        """
        Compare the first n elements of the result list and real essential protein list.
        Args:
            result (list): The list of results.
            n: The first n proteins chosen to be compared.
            real_essential_protein_file: The real essential protein file path.

        """
        df_essential = pd.read_csv(real_essential_protein_file)
        self.real_essential_protein_list = self._get_essential_protein(df_essential)

        count = 0

        for protein_tuple in result[:n]:
            protein_name, score = protein_tuple
            if protein_name in self.real_essential_protein_list:
                count = count + 1
        print(
            f"There're {count} essential proteins in the top {n} predicted by algorism."
        )
        return count
