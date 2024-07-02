import networkx as nx
import pandas as pd
import csv


class TEO:
    def __init__(self, ppi_file, gene_expression_file):
        # load ppi network data
        self.ppi_file = ppi_file
        df_ppi = pd.read_csv(ppi_file)
        G = nx.Graph()
        edges = df_ppi.apply(
            lambda row: (row["Protein A"], row["Protein B"]), axis=1
        ).tolist()
        G.add_edges_from(edges)
        self.G = G

        # load gene expression data
        self.gene_expression_filte = gene_expression_file
        df_expression = pd.read_csv(gene_expression_file, index_col=0)
        self.gene_expression_dict = self._create_expression_dict(df_expression)

        # load GO similarity data
        df_GO = pd.read_csv(ppi_file, index_col=[0, 1])
        self.GO_similarity_dict = self._create_GO_dict(df_GO)

        # the result dict using three kinds of GO term
        self.TEO_BP = self.TEO("BP")
        self.TEO_MF = self.TEO("MF")
        self.TEO_CC = self.TEO("tCC")

    def _get_essential_protein(self, df):
        essential_pro = (
            []
        )  # a list used to record all the essential proteins appear in the file
        for _, row in df.iterrows():
            pro = row.iloc[1]  # the name of the essential proteins
            essential_pro.append(pro)
        return essential_pro

    def _create_GO_dict(self, df):
        """
        Construct a dict according to the GO similarity data,
        in which the GO similarity under BP term can be accessed by 'BP',
        the GO similarity under MF term can be accessed by 'MF',
        and the GO similarity under tCC term can be accessed by 'tCC'
        """
        GO_dict = {}
        for index, row in df.iterrows():
            GO_BP = row.iloc[
                0
            ]  # read the GO term one by one and add them to the dictionary, the key is the protein name
            GO_MF = row.iloc[1]
            GO_CC = row.iloc[2]
            GO_dict[index] = {
                "BP": GO_BP,
                "MF": GO_MF,
                "tCC": GO_CC,
            }
            GO_dict[(index[1], index[0])] = GO_dict[index]
        return GO_dict

    def _create_expression_dict(self, df):
        """
        Construct a dict according to the gene expression data,
        in which the detailed data can be accessed by 'data',
        the average of the expression data by 'mean',
        and the standard deviation by 'std'
        """
        gene_expression_dict = {}
        for index, row in df.iterrows():
            data = row[:-2].tolist()
            mean = row.iloc[-2]
            std = row.iloc[-1]
            gene_expression_dict[index] = {"data": data, "mean": mean, "std": std}
        return gene_expression_dict

    def edge_clustering_coefficient(self, u, v):
        d_u = self.G.degree(u)  # the degree of a protein
        d_v = self.G.degree(v)
        d_min = min((d_u - 1), (d_v - 1))

        neighbors_u = set(
            self.G.neighbors(u)
        )  # the neighbors of a protein in a network
        neighbors_v = set(self.G.neighbors(v))
        triangle = 0
        for neighbor in neighbors_u:  # find the number of common neighbors
            if neighbor in neighbors_v:
                triangle += 1

        if d_min == 0 or triangle == 0:
            return 0
        return (triangle**3) / (d_min)

    def pearson_correlation_coefficient(self, u, v):
        # calculate pcc value of each protein pair
        if (u not in self.gene_expression_dict) or (v not in self.gene_expression_dict):
            return "NA"
        n = len(self.gene_expression_dict[u]["data"])
        mean_u, std_u = (
            self.gene_expression_dict[u]["mean"],
            self.gene_expression_dict[u]["std"],
        )
        mean_v, std_v = (
            self.gene_expression_dict[v]["mean"],
            self.gene_expression_dict[v]["std"],
        )
        list_u = self.gene_expression_dict[u]["data"]
        list_v = self.gene_expression_dict[v]["data"]
        pcc = 0
        for i in range(n):
            x_i = list_u[i]
            y_i = list_v[
                i
            ]  # x and y refers to the elements in the expression list created above
            add = ((x_i - mean_u) / std_u) * ((y_i - mean_v) / std_v)
            pcc += add
        pcc /= n - 1
        return abs(pcc)

    def GO_similarity(self, u, v, GO_term):
        if GO_term in [
            "BP",
            "MF",
            "tCC",
        ]:  # choose the GO value according to the GO term we choose
            go = self.GO_similarity_dict[(u, v)][GO_term]
            return go
        return "GO term error"

    def TEO(self, GO_term):
        # calculate the final TEO score of each protein in the PPIN
        if GO_term not in ["BP", "MF", "tCC"]:
            return "GO term error"
        TEO_score = {}  # record the final score of each protein
        for protein_a in self.G.nodes():
            a_TEO = 0
            for protein_b in self.G.neighbors(protein_a):
                if self.pearson_correlation_coefficient(protein_a, protein_b) == "NA":
                    continue
                ECC = self.edge_clustering_coefficient(protein_a, protein_b)
                GO_sim = self.GO_similarity(protein_a, protein_b, GO_term)
                PCC = self.pearson_correlation_coefficient(protein_a, protein_b)
                a_TEO += ECC * (
                    GO_sim + PCC
                )  # calculate the TEO score for protein A with every protein it connected to in the network
            TEO_score[protein_a] = a_TEO
        sorted_TEO_score = dict(
            sorted(TEO_score.items(), key=lambda item: item[1], reverse=True)
        )  # sort the result
        return sorted_TEO_score

    def first_n_comparison(self, n, GO_term, real_essential_protein_file):
        """
        Compare the first n results in the final outcome with the standard essential protein file.
        The number given out by the result refers to how many proteins in the top first n is included in the essential protein list.
        """
        df_essential = pd.read_csv(real_essential_protein_file)
        self.essential_protein_list = self._get_essential_protein(df_essential)
        # Evaluate the efficiency of the algorism by counting how many proteins with high teo score (top n) exist in the essential protein list
        count = 0
        if GO_term == "BP":
            top_TEO_score = list(self.TEO_BP.keys())[:n]
        elif GO_term == "MF":
            top_TEO_score = list(self.TEO_MF.keys())[:n]
        elif GO_term == "tCC":
            top_TEO_score = list(self.TEO_CC.keys())[:n]

        for ess_pro in top_TEO_score:
            if ess_pro in self.essential_protein_list:
                count += 1
        print(
            f"There're {count} essential proteins in the top {n} predicted by TEO_{GO_term} algorism."
        )
        return count

    def export_results_to_csv(self, GO_term, result_path):
        # get the TEO score result in a csv file
        if GO_term == "BP":
            score_list = list(self.TEO_BP.items())
        elif GO_term == "MF":
            score_list = list(self.TEO_MF.items())
        elif GO_term == "tCC":
            score_list = list(self.TEO_CC.items())

        with open(result_path, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Protein", f"TEO_Score_{GO_term}"])
            for key, value in score_list:
                writer.writerow([key, value])
