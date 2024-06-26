import networkx as nx
import pandas as pd
import csv


class TGSO:
    def __init__(
        self,
        ppi_file,
        gene_expression_file,
        subcellular_localization_file,
        i_score_file,
        alpha=0.3,
        max_iter=100,
        tol=1e-5,
    ):
        # load ppi network data
        self.ppi_file = ppi_file
        df_ppi = pd.read_csv(ppi_file)
        # construct a PPI netwrok
        G = nx.Graph()
        edges = df_ppi.apply(
            lambda row: (row["Protein A"], row["Protein B"]), axis=1
        ).tolist()
        G.add_edges_from(edges)
        self.G = G

        # load gene expression file
        self.gene_expression_file = gene_expression_file
        df_expression = pd.read_csv(gene_expression_file, index_col=0)
        self.gene_expression_dict = self._create_expression_dict(df_expression)

        # load subcellular localization file
        self.localization_file = subcellular_localization_file
        self.df_localization = pd.read_csv(subcellular_localization_file, index_col=0)

        # Load I_score file
        self.i_score_file = i_score_file
        self.i_score_dict = self._load_i_score(i_score_file)

        self.ADN = self.ADN()
        self.CEN = self.CEN()
        self.colo_sub = self.CLN()

        # the recommended scores are in default
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol

        self.P0 = self._initialize_scores()

        # get the final score for each protein and the iteration time
        self.protein_score, self.iter_time = self.calculate_P()

    def _get_essential_protein(self, df):
        essential_pro = []  # a list to record the essential proteins appear in the file
        for _, row in df.iterrows():
            pro = row.iloc[1]
            essential_pro.append(pro)
        return essential_pro

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

    def _load_i_score(self, file_path):
        """
        The oncology score of each protein in 100 different species.
        the I_score refers to the times each homologous gene appear in different species.
        """
        df_i_score = pd.read_csv(file_path, index_col=0)
        df_i_score = df_i_score.dropna()
        i_score_dict = df_i_score["Oscore"].to_dict()  # Oscore = i_score/total i_score
        return i_score_dict

    # construction of protein Aggregation Degree interactive Netwrok
    def ADN(self):
        ADN = {}
        for protein_a in self.G.nodes:
            NG_protein_a = list(
                self.G.neighbors(protein_a)
            )  # a list for all the neighbors proteins of protein A
            for (
                protein_b
            ) in NG_protein_a:  # traverse all the neighbor proteins of protein A
                if (protein_b, protein_a) in ADN:
                    ADN[(protein_a, protein_b)] = ADN[(protein_b, protein_a)]
                    continue
                NG_protein_b = list(self.G.neighbors(protein_b))
                numerator = 1
                for conj_pro in NG_protein_a:
                    if conj_pro in NG_protein_b:
                        numerator += 1  # calculate the number of common neighbors
                denominator = min(len(NG_protein_a), len(NG_protein_b))
                ADN[(protein_a, protein_b)] = numerator / denominator
        return ADN

    def pearson_correlation_coefficient(self, u, v) -> float:
        # calculate the pcc value of a protein pair
        if (u not in self.gene_expression_dict) or (v not in self.gene_expression_dict):
            return 0.0
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
            y_i = list_v[i]
            add = ((x_i - mean_u) / std_u) * ((y_i - mean_v) / std_v)
            pcc += add
        pcc /= n - 1
        return abs(pcc)

    # construction of protein co-expression interaction Network
    def CEN(self):
        CEN = {}
        for protein_a in self.G.nodes:
            NG_protein_a = list(self.G.neighbors(protein_a))
            for protein_b in NG_protein_a:
                if (protein_b, protein_a) in CEN:
                    CEN[(protein_a, protein_b)] = CEN[(protein_b, protein_a)]
                    continue
                pcc = self.pearson_correlation_coefficient(protein_a, protein_b)
                NG_protein_b = list(self.G.neighbors(protein_b))
                connection = (
                    pcc  # calculate the connection value based on the pcc value
                )
                for conj_pro in NG_protein_a:
                    if (
                        conj_pro in NG_protein_b
                    ):  # get the common neighbor of protein A and protein B
                        pcc_1 = self.pearson_correlation_coefficient(
                            protein_a, conj_pro
                        )
                        pcc_2 = self.pearson_correlation_coefficient(
                            protein_b, conj_pro
                        )
                        connection += pcc_1 * pcc_2
                CEN[(protein_a, protein_b)] = connection
        return CEN

    # construction of protein Co-Localization interaction Network
    def CLN(self):
        # calculate the colo_sub value according to the protein sublocalization information
        pro_localization = (
            {}
        )  # record the subcellular localizations possessed by each protein
        sub_pro = (
            {}
        )  # record the number of proteins appear in each subcellular localization
        go_terms = {
            "Nucleus": "GO:0005634",
            "Cytosol": "GO:0005829",
            "Cytoskeleton": "GO:0005856",
            "Peroxisome": "GO:0005777",
            "Lysosome": "GO:0005764",
            "Endoplasmic Reticulum": "GO:0005783",
            "Golgi Apparatus": "GO:0005794",
            "Plasma Membrane": "GO:0005886",
            "Endosome": "GO:0005768",
            "Extracellular Region": "GO:0005576",
            "Mitochondrion": "GO:0005739",
        }
        for index, row in self.df_localization.iterrows():
            protein_name = index
            go_term = row["GO_term"]
            for compartment, term in go_terms.items():
                if go_term == term:
                    if compartment not in sub_pro:
                        sub_pro[compartment] = []
                    sub_pro[compartment].append(protein_name)

        for compartment, pro_list in sub_pro.items():
            for protein in pro_list:
                if protein not in pro_localization:
                    pro_localization[protein] = []
                pro_localization[protein].append(compartment)

        sub_score = (
            {}
        )  # calculate the sub_score of each compartment in the go_terms dictionary
        total_score = 0
        for compartment in sub_pro:
            score = len(sub_pro[compartment])
            sub_score[compartment] = score
            total_score += score

        for compartment in sub_score:
            sub_score[compartment] /= total_score

        # calculate the S_score of each protein (the sum of all the sub_score this protein appears in)
        pro_S_score = {}
        for protein in pro_localization:
            S_total = 0
            for compartment in pro_localization[protein]:
                S_total += sub_score[compartment]
            pro_S_score[protein] = S_total

        # calculate the co-localization score of each protein pair
        colo_sub = {}
        for protein_a in self.G.nodes:
            for protein_b in self.G.neighbors(protein_a):
                if (protein_b, protein_a) in colo_sub:
                    colo_sub[(protein_a, protein_b)] = colo_sub[(protein_b, protein_a)]
                    continue

                if protein_a not in pro_localization:
                    localization_a = []
                else:
                    localization_a = pro_localization[protein_a]
                if protein_b not in pro_localization:
                    localization_b = []
                else:
                    localization_b = pro_localization[protein_b]
                intersection = list(set(localization_a) & set(localization_b))
                union = list(set(localization_a) | set(localization_b))
                len_intersection = len(
                    intersection
                )  # the number of the common sublocalization of protein A and B
                len_union = len(
                    union
                )  # the number of union sublocalization of protein A and B

                if protein_a not in pro_S_score:
                    S_a = 0
                else:
                    S_a = pro_S_score[protein_a]
                if protein_b not in pro_S_score:
                    S_b = 0
                else:
                    S_b = pro_S_score[protein_b]
                if len_union == 0:
                    colo_sub[(protein_a, protein_b)] = 0
                else:
                    colo_sub[(protein_a, protein_b)] = (
                        (len_intersection / len_union) * (S_a + S_b) / 2
                    )

        return colo_sub

    # construction of a comprehensive interaction between two proteins
    def LSG(self):
        # calculate the LSG value of each protein
        LSG = {}
        for protein_u in self.G.nodes:
            lsg_u = 0
            for protein_v in self.G.neighbors(protein_u):
                if (
                    protein_u,
                    protein_v,
                ) not in self.colo_sub:  # no information about the exact localization of the protein
                    co_sub = 0
                else:
                    co_sub = self.colo_sub[(protein_u, protein_v)]
                lsg_uv = self.ADN[(protein_u, protein_v)] * (
                    co_sub + self.CEN[(protein_u, protein_v)]
                )
                lsg_u += lsg_uv
            LSG[protein_u] = lsg_u
        return LSG

    def _initialize_scores(self):
        # the initial P score should be the P0 score
        P0 = self.i_score_dict.copy()
        return P0

    def calculate_P(self):
        LSG = self.LSG()
        P0 = self.P0.copy()
        LSG_total = sum(LSG.values())

        P = P0
        iter_time = 0
        for _ in range(self.max_iter):  # avoid too much iteration (error)
            P_new = {}
            iter_time += 1
            for protein_i in self.G.nodes:
                p_i = 0

                a = self.alpha
                for protein_j in self.G.nodes:
                    # calculate the PCIN score for each two proteins (not necessarily a interaction pair)
                    if protein_i != protein_j:
                        PCIN = min(LSG[protein_i], LSG[protein_j]) / LSG_total
                    else:
                        PCIN = LSG[protein_i] / LSG_total
                    if protein_j not in P:
                        p2 = 0
                    else:
                        p2 = P[protein_j]
                    if protein_j not in P0:
                        p0 = 0
                    else:
                        p0 = P0[protein_j]
                    p_j = (1 - a) * PCIN * p2 + a * p0
                    p_i += p_j
                P_new[protein_i] = p_i

            # Check for convergence
            diff = max(
                abs(P_new.get(protein, 0) - P.get(protein, 0))
                for protein in self.G.nodes
            )
            E = len(self.G.edges)
            if (diff / E) < self.tol:
                break
            P = P_new  # if not converge, update the P dict

        sorted_protein_score = dict(
            sorted(P.items(), key=lambda item: item[1], reverse=True)
        )
        return sorted_protein_score, iter_time

    def first_n_comparison(self, n, real_essential_protein_file):
        """
        Compare the first n results in the final outcome with the standard essential protein file.
        The number given out by the result refers to how many proteins in the top first n is included in the essential protein list.
        """
        # load essential protein data
        df_essential = pd.read_csv(real_essential_protein_file)
        self.essential_protein_list = self._get_essential_protein(df_essential)
        count = 0
        top_TGSO_score = list(self.protein_score.keys())[:n]
        for ess_pro in top_TGSO_score:
            if ess_pro in self.essential_protein_list:
                count += 1
        print(
            f"There're {count} essential proteins in the top {n} predicted by TGSO algorism. \nThe iteration has been repeated for {self.iter_time} times."
        )
        return count

    def export_results_to_csv(self, file_path):
        # Ensure the 'protein_score' dictionary is not empty
        if not self.protein_score:
            print("No scores to export. Please run the calculation first.")
            return

        # Open the file and write the headers and data
        with open(file_path, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Protein", "Score"])  # Writing header
            for protein, score in sorted(
                self.protein_score.items(), key=lambda item: item[1], reverse=True
            ):
                writer.writerow([protein, score])  # Writing each protein and its score
