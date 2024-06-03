<div align="center">
    <h1>
    CenProteo: Finding the Essential Proteins in a Protein Interaction Network
    </h1>
    <p>
    Project of BIO2502 Programming Languages for Bioinformatics, 2024 Spring, SJTU
    <br />
    <a href="https://github.com/xywawawa"><strong>xywawawa</strong></a>
    &nbsp;
    <a href="https://github.com/Cannizzaro-reaction"><strong>Cannizzaro-reaction</strong></a>
    &nbsp;
    <a href="https://github.com/Imiloin"><strong>Imiloin</strong></a>
    &nbsp;
    </p>
    <p>
    <a href="https://github.com/Imiloin/CenProteo"><img alt="Github Repository" src="https://img.shields.io/badge/Github-Repository-blue?logo=github&logoColor=blue"></a>
    <a href="https://github.com/Imiloin/CenProteo?tab=MIT-1-ov-file"><img alt="mit" src="https://img.shields.io/badge/License-MIT-red.svg"></a>
    </p>
</div>

在过去的几十年中，对于单一蛋白质的性质及功能方面的研究取得了很大进展。但是，蛋白质在生物体内很少单独发挥作用，因此了解蛋白质之间的相互作用对于揭示复杂分子机制至关重要。近年来，酵母双杂交系统（Yeast Two-Hybrid, Y2H），交叉链接质谱法（Cross-linking Mass Spectrometry, XL-MS）等高通量实验技术快速发展，使得越来越多蛋白质之间的相互作用被研究和发表，也积累了大量的相关实验数据，由此构建出蛋白质相互作用网络（PPIN） 。在 PPIN 中，关键蛋白具有特定的拓扑位置和功能角色，对维持网络的稳定性和功能具有重要影响。为了从 PPIN 中发现关键蛋白，出现了一系列如度中心性（Degree Centrality），介数中心性（Betweenness Centrality），聚类系数（Clustering Coefficient）等传统算法。

本项目构建了包 `CenProteo` ，实现了几种计算蛋白质网络中蛋白质的中心性，并进行排序从而寻找关键蛋白质的算法。



## Data Source & Preprocessing

pass



## Algorithms

蛋白质网络通常表示为一个无向图 $G=(V, E)$，节点 $u\in V$ 表示一个蛋白质，边 $(u,v) \in E$ 表示两个蛋白质之间的相互作用。我们用 $N$ 表示图中节点总数， $A$ 表示图的邻接矩阵。

根据使用的数据类型，算法可以大致分为以下几类：

#### 传统算法

仅使用网络拓扑数据（`CenProteo` 中实现的 `classical algortihms`）计算蛋白质的中心性：

+ DC（degree centrality）度中心性：一个节点 $u$ 的度中心性 $DC(u)$ 是其连接的边数。

    $$DC(u) = \sum_{v} a_{u,v}$$

+ BC（Betweenness Centrality）介数中心性：一个节点 $u$ 的介数中心性 $BC(u)$ 定义为通过节点 $u$ 的最短路径的平均比例。

    $$BC(u) = \sum_{s} \sum_{t} \frac{\rho(s, u, t)}{\rho(s, t)}, \quad s \neq t \neq u$$

    $\rho(s, t)$ 指的是 $s$ 和 $t$ 之间的最短路径数目， $\rho(s, u, t)$ 指的是 $s$ 和 $t$ 之间的最短路径中经过 $u$ 的数目。

+ EC（Eigenvector Centrality）特征向量中心性：一个节点 $u$ 的特征向量中心性 $EC(u)$ 定义为 $A$ 的主特征向量的第 $u$ 分量。

    $$EC(u) = \alpha_{\max}(u)$$

    $\alpha_{\max}$ 指的是 $A$ 的最大值对应的特征向量， $\alpha_{\max}(u)$ 指的是 $\alpha_{\max}$ 的第 $u$ 个分量。

+ SC（Subgraph Centrality）子图中心性：一个节点 $u$ 的子图中心性 $SC(u)$ 衡量的是节点 $u$ 参与的整个网络中子图的数量。

    $$SC(u) = \sum_{l=0}^{\infty} \frac{\mu_{l}(u)}{l!}$$

    $\mu_{l}(u)$ 指的是开始并结束于节点 $u$ 且长度为 $l$ 的环路数目。

+ IC（Information Centrality）信息中心性：一个节点 $u$ 的信息中心性 $IC(u)$ 衡量的是以节点 $u$ 结束的路径长度的调和平均值。

    $$IC(u) = \left[\frac{1}{N} \sum_{v} \frac{1}{I_{uv}}\right]^{-1},I_{uv} = (c_{uu} + c_{vv} - c_{uv})^{-1},C = (c_{uv}) = [D - A + J]^{-1}$$

    $D$ 为每个节点度的对角矩阵， $C$ 是改进的邻接矩阵， $J$ 是所有元素都为 1 的矩阵。

    在 `CenProteo` 中，为简化计算，信息中心性通过计算 `curent flow centrality` 来近似。

+ CC（Closeness Centrality）接近中心性：一个节点 $u$的接近中心性 $CC(u)$ 是从节点 $u$ 到网络中所有其他节点的图理论距离之和的倒数。

    $$CC(u) = \frac{N - 1}{\sum_{v} d(u, v)}$$

    $d(u,v)$ 指的是从节点 $u$ 到结点 $n$ 的距离。

+ NC（Neighbor Centrality）邻居中心性：节点 $u$ 的邻域中心性 $NC(u)$ 定义为节点 $u$ 邻居之间的边聚类系数（Edge Clustering Coefficient, ECC）之和。

    $$NC(u) = \sum_{v \in N_u} ECC(u, v),$$ 
    $$ECC(u, v) = \frac{z_{u, v}}{\min(d_u - 1, d_v - 1)},$$
    $$z_{u, v} = \sum_{w} A_{uw} A_{vw}.$$

    边聚类系数 $ECC(u, v)$ 表示节点 $u$ 和节点 $v$ 之间的共同邻居数 $z_{u, v}$ 与两者度数的最小值之比， $A_{uw}$ 和 $A_{vw}$ 分别表示节点 $u$ 和 $v$ 是否与节点 $w$ 相连。

   

#### 现代算法

使用网络拓扑数据和基因表达量等生物数据，主要实现以下几种方法：

+ TGSO algorithm：

  reference: Li S, Zhang Z, Li X, *et al*. An iteration model for identifying essential proteins by combining comprehensive PPI network with biological information. *BMC bioinformatics*, 2021, 22: 1-25.https://doi.org/10.1186/s12859-021-04300-7

+ JDC algorithm：
  
  **img here**
  reference: Zhong, J., Tang, C., Peng, W. *et al.* A novel essential protein identification method based on PPI networks and gene expression data. *BMC Bioinformatics* 22, 248 (2021). https://doi.org/10.1186/s12859-021-04175-8
  
+ TEO algorithm：
  
  reference: Zhang W, Xu J, Li Y, *et al*. Detecting essential proteins based on network topology, gene expression data, and gene ontology information. *IEEE/ACM transactions on computational biology and bioinformatics*, 2016, 15(1): 109-116.https://ieeexplore.ieee.org/document/7586077



## Installation

#### Clone this repo

```bash
git clone https://github.com/Imiloin/CenProteo.git
cd CenProteo
```

#### Setup

```bash
pip install -e .
```

## Usage

Pass
