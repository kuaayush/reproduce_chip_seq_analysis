# Codes I used for reproducing ChIP-Seq data analysis from a published research paper

I reproduced figures from a published research paper by following a tutorial by Dr. Ming “Tommy” Tang. Here, I reproduced the figures from a study titled “Genome-wide association between YAP/TAZ/TEAD and AP-1 at enhancers drives oncogenic growth”, conducted by Zanconato et al. in Italy back in 2015, which unravels how different protein groups team up to control gene activity in cancer cells.

You can also read a blog I wrote about it here: [My blog on Medium](https://medium.com/@iamaayushojha/reproducing-research-papers-and-surviving-chip-seq-4328ccfaff06?source=friends_link&sk=7ce4864e86c33af056fdea6449b20ba4)


_You can read the brief biology behind the data analysis here!_

# Completed a tutorial to reproduce figures from a published research paper

Reproducing the figures of any research paper is beneficial because: 
  - It helps you to be familiar with data types involved in different kinds of research
  - You will get to learn some biological concepts in depth.
  - You will learn what figures are better to use for what kind of data. 

Here, I followed the tutorial of Tommy Tang to reproduce figures from the paper titled “Genome-wide association between YAP/TAZ/TEAD and AP-1 at enhancers drives oncogenic growth”. This study was conducted in Italy, by Zanconato et. al, published in 2015. The study is focused on identifying how two groups of proteins work together to control the gene activity in cancer cells. 

**Hippo Pathway**

Hippo Pathway is a signalling pathway that controls cell proliferation and apoptosis to regulate organ size. It is a very important pathway involved in regulating tissue generation, apicobasal cell polarity, mechanotransduction, tissue specific stem cells, and  controlling cell contact inhibition. There are three main components if the Hippo Pathway:
  - **Hippo**: Hippo is a protein kinase
  - **MST and Lats**: These are highly conserved kinase cascades
  - **YAP/TAZ**: These are downstream transcriptional co-activators and is pathway's effector protein.
    - YAP means Yes associated protein
    - TAZ means Transcriptional co-activator with PDZ-binding motif

**General mechanism of Hippo Pathway**

When Hippo Pathway is turned off, YAP/TAZ is not phosphorylated,and hence can move from cytoplasm to the nucleus. Once it gets inside the nucleus, it binds with TEAD1-4 at the promoter region. This binding signals that the gene expression is “on. 
But, when Hippo Pathway is said to be “on”, MST1/2 and MAP4Ks get phosphorylated and are activated. Their  activation phosphorylates the downstream protein LATS1/2 binding with MOB1. The phosphorylated LATS1/2-MOB1 complex phosphorylates YAP/TAZ protein. Once the YAP/TAZ protein is phosphorylated, it cannot  move inside the nucleus (cytoplasmic retention is seen) and degrades in the cytoplasm.  As the YAP/TAZ did not bind with TEAD1-4 at the promoter region, the gene expression is off. 

**Role of Hippo Pathway in oncogenic growth**

YAP/TAZ-induced oncogenic growth is strongly enhanced by gain-of-AP1 and severely blunted by its loss. 

**Histone Modification**

Typical patterns of histone methylations exhibited at promoters, insulators, enhancers, and transcribed regions are identified. The mono methylation of H3K27, H3K9, H4K20, H3K79, and H2BK5 are all linked to gene activation, whereas trimethylation of H3K27, H3K9, and H3K79 are linked to repression. H2A.Z associates with functional regulatory elements, and CTCF marks boundaries of histone methylation domains. Chromosome banding patterns are correlated with unique patterns of histone modifications. Chromosome breakpoints detected in T cell cancers frequently reside in chromatin regions associated with H3K4 methylations. Our data provide new insights into the function of histone methylation and chromatin organization in genome function.

**Histone Modifications seen in the Hippo Pathway**

_**Key points about H3K4me1 and the Hippo pathway:**_

In the Hippo signaling pathway, H3K4me1 (histone H3 lysine 4 monomethylation) is associated with the transcriptional activation of target genes by the pathway's effector protein, YAP/TAZ, indicating a role in regulating chromatin structure and promoting gene expression when the Hippo pathway is inhibited; essentially, increased H3K4me1 levels at YAP/TAZ target genes can lead to enhanced transcription of those genes, often linked to cell proliferation and growth. [1, 2, 3, 4, 5]
  - Histone modification: H3K4me1 is a histone modification where a single methyl group is added to lysine 4 on histone H3, typically associated with active transcription. [1, 2, 3]
  - YAP/TAZ interaction: Studies have shown that YAP/TAZ, the central transcriptional coactivators of the Hippo pathway, can recruit histone methyltransferases that add H3K4me1 marks to their target gene promoters, thereby promoting transcription. [2, 3, 4]
  - Regulation of cell growth: By influencing the expression of genes involved in cell proliferation and differentiation, the interplay between H3K4me1 and the Hippo pathway plays a crucial role in controlling organ size and tissue development. [1, 2, 5, 6]
  - Cancer implications: Aberrant activation of the Hippo pathway, often due to decreased YAP/TAZ activity, can lead to tumorigenesis, and dysregulation of H3K4me1 at YAP/TAZ target genes is implicated in cancer progression. [1, 4, 5] 

**Sources:**

  - [1] https://academic.oup.com/nar/article/47/5/2349/5289487
  - [2] https://www.science.org/doi/10.1126/sciadv.abe8159
  - [3] https://pmc.ncbi.nlm.nih.gov/articles/PMC6562424/
  - [4] https://pmc.ncbi.nlm.nih.gov/articles/PMC4644679/
  - [5] https://pubmed.ncbi.nlm.nih.gov/34173335/
  - [6]https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/hippo-signaling-pathway

_**Key points about H3K4me3 and the Hippo pathway:**_

In the Hippo pathway, H3K4me3 acts as a marker of active transcription, indicating that when the Hippo pathway is inactive, the target genes associated with cell growth and proliferation are actively expressed due to the presence of this histone modification at their transcription start sites, essentially promoting gene expression by making chromatin more accessible to the transcription machinery. [1, 2, 3, 4, 5] 
  - Active chromatin mark: H3K4me3 (trimethylation of histone H3 at lysine 4) is considered a key "active" histone mark, signifying regions of DNA that are actively transcribed. [1, 5, 6]
  - YAP/TAZ as effectors: The Hippo pathway primarily regulates cell growth through the activity of its downstream effectors, YAP and TAZ, which act as transcriptional coactivators when the pathway is inhibited. [2, 3, 4]
  - Regulation of gene expression: When the Hippo pathway is active, YAP/TAZ are inhibited, preventing them from binding to DNA and activating target genes; conversely, when the Hippo pathway is inactive, YAP/TAZ are free to bind to DNA and promote transcription, often marked by increased H3K4me3 at the target gene promoters. [1, 3, 4]
  - Mechanism: The Hippo pathway can influence the level of H3K4me3 by regulating the activity of histone methyltransferases, which add the methyl group to histone H3 at lysine 4, thereby facilitating transcription. [1, 7, 8] 


**Sources:**

  - [1] https://pmc.ncbi.nlm.nih.gov/articles/PMC10240540/
  - [2] https://www.sciencedirect.com/science/article/pii/S2211124715002685
  - [3] https://academic.oup.com/nar/article-pdf/47/5/2349/28042012/gky1317.pdf
  - [4] https://academic.oup.com/nar/article/51/9/4266/7067949
  - [5] https://pmc.ncbi.nlm.nih.gov/articles/PMC4137894/
  - [6] https://pmc.ncbi.nlm.nih.gov/articles/PMC8714964/
  - [7] https://en.wikipedia.org/wiki/H3K4me3
  - [8] https://www.sciencedirect.com/science/article/pii/S2211124714004859
    

_**Key Players in the paper I reproduced:**_
  - _**YAP/TAZ**_: Proteins that are important for regulating organ size and are linked to cancer when they malfunction.
  - _**TEAD Factors**_: Proteins that help YAP/TAZ attach to DNA.
  - _**AP-1 (Activator Protein-1)**_: A group of proteins (made of JUN and FOS) that also binds DNA and helps control gene activity.

Also, 
  - promoters are YAP/TAZ/TEAD4 peaks overlapping with H3K4me3 peaks
  - active enhancers are YAP/TAZ/TEAD4 peaks overlapping with H3K4me1 and H3K27ac peaks
  - inactive enhancers are YAP/TAZ/TEAD4 peaks overlapping with H3K4me1 but not H3K37ac peaks

**How These Proteins Work Together:**

This study used a technique called ChIP-seq on breast cancer cells to see where these proteins attach to the DNA. They found that YAP/TAZ don’t act alone. Instead, they join forces with TEAD and AP-1. The DNA regions where they bind have specific signals (motifs) recognized by both TEAD and AP-1. When YAP/TAZ, TEAD, and AP-1 come together, they form a complex that turns on genes needed for a cell to progress through its cycle, particularly genes that control when a cell starts copying its DNA (S-phase) and when it divides (mitosis). Interestingly, the activation of these genes mostly happens through distant regulatory regions called enhancers. These enhancers loop around to contact the gene’s promoter (the region where gene activation starts), highlighting the importance of the 3D structure of DNA in gene regulation.

Ultimately, the study shows that increasing AP-1 levels can boost the cancer-promoting effects of YAP/TAZ, while losing AP-1 diminishes these effects. In animal models (mice), removing YAP/TAZ was able to prevent AP-1-driven skin tumor formation, demonstrating that both are crucial for tumor growth.

**What I learned through this process**

  - I learned about the ChIP-seq technique, its significance, and what kind of biological questions can be answered through chip seq experiments.
  - I became familiar with ChIP-seq data. I examined different file formats involved in chip- seq, from .bedGraph files to .BigWig files, and learned their individual importance, and the information they carried.
  - I had almost forgotten to code in R lately. So, analyzing and visualizing the ChIP-seq data using different R packages sharpened my R skills, and I also got to explore and work with packages that I was unfamiliar with ( diplyr, tidyr, EnhancedHeatMap).
  - I used the IGV genome browser for the first time, and realized it's such an important tool as it helps us to visualize the genomic regions that are regulatory, and have important functions. Before this, I have only visualized whole genomes, but we all know whole genomes are not used for any activity. Only certain genomic regions are used (called exons), and these regions (called genes altogether) also don't work alone. They need regulatory proteins to carry out their function. 

**Reference paper:**

1. Zanconato F, Forcato M, Battilana G, Azzolin L, Quaranta E, Bodega B, Rosato A, Bicciato S, Cordenonsi M, Piccolo S. Genome-wide association between YAP/TAZ/TEAD and AP-1 at enhancers drives oncogenic growth. Nat Cell Biol. 2015 Sep;17(9):1218-27. doi: 10.1038/ncb3216. Epub 2015 Aug 10. PMID: 26258633; PMCID: PMC6186417.

