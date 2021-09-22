## Step 1: Load Sleuth and accessory libraries 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("rlang")
BiocManager::install("pachterlab/sleuth")
BiocManager::install(c("rhdf5"))
library("sleuth")

# Additional libraries needed for plotting, etc
install.packages(c("gridExtra","cowplot"))
library("gridExtra")
library("cowplot")
BiocManager::install("biomaRt")
library("biomaRt")

## Step 2: Load experimental design and label kallisto outputs with metadata
# We need to provide Sleuth with our sample names:
sample_id <- dir(file.path("~/Downloads/sleuth-main/Analysis"))
sample_id

# We also need to get the file paths to our results files. 
kal_dirs <- file.path("~/Downloads/sleuth-main/Analysis", sample_id)
kal_dirs

# We also need a table that provides more meaningful names for describing our experiment...
s2c <- read.table(file.path("~/Downloads/sleuth-main/kallisto_demo.tsv"),
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  sep = "\t")

# We will add our file paths to the table
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# And take a look at the table we just created
s2c

## Step 3: Load gene names from Ensembl 
# We want the `athaliana_eg_gene` data set so we load that into R. 
plants_mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org" )

# Now we want to get specific attributes from the list of genes we can import from biomart
listAttributes(plants_mart)

# We can choose whichever of these we'd like to use. Let's get transcript ids, gene ids, a description, and gene names.
t2g <- getBM(attributes = c("ensembl_transcript_id", 
                            "ensembl_gene_id", 
                            "description",
                            "external_gene_name"), 
             mart = plants_mart)

# We need to make sure the `ensembl_transcript_id` column is named `target_id`
ttg <- dplyr::rename(t2g, target_id= ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

##Step 4: Prepare data for Sleuth 
# first we need to alter our experimental design so that we consider the full transcriptome sample to be the "control" to compare to...
s2c$genotype_variation_s <- as.factor(s2c$genotype_variation_s)
s2c$genotype_variation_s <- relevel(s2c$genotype_variation_s, ref = "wild type")

# Now we need to tell Sleuth both about the Kallisto results and the gene names (and gene descriptions/metadata) we obtained from biomaRt. 
# The `sleuth_prep` function does this. 
so <- sleuth_prep(s2c,
                  full_model = ~genotype_variation_s,
                  target_mapping = ttg,
                  read_bootstrap_tpm=TRUE,
                  extra_bootstrap_summary = TRUE)

##Step 5: Initial data exploration
### Examine Sleuth PCA
# Next, we should check to see if our samples (and replicates) cluster on a PCA (as should expect) or if there are outliers. 
# When we plot by condition, we'd expect that similar colors group together. 
ggplot2::theme_set(theme_cowplot())
plot_pca(so, color_by = 'genotype_variation_s', text_labels = TRUE)

# Let's try plotting by treatment
plot_pca(so, color_by = 'treatment_s', text_labels = TRUE)

# We can also see genes involved in the the 1st PC by looking at the loadings
# (primary genes whose linear combinations define the principal components)
plot_loadings(so, pc_input = 1)

# Let's see how this "influential" gene (at least as far as PCA tells us) looks by condition
plot_bootstrap(so, 'AT2G34420.1', color_by = 'genotype_variation_s') 

# Let's see how this "influential" gene (at least as far as PCA tells us) looks by treatment
plot_bootstrap(so, 'AT2G34420.1', color_by = 'treatment_s') 

##Step 6: Modeling, testing, and results exploration 
### Differential expression testing with Sleuth 
# Now we need to run a few functions that will test for differential expression (abundance). 
# First we will create a model
so <- sleuth_fit(so, ~genotype_variation_s, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# Now we can get the results of this analysis
full_results <- sleuth_results(so, 'reduced:full', 'lrt',
                               show_all = FALSE)
head(full_results)

# Let's add  Wald test
wald_test <- colnames(design_matrix(so))[2]
so <- sleuth_wt(so, wald_test)

# And start a Shiny Browser
sleuth_live(so)

# NIC3 (https://pubmed.ncbi.nlm.nih.gov/30692220/#&gid=article-figures&pid=figure-5-uid-4)
plot_bootstrap(so, "AT5G23220.1", units = "tpm", color_by = "genotype_variation_s")
