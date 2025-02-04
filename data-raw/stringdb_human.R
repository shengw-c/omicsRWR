## code to prepare `stringdb_human` dataset goes here
stringdb_score_threshold = 900

## download data from StringDB
stringdb = read_delim("https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz", delim = " ")
stringdb_info = read_delim("https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz", delim = "\t")
proteins = stringdb_info$preferred_name

stopifnot(!any(duplicated(proteins)))
names(proteins) = stringdb_info$`#string_protein_id`
stringdb = stringdb %>%
  mutate(protein1 = proteins[protein1],
         protein2 = proteins[protein2])
stopifnot(all(complete.cases(stringdb)))

sub_stringdb = stringdb %>% filter(combined_score > stringdb_score_threshold) %>% mutate(weight = combined_score / 1000)
net = graph_from_data_frame(sub_stringdb, directed = FALSE) %>% igraph::simplify(edge.attr.comb = "max")

## keep largest connect community
stringdb_human = largest_component(net)

## column-wise normalization and convert to sparse matrix
stringdb_human_norm = norm_W(stringdb_human)

### proc some file for example RWR input
testInput = read.table("data-raw/test_input.tsv", sep = "\t", stringsAsFactors = F, header = TRUE) %>%
  rename(score = score_rna) %>% select(gene_name, score)

testInput = tibble(gene_name = V(stringdb_human)$name) %>%
  left_join(testInput, by = "gene_name", relationship = "one-to-one") %>%
  replace(is.na(.), 0)

usethis::use_data(stringdb_human_norm, testInput, overwrite = TRUE)
