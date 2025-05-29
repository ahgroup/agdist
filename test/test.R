# Quickly testing the thing to make sure it works
# TODO convert to formal unit tests

seq_data <- readRDS(here::here("test/aligned-sequence-data.Rds")) |>
  dplyr::filter(subtype == "H1N1") |>
  dplyr::select(strain_type, strain_subtype, strain_name = full_name, protein_sequence)

my_sequences <- seq_data$protein_sequence
names(my_sequences) <- seq_data$strain_name
head(my_sequences)

R_scripts <- list.files(here::here("R"), full.names = TRUE)
for (f in R_scripts) source(f)

my_map <- readRDS(here::here("test/A(H1N1).Rds"))
tree <- readRDS(here::here("test/ml-tree-list.Rds"))
my_tree <- tree$H1N1

calculate_distance(my_sequences, "sdvfasv") # errors
calculate_distance(my_sequences, "hamming")
calculate_distance(my_sequences, "dl")
calculate_distance(my_sequences, "p-epitope") # errors
calculate_distance(my_sequences, "p-epitope", subtype = "A(H1N1)")
calculate_distance(my_sequences, "p-epitope", subtype = "A(H1N1)", mode = "anderson")
calculate_distance(my_sequences, "p-all-epitope") # errors
calculate_distance(my_sequences, "p-all-epitope", subtype = "A(H1N1)")
calculate_distance(my_sequences, "cophenetic") # errors
calculate_distance(my_sequences, "cophenetic", tree = my_tree)
calculate_distance(my_sequences, "cartographic") # errors
calculate_distance(my_sequences, "cartographic", map = my_map)
calculate_distance(my_sequences, "temporal")

out <- calculate_distance(my_sequences, "hamming")
tidy_dist_mat(out)
tidy_dist_mat(out, "col1", "col2", TRUE)
