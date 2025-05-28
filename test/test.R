# Quickly testing the thing to make sure it works
# TODO convert to formal unit tests

seq_data <- readRDS(here::here("test/aligned-sequence-data.Rds")) |>
  dplyr::filter(subtype == "H1N1")

seq_vec <- seq_data$protein_sequence
names(seq_vec) <- seq_data$full_name

str(seq_vec)
head(seq_vec)

R_scripts <- list.files(here::here("R"), full.names = TRUE)
for (f in R_scripts) source(f)

ac_map <- readRDS(here::here("test/A(H1N1).Rds"))
tree <- readRDS(here::here("test/ml-tree-list.Rds"))
my_tree <- tree$H1N1

calculate_distance(seq_vec, "sdvfasv") # errors
calculate_distance(seq_vec, "hamming")
calculate_distance(seq_vec, "dl")
calculate_distance(seq_vec, "p-epitope") # errors
calculate_distance(seq_vec, "p-epitope", subtype = "A(H1N1)")
calculate_distance(seq_vec, "p-epitope", subtype = "A(H1N1)", mode = "anderson")
calculate_distance(seq_vec, "p-all-epitope") # errors
calculate_distance(seq_vec, "p-all-epitope", subtype = "A(H1N1)")
calculate_distance(seq_vec, "cophenetic") # errors
calculate_distance(seq_vec, "cophenetic", tree = my_tree)
calculate_distance(seq_vec, "cartographic") # errors
calculate_distance(seq_vec, "cartographic", map = ac_map)
calculate_distance(seq_vec, "temporal")

out <- calculate_distance(seq_vec, "hamming")
tidy_dist_mat(out)
tidy_dist_mat(out, "col1", "col2", TRUE)
