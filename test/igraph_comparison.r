# igraph_comparison.r

# install.packages(c("igraph", "readr", "dplyr", "microbenchmark"))

library(igraph)
library(readr)
library(dplyr)
library(tibble)
library(microbenchmark)

dir.create("data/r", showWarnings = FALSE, recursive = TRUE)
meta <- readr::read_csv("data/_meta.csv", show_col_types = FALSE)

RES <- 1.0

summ <- tibble(
  name=character(), n=integer(), m=integer(),
  r_time_ms=double(), communities=integer(), modularity=double(),
  stringsAsFactors = FALSE
)

set.seed(08540)

for (i in 1:nrow(meta)) {
    name <- meta$name[i]
    df <- readr::read_csv(paste0("data/", name, ".csv"), show_col_types = FALSE)
    n <- max(c(df$src, df$dst))
    g <- graph_from_data_frame(
        df, directed = FALSE, vertices = tibble(name = 1:n)
    )

    mb <- microbenchmark(
        res = cluster_leiden(g, resolution = RES, objective_function = "modularity"),
        times = 5
    )
    cl <- cluster_leiden(g, resolution = RES, objective_function = "modularity")

    memb <- membership(cl)
    out <- tibble(node = as.integer(names(memb)), community = as.integer(memb))
    out <- out[order(out$node), ]
    readr::write_csv(out, paste0("data/r/", name, "_r_partition.csv"))

    summ <- rbind(
        summ, tibble(
            name = name,
            n = vcount(g),
            m = ecount(g),
            r_time_ms = median(mb$time) / 1e6,
            communities = length(sizes(cl)),
            mod = modularity(g, membership(cl))
        )
    )

  cat("R done: ", name, " | Q=", round(modularity(g, membership(cl)), 6),
      " | k=", length(sizes(cl)),
      " | ", round(median(mb$time)/1e6, 1), " ms\n", sep="")
}

readr::write_csv(summ, "test/r/_summary_r.csv")
cat("Wrote R summary to test/r/_summary_r.csv\n")
