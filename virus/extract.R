library(tidyverse)
library(readr)

df = read_tsv("virushostdb.tsv") %>%
  filter(`host tax id`==9606) %>%
  separate_rows(`refseq id`, sep="[ ,]+")

writeLines(con="viralhost.table.sql", c("
  DROP TABLE IF EXISTS viralhost;
  CREATE TABLE viralhost
  (
    virus_tax_id int NOT NULL,
    virus_name varchar(255),
    virus_lineage varchar(255),
    refseq_id varchar(32),
    pmid varchar(32),
    index(refseq_id)
  );
  INSERT INTO viralhost (virus_tax_id, virus_name, virus_lineage, refseq_id, pmid) VALUES ",
  paste(
    df %>%
      mutate(line=paste0("(\"", `virus tax id`, "\",\"", `virus name`, "\",\"", `virus lineage`, "\",\"", `refseq id`, "\",\"", ifelse(is.na(pmid), "", pmid), "\")")) %>%
      pull(line),
    collapse=",\n"),
  ";"
))
