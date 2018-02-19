library(tidyverse)
library(assertthat)
options(stringsAsFactors = FALSE)

# lowercase = reverse complement
base_edges <- data.frame(
	from=c("A", "B", "C", "D", "e", "d", "c", "b"),
	to=  c("B", "C", "D", "E", "d", "c", "b", "a")
)

breakpoints <- data.frame(
		# Breakend ordering
		bp1_be1=c(1, 1, 1),
		bp1_be2=c(4, 3, 2),
		bp2_be1=c(2, 2, 3),
		bp2_be2=c(3, 4, 4),
		desc=c("Enclosing", "Overlapping", "Neighbouring")
	) %>% merge(expand.grid(
		bp1_dir1=c("+", "-"),
		bp1_dir2=c("+", "-"),
		bp2_dir1=c("+", "-"),
		bp2_dir2=c("+", "-"))) %>%
	mutate(
		bp1=ifelse(bp1_dir1=="+" & bp1_dir2=="-", "DEL",
			ifelse(bp1_dir1=="-" & bp1_dir2=="+", "DUP",
			ifelse(bp1_dir1=="-" & bp1_dir2=="-", "INV-", "INV+")))) %>%
	mutate(
		bp2=ifelse(bp2_dir1=="+" & bp2_dir2=="-", "DEL",
			ifelse(bp2_dir1=="-" & bp2_dir2=="+", "DUP",
			ifelse(bp2_dir1=="-" & bp2_dir2=="-", "INV-", "INV+"))))
#' Edges induced by the given breakpoint
get_edges_for <- function(pos1, pos2, dir1, dir2) {
	flip <- c("a"="A", "b"="B", "c"="C", "d"="D", "e"="E", "A"="a", "B"="b", "C"="c", "D"="d", "E"="e")
	from=ifelse(dir1 == "+", c("A", "B", "C", "D")[pos1], c("b", "c", "d", "e")[pos1])
	to=ifelse(dir2=="-", c("B", "C", "D", "E")[pos2], c("a", "b", "c", "d")[pos2])
	data.frame(
		from=c(from, flip[to]),
		to=c(to, flip[from]))
}
#' Resultant edge graph for the two given breakpoints
get_graph_edges_for <- function(df) {
	base_edges %>% mutate(is_bp1=0, is_bp2=0) %>%
		bind_rows(get_edges_for(df$bp1_be1, df$bp1_be2, df$bp1_dir1, df$bp1_dir2) %>%
			mutate(is_bp1=1, is_bp2=0)) %>%
		bind_rows(get_edges_for(df$bp2_be1, df$bp2_be2, df$bp2_dir1, df$bp2_dir2) %>%
			mutate(is_bp1=0, is_bp2=1))
}
traverse <- function(edges) {
	step <- data.frame(
		start=c("A", "e"),
		current=c("A", "e"),
		path=c("A", "e"),
		end=NA_character_,
		used_bp1=0,
		used_bp2=0
	)
	paths <- NULL
	while(nrow(step) > 0) {
		step <- step %>% left_join(edges, by=c("current"="from")) %>%
			mutate(end=ifelse(to %in% c("a", "E"), to, end),
				path=paste0(path, to),
				current=to,
				used_bp1=used_bp1 + is_bp1,
				used_bp2=used_bp2 + is_bp2) %>%
			select(-to, -is_bp1, -is_bp2) %>%
			filter(used_bp1 <= 1 & used_bp2 <= 1)
		paths <- paths %>% bind_rows(step %>% filter(!is.na(step$end) & (used_bp1 > 0 | used_bp2 > 0)))
		step <- step %>% filter(is.na(step$end))
	}
	return(paths %>% mutate(
		has_telomere=start == "A" | end == "A" | start == "a" | end == "a",
		has_centromere=start == "E" | end == "E" | start == "e" | end == "e"
	))
}
analyse_paths <- function(paths) {
	# Single consistent path
	replication_paths <- paste0((paths %>%
			filter(start=="A") %>%
			filter(has_telomere & has_centromere & used_bp1 > 0 & used_bp2 > 0))$path,
		collapse = ", ")
	# Pairs of paths that are consistent
	paired_paths <- full_join(
		paths %>% filter(used_bp1==1 & used_bp2==0) %>% mutate(hack=1),
		paths %>% filter(used_bp1==0 & used_bp2==1) %>% mutate(hack=1),
		by="hack")
	normal_paths <- paired_paths %>%
		filter(has_telomere.x & has_centromere.x & has_telomere.y & has_centromere.y) %>%
		filter(start.x == "A" & start.y == "A")
	normal_paths_str = paste0(
		paste0(unique(normal_paths$path.x), collapse=","),
		"/",
		paste0(unique(normal_paths$path.y), collapse=","))

	inv_paths <- paired_paths %>%
		filter((has_telomere.x & !has_centromere.x & !has_telomere.y & has_centromere.y) |
				(!has_telomere.x & has_centromere.x & has_telomere.y & !has_centromere.y))
	inv_paths_str = paste0(
		paste0(unique(inv_paths$path.x), collapse=","),
		"/",
		paste0(unique(inv_paths$path.y), collapse=","))
	return(data.frame(replication_paths=replication_paths,
		alt_hap_paths=normal_paths_str,
		inv_hap_paths=inv_paths_str
		))
}

results <- breakpoints %>%
	mutate(row_number_hack=row_number()) %>%
	group_by(row_number_hack) %>%
	do({
		paths <- traverse(get_graph_edges_for(.))
		bind_cols(., analyse_paths(paths))
	}) %>%
	ungroup() %>%
	select(-row_number_hack) %>%
	arrange(bp1, bp2, desc)



