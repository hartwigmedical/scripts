source("libSVPaper.R")


bind_rows(
  sva_svs %>%
    mutate(length=InexactHOEnd - InexactHOStart, type="Inexact homology") %>%
    group_by(type, length) %>%
    summarise(count=n()),
  sva_svs %>%
    mutate(length=str_length(Homology), type="Exact homology") %>%
    replace_na(list(length=0)) %>%
    group_by(type, length) %>%
    summarise(count=n())) %>%
  ggplot() +
  geom_line(data=data.frame(length=0:50) %>% mutate(count=662929 * (0.25**length)) %>% filter(count > 100 & count < 750000), colour="grey") +
  geom_line(data=data.frame(length=0:50) %>% mutate(count=10**6 * ((1/1.93)**length)) %>% filter(count > 100 & count < 750000), colour="grey") +
  geom_line(data=data.frame(length=13:50) %>% mutate(count=10**3 * ((1/1.07)**length)), colour="orange") +
  geom_line(data=data.frame(length=13:50) %>% mutate(count=10**2.7 * ((1/1.04)**length)), colour="grey") +
  geom_line(data=data.frame(length=13:50) %>% mutate(count=10**3.2 * ((1/1.11)**length)), colour="grey") +
  aes(x=length, y=count, colour=type) +
  geom_point() +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 50)) +
  labs(title="PCAWG patterns of SVs Fig 5d")

bind_rows(
  sva_svs %>%
    filter(ChrEnd != 0) %>%
    mutate(length=InexactHOEnd - InexactHOStart) %>%
    mutate(length=ifelse(length==0, -str_length(sva_svs$InsertSeq), length)) %>%
    replace_na(list(length=0)) %>%
    group_by(length, ResolvedType) %>%
    summarise(count=n()) %>%
    group_by(ResolvedType) %>%
    filter(sum(count) > 10000)) %>%
  ggplot() +
  aes(x=length, y=count, colour=ResolvedType, linetype=ResolvedType) +
  geom_vline(xintercept=0, colour="grey") +
  geom_line() +
  scale_y_log10() +
  scale_x_continuous(limits = c(-100, 50)) +
  labs(title="PCAWG patterns of SVs Fig 5d faceted by ResolvedType")
# What complex SV patterns does SSA repair explain?





