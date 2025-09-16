library(dplyr)
library(dplyr)
library(readr)
library(stringr)

repeatmasker_out <- "susScr11.fa.out.gz"

# 把 reference 三列先读成字符，避免括号导致 NA
df <- read_table2(
  repeatmasker_out, skip = 3,
  col_names = c("SW_score","perc_div","perc_del","perc_ins",
                "query","q_begin","q_end","q_left","strand",
                "repeat_id","repeat_class",
                "r_begin","r_end","r_left","id_graphic"),
  col_types = cols(
    .default = col_character(),
    SW_score = col_double(),
    perc_div = col_double(),
    perc_del = col_double(),
    perc_ins = col_double(),
    q_begin  = col_integer(),
    q_end    = col_integer()
  )
)

df2 <- df %>%
  mutate(
    r1 = r_begin, r2 = r_end, r3 = r_left,
    r_begin_n = case_when(
      strand == "+" ~ as.numeric(r1),
      strand == "C" ~ as.numeric(r3)
    ),
    r_end_n = as.numeric(r2),
    r_left_n = case_when(
      strand == "+" ~ as.numeric(str_remove_all(r3, "[()]")),
      strand == "C" ~ as.numeric(str_remove_all(r1, "[()]"))
    ),
    aln_length   = abs(r_end_n - r_begin_n) + 1,
    total_length = r_end_n + r_left_n,
    coverage = aln_length / total_length
  )

#如果要筛高覆盖（例：>0.8），可用：
df2_high <- df2 %>% filter(coverage ==1)

best_te <- df2_high %>%
  filter(str_detect(repeat_class, "^(SINE|LINE|LTR|DNA)/"))

write.table(best_te,file='TE_ref.bed',quote=F,row.names=F,sep='\t')

