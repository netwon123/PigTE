library(dplyr)
library(readr)
library(stringr)

# 替换为你的 RepeatMasker .out 文件路径
repeatmasker_out <- "panel_ins.fa.out"

# 读取 RepeatMasker 输出文件（跳过前三行注释）
df <- read_table2(repeatmasker_out, skip = 3,
                  col_names = c("SW_score", "perc_div", "perc_del", "perc_ins",
                                "query", "q_begin", "q_end", "q_left", "strand",
                                "repeat_id", "repeat_class",
                                "r_begin", "r_end", "r_left", "id_graphic"))

# 转换数值类型 + 添加匹配长度
df <- df %>%
  mutate(q_begin = as.integer(q_begin),
         q_end = as.integer(q_end),
         SW_score = as.numeric(SW_score),
         match_len = q_end - q_begin)

# 对每个 query（SV）选出：最长 match_len，其次最大 SW_score
best_te <- df %>%
  group_by(query) %>%
  filter(match_len == max(match_len)) %>%
  filter(SW_score == max(SW_score)) 

best_te <- best_te %>%
  mutate(
    q_begin = as.integer(q_begin),
    q_end = as.integer(q_end),
    q_left = as.integer(str_remove_all(q_left, "[()]")),
    match_len = q_end - q_begin,
    sv_length = q_end + q_left,
    match_span = match_len / sv_length
  )

df_te <- best_te %>%
  filter(match_span > 0.8)




# 输出信息（可根据需要修改）
final_result <- df_te %>%
  select(sv_id = query, sv_start = q_begin, sv_end = q_end,
         te_id = repeat_id, te_class = repeat_class,
         match_len, SW_score)

# 写出结果
write_tsv(final_result, "Best_TE_annotation_per_SV_non_ref-TE.tsv")



