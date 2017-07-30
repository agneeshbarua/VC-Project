library(dplyr)
library(readr)
library(readxl)
library(tibble)
library(tidyr)
library(DBI)
library(reshape2)
library(xlsx)

#give name to data set using table
Snake2 = read_csv("Snakes.csv")
Snake2

#define what data.frame will be made of
Snake2 = tbl_df(Snake2)

#Chain process to selec stuff needed
Snake2 = Snake2 %>% select(Snake, References, `3FTx`, BPP, CRISP, CTL, GF, LAAO, OHA, PLA2, SVMP, SVSP)
Snake2

#making new csv using the changes data.frame
write_csv(Snake2, "Snake no NA.csv")
