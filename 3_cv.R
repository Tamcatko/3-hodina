
# Zde si dopisujte potrebne knihovny. 
library(readxl)
library(dplyr)
library(tidyr)

#############################
## data qPCR ----

# Nactete data ze souboru data_examples.xlsx, list "qPCR". 
d <- read_excel("data_examples.xlsx", sheet = "qPCR")
# Spocitejte sumarizace pro vsechny sloupce. 
summary(d)
# Prevedte sloupce na faktory, kde to dava smysl. 
Gene = as.factor(d$Gene)
Condition =as.factor(d$Condition)
Replicate = as.factor(d$Replicate)
# Vypiste pouze radky s kondici wt 
filter(d, Condition == "wt")
# Vypiste pouze radky s kondici wt pro geny Fzd2 a Wnt9a
filter(d, Condition == "wt", Gene == "Fzd2" |  Gene == "Wnt9a")

# Spocitejte prumery pro jednotlive kondice u jednotlivych genu. 
d %>% group_by(Condition, Gene) %>% summarise(mean(Expression))

# Prumery vypiste v poradi genu Prickle1, Fzd2, Wnt9a a kondice wt, #1, #2, #3. 

d %>%  group_by(Condition, Gene) %>%  summarise(priemer=mean(Expression)) %>%
  arrange(factor(Gene, levels = c("Prickle1", "Fzd2", "Wnt9a")), factor(Condition, levels = c("wt", "#1", "#2", "#3")))

# ÚKOL NAVÍC: Muzete se pokusit vykreslit jednoduchy graf s prumery. 

# Vyjadrete expresi relativne vuci hodnote wt v kazdem replikatu. 

d %>% group_by(Gene, Replicate) %>% mutate(Expression_rel = Expression / Expression[Condition == "wt"]) %>% ungroup()

# Pro tyto hodnoty vyjadrete prumer. 

d %>% group_by(Gene, Replicate) %>% mutate(Expression_rel = Expression / Expression[Condition == "wt"]) %>% ungroup() %>%  summarise(mean(Expression_rel))

#############################
## data WB quantification ----

# Nactete data ze souboru data_examples.xlsx, list "WB quantification". 


# Vyberte vsechny sloupce obsahujici retezec "Condition".  

de <- read_excel("data_examples.xlsx", sheet = "WB quantification")
de %>%  select(contains("Condition"))
str(de)

# Spocitejte prumery a sd pro vsechny sloupce zacinajici na "Condition". 

de %>% select(grep("Condition", names(de))) %>% summarise_all(list(mean = mean, median = median))

# Prevedte data do long formatu. 

tab<- de %>% pivot_longer( cols = 2:7)

# Spocitejte prumery a sd pro jednotlive kondice. 

tab %>% group_by(name) %>% summarise(mean(value), sd(value))

# Seradte je dle prumeru sestupne. 

tab %>% group_by(name) %>% summarise(m = mean(value), sd(value)) %>% arrange(desc(m))

# Vypiste vsechny kondice, u kterych je prumer WB nad 1 a sd pod 0.2.

tab %>% group_by(name) %>% summarise( m = mean(value), s = sd(value)) %>% filter(m > 1 & s < 0.2)

# Vytvorte (v dlouhych datech) novou kategorialni velicinu vyjadrujici, zda je hodnota WB mensi ci vetsi nez globalni prumer. 

tab <- tab %>% select(-"pomer")
tab <- tab %>% mutate(pomer = case_when(value > mean(value) ~ "vacsia",value < mean(value) ~ "mensia") )


# Kolik je hodnot nad globalnim prumerem u jednotlivych kondic? 


tab %>%  filter(pomer== "vacsia") %>% group_by(name) %>% count(pomer)

#############################
### DATA ELOVL ----

# Nactete data ELOVLs.xlsx, list data_wide. 

kek <- read_excel("ELOVLs.xlsx", sheet = "data_wide")
str(kek)
kek
# Prevedte je na long format parove usporadany, tzn. tumor, non-tumor vedle sebe.

kak<- kek %>% pivot_longer(2:11) %>% 
  separate(name, sep = "_", into = c("gene", "stadium")) %>%
  pivot_wider(names_from = stadium, values_from = value)

# Vypocitejte promennou ratio jako podil tumot/non.tumor. 

kak <- kak %>%  mutate(ratio = tumor / non.tumor)

# Pridejte k datum prumernou hodnotu ratio pro kazdy gen. 

kak <- kak %>%  group_by(gene) %>% mutate(ratio.gene = mean(ratio, na.rm = TRUE)) %>% ungroup() 

# Vypocitejte prumer a sd z ratio pro kazdy gen samostatne, priradte do samostatne tabulky   

tab_new <- kak %>% group_by(gene) %>% summarise( mean_ratio = mean(ratio.gene), sd_ratio = sd(ratio.gene)) %>% ungroup()

# a napojte ji do puvodnich dat v dlouhem formatu. 

kak <- kak %>% left_join(tab_new)

## NORMALIZACE DAT SPECIFICKY PRO JEDNOTLIVE GENY (TJ. PRUMER A SD BUDOU ROZDILNE PRO RUZNE GENY)

# Transformujte ratio pomoci dekadickeho logaritmu.

kak <- kak %>%  mutate(log_ratio = log10(ratio))

# Normalizujte = odectete prumer a podelte smerodatnou odchylkou prislusnou danemu genu. 

kak <- kak %>% group_by(gene) %>% mutate (dalsi= ( ratio - mean) / sd_ratio) %>% ungroup()

# Vypoctete prumer a smerodatnou odchylku normalizovanych hodnot pro kazdy gen. 

kak %>%  group_by(gene) %>% 

## UKOL NAVIC: 
# Provedte normalizaci znovu tak, ze gen normalizujete jen pokud prislusny gen obsahuje alespon 22 valinich hodnot. 
