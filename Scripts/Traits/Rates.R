# Load packages and data ----------------------------------------------------
library(tidyverse)


# Read the traits data -----------------------------------------------------

MorphTraits = readRDS(file = "Intermediate/Morphology/AllTraits.rds")



# Computer rate traits (growth and development) ---------------------------

c = filter(MorphTraits, Treatment == "C")
c1 = filter(c, Timepoint == 1)
c3 = filter(c, Timepoint == 3)

s = filter(MorphTraits, Treatment == "S", Timepoint <= 3)
s = bind_rows(s, filter(c1, XP == "PSD") %>% mutate(Treatment = "S"))
s3 = filter(s, Timepoint == 3)

ps = filter(MorphTraits, Treatment == "S", Timepoint > 3)
ps = bind_rows(ps, s3) %>% mutate(Treatment = "PS")

psd = filter(MorphTraits, Treatment == "PSD", Timepoint > 3)
psd = bind_rows(psd, s3) %>% mutate(Treatment = "PSD")

d1 = filter(MorphTraits, Treatment == "D", XP == "HTD")
d1 = bind_rows(d1, filter(c1, XP == "HTD")) %>% mutate(Treatment = "D")

d2 = filter(MorphTraits, Treatment == "D", XP == "PSD")
d2 = bind_rows(d2, filter(c3, XP == "PSD")) %>% mutate(Treatment = "D")

ht = filter(MorphTraits, Treatment == "HT")
ht1 = filter(ht, Timepoint == 1)

htd = filter(MorphTraits, Treatment == "HTD")
htd = bind_rows(htd, ht1) %>% mutate(Treatment = "HTD")

MorphTraitsRGR = bind_rows(c, s, ps, psd, d1, d2, ht, htd)


# Save the data -----------------------------------------------------------

saveRDS(MorphTraitsRGR, file = "Intermediate/Rates/AllTraits.rds")

