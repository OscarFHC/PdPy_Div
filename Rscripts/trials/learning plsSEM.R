data("spainfoot")

Attack = c(0, 0, 0)
Defense = c(0, 0, 0)
Success = c(1, 1, 0)
# path matrix created by row binding
foot_path = rbind(Attack, Defense, Success)
# add column names (optional)
colnames(foot_path) = rownames(foot_path)

foot_blocks = list(1:4, 5:8, 9:12)

innerplot(foot_path)

foot_modes = c("A", "A", "A")

data("spainfoot")
head(spainfoot)
foot_pls = plspm(spainfoot, foot_path, foot_blocks, modes = foot_modes)
foot_pls$inner_model

summary(foot_pls)

head(HNF_Bac_A)
Envi_pca <- prcomp(HNF_Bac_A[,c("Temp", "Sal", "PAR", "NO2", "NO3", "DIN", "PO3")], center = TRUE,scale. = TRUE)
summary(Envi_pca)
