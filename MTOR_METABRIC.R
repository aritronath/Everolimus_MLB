

##### Function to create risk groups from ENDORSE predicted risk ----------------------
MakeRiskGrp <- function (risk.pred, ngrps=3, thresh=c(1,2)) {
  risk.grps <- array(length(risk.pred))
  
  if (ngrps == 2) {
    risk.grps[which(risk.pred <= 0.5*(thresh[1] + thresh[2]))] <- "Low Risk"
    risk.grps[which(risk.pred > 0.5*(thresh[1] + thresh[2]))] <- "High Risk"
    risk.grps <- factor(risk.grps, levels=c("Low Risk", "High Risk"))
  }
  
  if (ngrps == 3) {
    risk.grps[which(risk.pred <= thresh[1])] <- "Low Risk"
    risk.grps[which(risk.pred > thresh[1] & risk.pred < thresh[2])] <- "Medium Risk"
    risk.grps[which(risk.pred >= thresh[2])] <- "High Risk"
    risk.grps <- factor(risk.grps, levels=c("Low Risk", "Medium Risk", "High Risk"))
  }
  
  return(risk.grps)
}

#### Match gene expression matrices by rownames and perform ComBat batch correction ---------------------
library(sva)
DoCombat <- function(mat1, mat2) {
  
  comG <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[comG, ]
  mat2 <- mat2[comG, ]
  
  combined_mat <- cbind(mat1, mat2)
  batch_labels <- c(rep("A", ncol(mat1)), rep("B", ncol(mat2)))
  
  corrected_mat <- ComBat(combined_mat, batch = batch_labels)
  
  return(corrected_mat)
}

########### Match matrices and train the ENDORSE model on batch-corrected data ####-------------------------
library(GSVA)
library(survival)

metabric.s <- readRDS(file="~/Desktop/Biomarkers/Scripts/metabric.s.RDS")
metabric.os <- readRDS(file="~/Desktop/Biomarkers/Scripts/metabric.os.RDS")
metabric.os.event <- readRDS(file="~/Desktop/Biomarkers/Scripts/metabric.os.event.RDS")
ENDORSE <- readRDS(file="~/Desktop/Biomarkers/Scripts/ENDORSE.RDS")
H_ESTR_EARLY <- readRDS(file="~/Desktop/Biomarkers/Scripts/H_ESTR_EARLY.RDS")

MBS.ENDORSE <- gsva(metabric.s, 
                    gset.idx.list = list("ENDORSE"=ENDORSE, "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY), 
                    method="ssgsea", kcdf="Gaussian")

df <- data.frame('Time'=metabric.os, 
                 'Event'=as.numeric(metabric.os.event), 
                 t(MBS.ENDORSE))

cfit.6 = coxph(Surv(Time, Event) ~ ENDORSE + HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df)

NEWDAT.ssgsea <- t(MBS.ENDORSE)

# predict
NEWDAT.pred <- predict(cfit.6, data.frame(NEWDAT.ssgsea), type="risk") 
NEWDAT.pred.grp <- MakeRiskGrp(NEWDAT.pred)

print(paste("ENDORSE score = ", round(NEWDAT.pred, 2)))
print(paste("ENDORSE risk = ", as.character(NEWDAT.pred.grp)))

######  Predict MTOR response ####################------------------------------
MTOR.sig.genes <- as.character(readRDS(file="~/Desktop/Biomarkers/Scripts/MTOR.sig.genes.RDS"))
MTOR.pre.train <- readRDS(file="~/Desktop/Biomarkers/Scripts/MTOR.pre.train.RDS")
MTOR.resp.pre.train <- readRDS(file="~/Desktop/Biomarkers/Scripts/MTOR.resp.pre.train.RDS")

temp <- DoCombat(t(MTOR.pre.train), metabric.s)
train <- match(rownames(MTOR.pre.train), colnames(temp))
test <- c((length(train) + 1) : ncol(temp))

# Train model on batch corrected data 
library(caret)

# List of seeds for reproducibility 
seeds_list <- split(seq(1, length(train)*90, 1), c(1:length(train)))
seeds_list[[length(train)+1]] <- 1234

MTOR.classifier <- intersect(MTOR.sig.genes, rownames(temp))

message(paste("Found", length(MTOR.classifier), "of", length(MTOR.sig.genes), "in the integrated dataset"))

fitControl <- trainControl(method = "LOOCV", seeds = seeds_list)
rfFit.full <- train(as.factor(MTOR.resp.pre.train) ~ ., 
                    data = cbind(MTOR.resp.pre.train, t(temp[MTOR.classifier, train])), 
                    method = "rf", tuneGrid=data.frame(mtry=seq(1, 90, 1)),
                    trControl = fitControl,
                    verbose = T)

nd <- temp[, test]
rownames(nd) <- rownames(temp)

pred.MTOR <- predict(rfFit.full, newdata = t(nd), type = "prob") #1=Non-responder 2=Responder

print(paste("Probability of MTOR response =", round(pred.MTOR$`2`, 2)))

# background probabilities 
#k <- which(MTOR.resp.pre.train == '2')
pred.Train <- predict(rfFit.full, t(temp[MTOR.classifier, train]), type="prob")
summary(pred.Train$`2`)

plot(density(pred.Train$`2`))
plot(density(pred.MTOR$`2`))

cor.test(NEWDAT.pred, pred.MTOR$`2`)

boxplot(pred.MTOR$`2` ~ NEWDAT.pred.grp)

pred.MTOR$`2`

n1 <- data.frame(table(NEWDAT.pred.grp))
n2 <- data.frame(table(NEWDAT.pred.grp[which(pred.MTOR$`2` > 0.75)]))

n2$Freq*100/n1$Freq
#9.701493 15.109890 40.298507

df.MTOR <- data.frame(df, "ENDORSE_score" = NEWDAT.pred,
                      "ENDORSE_cat" = as.factor(NEWDAT.pred.grp),
                      "Pred.mTOR" = pred.MTOR$`2`,
                      "Pred.mTOR_cat" = ifelse(pred.MTOR$`2` >= 0.75, "Responder", "Non-responder")
                      )

summary(aov(Pred.mTOR ~ ENDORSE_cat, data=df.MTOR))
TukeyHSD(aov(Pred.mTOR ~ ENDORSE_cat, data=df.MTOR))

#Medium Risk-Low Risk  0.02131554 0.003109576 0.0395215 1.679969e-02
#High Risk-Low Risk    0.11084080 0.077636142 0.1440455 7.582823e-14
#High Risk-Medium Risk 0.08952526 0.056073919 0.1229766 1.599518e-09


pdf(file="~/Desktop/Biomarkers/mTOR/Figure_5_ENDORSE_mtor.pdf", width=7, height=7)
ggplot(df.MTOR, aes(x=ENDORSE_cat, y=Pred.mTOR, fill=ENDORSE_cat)) +
  geom_boxplot() +
  geom_jitter(col = "grey25", width=0.15) +
  theme_classic(base_size = 14) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x="ENDORSE risk group", y="Probability of response")
dev.off()

pb <- txtProgressBar(min=1, max=nrow(metabric.s), char="+")
ThePvals <- matrix(nrow=nrow(metabric.s), ncol=4)
for (i in 1:nrow(metabric.s)) {
  aov.temp <- TukeyHSD(aov(metabric.s[i,] ~ df.MTOR$Pred.mTOR_cat))
  ThePvals[i, ] <- aov.temp$`df.MTOR$Pred.mTOR_cat`
  setTxtProgressBar(pb, i)  
}
close(pb)

msig_up <- which(ThePvals[,4] < 0.05 & ThePvals[,1] >= log(2))
msig_dn <- which(ThePvals[,4] < 0.05 & ThePvals[,1] <= log(0.5))

msig_g_up <- rownames(metabric.s)[msig_up]
msig_g_dn <- rownames(metabric.s)[msig_dn]

library(gprofiler2)
gp_msig_up <- gost(query = msig_g_up, 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = F, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = FALSE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = c("GO:BP", "REAC", "KEGG"), as_short_link = FALSE)


gp_msig_dn <- gost(query = msig_g_dn, 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = F, exclude_iea = FALSE, 
                   measure_underrepresentation = FALSE, evcodes = FALSE, 
                   user_threshold = 0.05, correction_method = "g_SCS", 
                   domain_scope = "annotated", custom_bg = NULL, 
                   numeric_ns = "", sources = c("GO:BP", "REAC", "KEGG"), as_short_link = FALSE)

gop3 <- gostplot(gp_msig_up, capped = FALSE, interactive = F)  #Higher in responders
gop4 <- gostplot(gp_msig_dn, capped = FALSE, interactive = F)  #Higher in non-responders

write.table(apply(gp_msig_up$result, 2, as.character), sep="\t", file="~/Desktop/Biomarkers/mTOR/gp_msig_up.txt", quote=F, row.names=F)
write.table(apply(gp_msig_dn$result, 2, as.character), sep="\t", file="~/Desktop/Biomarkers/mTOR/gp_msig_dn.txt", quote=F, row.names=F)

publish_gostplot(gop3, highlight_terms = c("GO:0006955", "GO:0001775", "GO:0048583",  
                                           "KEGG:04062", "KEGG:04660", "REAC:R-HSA-1280215", "REAC:R-HSA-9664873"), 
                 width = 7.5, height = 6.5, filename = "~/Desktop/Biomarkers/mTOR/gp_msig_up.pdf")

publish_gostplot(gop4, highlight_terms = c("REAC:R-HSA-9018519", "GO:0045880", "GO:0007224"), 
                 width = 7.5, height = 5.5, filename = "~/Desktop/Biomarkers/mTOR/gp_msig_dn.pdf")

upload_GMT_file(gmtfile = "~/Desktop/msigdb/h.all.v7.4.symbols.gmt")

gp_msig_up.h <- gost(query = msig_g_up, 
                   ordered_query = FALSE, 
                   multi_query = FALSE, significant = F, exclude_iea = FALSE, 
                   measure_underrepresentation = FALSE, evcodes = FALSE, 
                   user_threshold = 0.05, correction_method = "g_SCS", 
                   organism = "gp__cKdl_eD9N_kI0")


gp_msig_dn.h <- gost(query = msig_g_dn, 
                     ordered_query = FALSE, 
                     multi_query = FALSE, significant = F, exclude_iea = FALSE, 
                     measure_underrepresentation = FALSE, evcodes = FALSE, 
                     user_threshold = 0.05, correction_method = "g_SCS", 
                     organism = "gp__cKdl_eD9N_kI0")

gop5 <- gostplot(gp_msig_up.h, capped = FALSE, interactive = F)  #Higher in responders
gop6 <- gostplot(gp_msig_dn.h, capped = FALSE, interactive = F)  #Higher in non-responders

write.table(apply(gp_msig_up.h$result, 2, as.character), sep="\t", file="~/Desktop/Biomarkers/mTOR/gp_msig_up.h.txt", quote=F, row.names=F)
write.table(apply(gp_msig_dn.h$result, 2, as.character), sep="\t", file="~/Desktop/Biomarkers/mTOR/gp_msig_dn.h.txt", quote=F, row.names=F)

publish_gostplot(gop5, highlight_terms = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",  
                                           "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 
                 width = 7.5, height = 5.5, filename = "~/Desktop/Biomarkers/mTOR/gp_msig_up_h.pdf")

publish_gostplot(gop6, highlight_terms = c("HALLMARK_ESTROGEN_RESPONSE_EARLY"), 
                 width = 7.5, height = 5.5, filename = "~/Desktop/Biomarkers/mTOR/gp_msig_dn_h.pdf")
