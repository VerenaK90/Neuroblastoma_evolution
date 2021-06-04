## Reproduce Fig. 6
##############################################################################################################################################
### Learn the dynamics of tumor progression
##############################################################################################################################################
### Load libraries and settings

## if you also want to re-run Mobster, uncomment the following line.
#source(paste0(custom.script.directory, "Mobster.R"))

load(paste0(rdata.directory, "MRCA_timing.RData"))
load(paste0(rdata.directory,"Clonal_mutations_different_ploidies.RData"))
load(paste0(rdata.directory,"Vafs_all_tumors.RData"))

##### take the mutation rate from the fit of tumor initiation
fits <- read.csv(paste0(output.directory, "Figure5/Expansion_decay_continuous_evol.csv"))
mutation.rate <- c(mean(fits$par_mu*2), sd(fits$par_mu*2))

source(paste0(custom.script.directory, "Compute_evolutionary_parameters_from_growth_model.R"))

## store source data
wb <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure6/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}


##############################################################################################################################################
## Figure 5a: effective mutation rate


to.plot <- data.frame(Age=subset$Age/365, mu.eff.mean = effective.mutation.rates[1,],
                      Division.rate=division.rate[1,],
                      Subtype=subset$Telomere.maintenance.mechanism,
                      Effective.division.rate=(1-deltas[1,]),
                      Location=subset$Location)

to.plot$Location[to.plot$Location %in% c("Relapse tumor", "Relapse metastasis")] <- "Relapse"

to.plot$Location <- factor(to.plot$Location, levels=c("Primary", "Metastasis", "Relapse"))
to.plot <- to.plot[to.plot$Subtype %in% c("MNA", "TERT", "ALT"),]
to.plot$Subtype <- factor(to.plot$Subtype, levels=c("MNA", "TERT", "ALT"))

addWorksheet(wb, "a")
writeData(wb, "a", to.plot) 

pdf(paste0(panel.directory, "Figure_6a.pdf"), width = 5, height=4)

p <- ggplot(to.plot, aes(x=Subtype, y=mu.eff.mean, fill=Location)) + 
  geom_boxplot()+ 
  scale_fill_manual(values=time.colors) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(name = "") + scale_y_log10()

print(p)


dev.off()

##############################################################################################################################################
## Figure 6b: relative loss rates

pdf(paste0(panel.directory, "Figure_6b.pdf"), width = 5, height=4, useDingbats = F)

## to summarize the data appropriately, I compute for each subclass the standard error of the mean

to.plot <- data.frame(delta.mean = sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
  mean(deltas[1,subset$Telomere.maintenance.mechanism==x & subset$Location %in% c("Primary", "Metastasis")])}), 
  delta.sd =sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
    sqrt(sum((deltas[2,subset$Telomere.maintenance.mechanism==x& subset$Location %in% c("Primary", "Metastasis")])^2))/sum(subset$Telomere.maintenance.mechanism==x& subset$Location %in% c("Primary", "Metastasis"))}),
  Subtype=unique(subset$Telomere.maintenance.mechanism))

to.plot <- to.plot[to.plot$Subtype %in% c("MNA", "TERT", "ALT"),]

to.plot$Subtype <- factor(to.plot$Subtype, levels=c("MNA", "TERT", "ALT"))

addWorksheet(wb, "b")
writeData(wb, "b", to.plot) 


p <- ggplot(to.plot, aes(x=Subtype, y=delta.mean, ymin=delta.mean - delta.sd, ymax=delta.mean+delta.sd)) + 
  geom_pointrange()+ scale_y_continuous(limits=c(0,1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)
dev.off()


##############################################################################################################################################
## Figure 6c: division rate

pdf(paste0(panel.directory, "Figure_6c.pdf"), width = 5, height=4, useDingbats = F)

## to summarize the data appropriately, I compute for each subclass the standard error of the mean


to.plot <- data.frame(division.rate.mean = sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
  mean(division.rate[1,subset$Telomere.maintenance.mechanism==x & subset$Location %in% c("Primary", "Metastasis")])}), 
  division.rate.sd =sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
    sqrt(sum((division.rate[2,subset$Telomere.maintenance.mechanism==x & subset$Location %in% c("Primary", "Metastasis")])^2))/sum(subset$Telomere.maintenance.mechanism==x & subset$Location %in% c("Primary", "Metastasis"))}),
  Telomere.maintenance.mechanism=unique(subset$Telomere.maintenance.mechanism))

to.plot <- to.plot[to.plot$Telomere.maintenance.mechanism %in% c("MNA", "TERT", "ALT"),]

to.plot$Telomere.maintenance.mechanism <- factor(to.plot$Telomere.maintenance.mechanism, levels=c("MNA", "ALT", "TERT"))

addWorksheet(wb, "c")
writeData(wb, "c", to.plot) 

p <- ggplot(to.plot, aes(x=Telomere.maintenance.mechanism, y=division.rate.mean, ymin=division.rate.mean - division.rate.sd,
                    ymax=division.rate.mean+division.rate.sd)) + 
  geom_pointrange()+ scale_y_continuous(limits=c(0, 1.1*max(to.plot$division.rate.mean + to.plot$division.rate.sd)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)

dev.off()

##############################################################################################################################################
saveWorkbook(wb, file = paste0(panel.directory, "Source_data_Fig.6.xlsx"), overwrite=T)
