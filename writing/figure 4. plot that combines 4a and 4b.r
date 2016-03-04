
# The color-blind friednly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#check power using true, imputed founder carrier, minimum offspring carrier
#three versions -- command and rare and super rare #2 for common, 7 for rare, 39 for super rare
result <- read.csv("figure 4a. fouder imputatuion estimate nontransmitted by population f see power.csv", header=T)
result <- result %>% gather(key="method", value="p.value", -c(1:8))
result.plot <- result %>% group_by(f, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05, na.rm=T))
result.plot$method <- factor(result.plot$method, levels=c("result.trap.3c.noimpute","result.trap.3c.impute.founder.pop.f","result.trap.3c.impute.founder.pop.f.off.10",
                                                          "result.trap.3c.impute.founder.pop.f.off.20","result.trap.3c.impute.founder.pop.f.off.50",
                                                          "result.trap.3c.impute.founder.pop.f.off.n10","result.trap.3c.impute.founder.pop.f.off.n20",
                                                          "result.trap.3c.impute.founder.pop.f.off.n50","result.trap.3c.impute.founder.sample.f"),
                             labels=c("w.o. imputation", "pop.f", "pop.f off 10%",
                                      "pop.f off 20%", "pop.f off 50%", "pop.f off -10%"
                                      , "pop.f off -20%", "pop.f off -50%", "sample.f"))
result.plot$f <- factor(result.plot$f, levels=c(0.01, 0.05, 0.20), 
                        labels=c("f = 0.01", "f = 0.05", "f = 0.20"))

#only include 20%off
pd <- position_dodge(0.0)
plot1 <- filter(result.plot, !grepl("sample|conditional|10%|5%|50%", method), grepl("0\\.01|0\\.2", f)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_grid(~f, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="odds ratio r") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette) +
  ggtitle('(a)') +
  theme(plot.title = element_text(hjust = 0))
plot1

#4.b

# The color-blind friednly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#check power using true, imputed founder carrier, minimum offspring carrier
#three versions -- command and rare and super rare #2 for common, 7 for rare, 39 for super rare
result <- read.csv("figure 4b. fouder imputatuion estimate nontransmitted by population f see power.csv", header=T)
result <- result %>% gather(key="method", value="p.value", -c(1:8))
result.plot <- result %>% group_by(f, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05, na.rm=T))
result.plot$method <- factor(result.plot$method, 
                             levels = c("result.trap.3c.noimpute", "result.trap.3c.impute.founder.assumed.f","result.trap.3c.impute.founder.pop.f"),
                             labels=c("w.o. imputation", "assumed.f" ,"pop.f"))
result.plot$f <- factor(result.plot$f, levels=c("0.005n3", "0.005n7"), 
                        labels=c("3 snps", "7 snps"))

#everything
pd <- position_dodge(0.0)
plot2 <- filter(result.plot, !grepl("conditional", method)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="odds ratio r") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85')) +
  ggtitle('(b) f = 0.005') +
  theme(plot.title = element_text(hjust = 0))
plot2

library(gridExtra)
grid.arrange(plot1, plot2, nrow=2, ncol=1)
