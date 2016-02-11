
# The color-blind friednly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#check power using true, imputed founder carrier, minimum offspring carrier
f_1 <- 0.01
f_2 <- 0.05
f_3 <- 0.20

mu_1 <- seq(1,3, length.out = 1000)
mu_2 <- seq(1,2.2, length.out = 1000)
mu_3 <- seq(1,1.6, length.out = 1000)

result.plot <- rbind(data.frame(f=f_1, S=1, mu=mu_1, f_untrans=E_f_untrans(f=f_1, mu=mu_1, sigma2 = 0, S = 1)), 
                data.frame(f=f_2, S=1, mu=mu_2, f_untrans=E_f_untrans(f=f_2, mu=mu_2, sigma2 = 0, S = 1)),
                data.frame(f=f_3, S=1, mu=mu_3, f_untrans=E_f_untrans(f=f_3, mu=mu_3, sigma2 = 0, S = 1)),
                data.frame(f=f_1, S=2, mu=mu_1, f_untrans=E_f_untrans(f=f_1, mu=mu_1, sigma2 = 0, S = 2)), 
                data.frame(f=f_2, S=2, mu=mu_2, f_untrans=E_f_untrans(f=f_2, mu=mu_2, sigma2 = 0, S = 2)),
                data.frame(f=f_3, S=2, mu=mu_3, f_untrans=E_f_untrans(f=f_3, mu=mu_3, sigma2 = 0, S = 2)))
result.plot <- result.plot %>% group_by(f,S,mu)
result.plot$f <- factor(result.plot$f, labels=c("f=0.01", "f=0.05", "f=0.20"))
result.plot$S <- factor(result.plot$S)
#using weight in TRAP
pd <- position_dodge(0.0)
filter(result.plot) %>% ggplot(aes(x=mu, y=f_untrans, ymax=max(f_untrans), group=S, col=S)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk r", y="f for non-transmitted") +
#   scale_y_continuous(limits=c(0,1)) +
  ggtitle("allele frequency for non-transmitted chromosome by population allele frequency") +
  theme_gray(base_size = 20) +
  theme(legend.position="bottom", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette)

