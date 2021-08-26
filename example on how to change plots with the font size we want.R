library(ggpubr)
# ToothGrowth
data("ToothGrowth")
head(ToothGrowth)
# mtcars 
data("mtcars")
mtcars$name <- rownames(mtcars)
mtcars$cyl <- as.factor(mtcars$cyl)
head(mtcars[, c("name", "wt", "mpg", "cyl")])

# Box plot (bp)
bxp <- ggboxplot(ToothGrowth, x = "dose", y = "len",
                 color = "dose", palette = "jco")
bxp
# Dot plot (dp)
dp <- ggdotplot(ToothGrowth, x = "dose", y = "len",
                color = "dose", palette = "jco", binwidth = 1)
dp


# Bar plot (bp)
bp <- ggbarplot(mtcars, x = "name", y = "mpg",
                fill = "cyl",               # change fill color by cyl
                color = "white",            # Set bar border colors to white
                palette = "jco",            # jco journal color palett. see ?ggpar
                sort.val = "asc",           # Sort the value in ascending order
                sort.by.groups = TRUE,      # Sort inside each group
                x.text.angle = 90           # Rotate vertically x axis texts
)
bp + font("x.text", size = 8)
# Scatter plots (sp)
sp <- ggscatter(mtcars, x = "wt", y = "mpg",
                add = "reg.line",               # Add regression line
                conf.int = TRUE,                # Add confidence interval
                color = "cyl", palette = "jco", # Color by groups "cyl"
                shape = "cyl"                   # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = cyl), label.x = 3)       # Add correlation coefficient
sp

ggarrange(sp, bp + font("x.text", size = 30),
          ncol = 1, nrow = 2)
ggarrange(sp, bp + font("x.text", size = 4),
          ncol = 1, nrow = 2)
ggarrange(sp + font("x.text", size = 30), bp + font("x.text", size = 10),
          ncol = 1, nrow = 2)
           