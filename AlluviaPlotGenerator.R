library(ggalluvial)
library(Peptides)
library(runner)
library(tibble)
library(ggpubr)

# TestSequence <- "MWPTRRLVTIKRSGVDGPHFPLSLSTCLFGRGIECDIRIQLPVVSKQHCKIEIHEQEAILHNFSSTNPTQVNGSVIDEPVRLKHGDVITIIDRSFRYENMWPTRRLVTIKRSGVDGPHFPLSLSTCLFGRGIECDIRIQLPVVSKQHCKIEIHEQEAILHNFSSTNPTQVNGSVIDEPVRLKHGDVITIIDRSFRYEN"
TestSequence <- subset(human_data,`ACC#`=="O43521")$Sequence

testList<-runner(strsplit(TestSequence,split = "")[[1]],f = function(x){return(aaComp(paste(x,collapse = "")))},k=20)
testMat<-do.call(rbind,testList)
testMat<-cbind(testMat,rep(1:length(testList),each=9))
testData <- as.data.frame(testMat)
colnames(testData)<-c("Number","Percentage","Position")
testData <- rownames_to_column(testData,var = "Type")

testData$Type <- trimws(testData$Type,which = "right",whitespace = "\\.[0-9]+")


ggplot(data = testData,
       aes(x = Position, y = Percentage, alluvium = Type)) +
  geom_alluvium(aes(fill = Type, colour = Type),
                alpha = .75, decreasing = FALSE) +
  # scale_x_continuous(breaks = seq(2003, 2013, 2)) +
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=19),legend.position = "none") +
  scale_x_continuous(limits = c(0, nchar(TestSequence)),breaks = c(0,50,100,150,200),expand = c(0.005,0.005))+ xlab("Positions") +
  scale_y_continuous(expand = c(0.05,0.05)) + ylab("% Units") +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  scale_color_brewer(type = "qual", palette = "Set3")