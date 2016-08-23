library(ggplot2)

#July 18, 2016
#Rocio Dominguez Vidana, PhD
#Springboard Foundations of Data Science
#ggplot2 code for Titanic dataset exercise


#get data from data wrangling exercise.2
load("titanic.RData")
data.clean<-data.clean[-1310,] #get rid of last line (blank line on csv)

ggplot(data.clean, aes(x = factor(pclass), fill = factor(sex))) + 
  geom_bar(position = "dodge")


ggplot(data.clean, aes(x = factor(pclass), fill = factor(sex))) + 
  geom_bar(position = "dodge") + 
  facet_grid(. ~ survived)

ggplot(data.clean, aes(x = factor(pclass), y = age, col = factor(sex))) + 
  geom_jitter(size = 3, alpha = 0.5, position = position_jitter(0.5, 0)) + 
facet_grid(. ~ survived)
