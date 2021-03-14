# R connect to a database

install.packages("RMySQL")
library(DBI)
con<- dbConnect(RmySQL::MySQL(),
                dbname = "Userdata.sql")
summary(con)
dbListTables(con) #list of table in the database

dbReadTable(con,"Userdata") #read table

dbGetQuery(con,"select * from Userdata")
