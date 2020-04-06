library("RMySQL")
library(dplyr)  


# Database info  -------------------------------------------------------


db_pr2 <- list( user='pr2_read', 
                  password='XjerZ1yX5d7tVjkQTYCW', 
                  dbname='pr2_database', 
                  host='62.210.222.245')

# Database function -------------------------------------------------------

get_query <- function(database_info, query)  { 
    db <- dbConnect(MySQL(),  user=database_info$user, 
                              password=database_info$password, 
                              dbname=database_info$dbname, 
                              host=database_info$host)
    table <- dbGetQuery(db,query )
    dbDisconnect(db)
    return(table)
}


# PR2 read function -------------------------------------------------------

pr2_read <- function() {

# Read the PR2 full tables including removed records

  pr2_main <- get_query(db_pr2, "select * from pr2_main") 
  pr2_seq <- get_query(db_pr2, "select * from pr2_sequences") 
  pr2_taxo <- get_query(db_pr2, "select * from taxo") 
  pr2_metadata <- get_query(db_pr2, "select * from pr2_metadata") 
 
# Join the tables and keep only sequences that are not removed 
  
  pr2 <- left_join(pr2_main, pr2_taxo, by = c("species"="species"))
  pr2 <- left_join (pr2, pr2_seq) 
  pr2 <- left_join (pr2, pr2_metadata)
  return(pr2)
}

# Load PR2 into dataframe -------------------------------------------------

pr2 <- pr2_read()