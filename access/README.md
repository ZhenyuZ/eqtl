# Functions to access TCGA data
## GetTCGACredential.r
Ask for TCGA username and password, and return a data frame with $username as username and $password as password <br>
`GetTCGACredential()`

## GetTCGATable.r
Input local file name or url of a TCGA table file (sdrf or txt), and return the content as data frame <br>
`GetTCGATable(url)`
