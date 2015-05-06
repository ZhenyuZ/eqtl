# Functions to access TCGA data
## GetSecureEntry.r
Ask the user to input TCGA username and password, and return TCGA secure site entry point (with username and password embeded in the URL) <br>
`GetSecureEntry()`

## GetTCGATable.r
Input local file name or url of a TCGA table file (sdrf or txt), and return the content as data frame <br>
`GetTCGATable(url)`
