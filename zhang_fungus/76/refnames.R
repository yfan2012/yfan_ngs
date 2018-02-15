library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))


##Read some name/species info from a google sheet and change the assembly names accordingly

gs_auth(token = "~/.ssh/googlesheets_token.rds")
names=gs_url('https://docs.google.com/spreadsheets/d/1a5ztAbbYKFLM52OvW8ikt7dzqZGJtbxQCOtJaR9V2fc/edit#gid=0')
data=gs_read(names)

data$id=gsub(' ', '' ,data$samp)
data$new=paste(data$name, data$id, sep='_')

gs_edit_cells(names, input=data, anchor='A1')


##rename assemblies of reference 
for (i in 1:length(data$id)) {
    parsnp='/dilithium/Data/NGS/Aligned/170928_fungus76/assemble/parsnp_assemblies'
    system(paste0('mv ', parsnp, '/', data$id[i], '.fasta ', parsnp, '/', data$new[i], '.fasta'))
}


