source('animate_functions.R')

dbxdir='~/gdrive/dunlop/razan_animation'
deadfile=file.path(dbxdir, 'gadXrecA-combdeadcells.csv')
livefile=file.path(dbxdir, 'gadXrecA-comblivecells.csv')

plot=all.steps.plot(livefile, deadfile, 512)
