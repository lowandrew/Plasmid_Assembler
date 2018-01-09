# Customizing the Plasmid Database

The plasmid database mentioned in the [Installation](installation.md) section should be suitable for most purposes, as it contains 
more than 9000 RefSeq plasmids covering all Bacteria. However, should you wish to add more sequences or use a completely different database,
this is entirely possible.

#### Adding on to the Default Database

Adding on to the default database is simple. To do so, just append any other sequences you would like PlasmidExtractor to search for to the end of the `plasmid_db.fasta` file
included in the database download. This is particularly easily done from the command line. If, for example, you have a new sequence you would like to add to the database called
`new_sequence.fasta` in the same folder as `plasmid_db.fasta`, a simple `cat` command can append for you:

`cat new_sequence.fasta >> plasmid_db.fasta` 

#### Using A Different Database

Using a different database is also possible and easy to do. A valid plasmid database is any FASTA-formatted (uncompressed) file with one or more sequences within it, with unique names for each.
To use a different database, just change the file you're pointing to with the `-p` option when calling PlasmidExtractor. 

# Using Custom AMR/Virulence/Incompatibility Group Databases

The defaults databases for AMR/Virulence/Incompatibility Group Detection included with PlasmidExtractor
are derived from [CGE Databases](https://bitbucket.org/account/user/genomicepidemiology/projects/CGE). Given that there are approximately 8 million different AMR databases out there, it's entirely possible that you have your own favorite that you would like to use. With PlasmidExtractor, you are able to do this with a fairly simple procedure.

#### Adding on to the Defaults

If you want to add on to the default database, you can just add your sequences of interest to the relevant folder, with the caveat that the file extension needs to be `.tfa`. So, for example, if you have 
a FASTA file with new AMR genes called `amr.fasta` that you would also like searched, you would move (or copy) the file to `databases/resistance_db` and rename it to `amr.tfa`. The following command would do this for you:

`cp amr.fasta databases/resistance_db/amr.tfa`

#### Creating Your Own Database

In the event you don't want anything to do with the included AMR databases, you don't have to use them at all! You can just delete the file in whichever folder you want to use your own databases for, 
and then put the FASTA files for your own database (with the `.tfa` extension!) into that folder, and you'll be good to go. Note that the names of the folders for each type of detection must remain the same, and at least one target file must be in each of the folders, or you will have an error.
