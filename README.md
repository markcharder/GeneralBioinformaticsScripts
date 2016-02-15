# GeneralReformattingScripts

This repository includes the following kinds of scripts:

- Scripts used to reformat typical sequence data files such as fasta files and gff files.
- Scripts to check validity of typical sequence data files.

Most of the reformatting scripts are semi-portable, i.e. they may require some pre-emptive formatting by the user.
They should take care of the tricky bits that would be impractical to attempt to do manually.
Note: checkStops.pl is for use on a peptide sequence file that should include stops, otherwise it will report missing stops for all sequences within the file.
