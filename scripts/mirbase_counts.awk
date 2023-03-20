BEGIN {FS = "\t"}
{printf $1"\t"$2"\t"}{for (i = 3; i <= NF; i += 4) printf ("%s%c", $i, i + 4 <= NF ? "\t" : "\n");}
