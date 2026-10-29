open (IN, "$ARGV[0]") || die "nope\n";
while (<IN>)
{
    chomp;
    @fields = split /\t/;
    $P3_ST = $fields[1]-$ARGV[4];
    $LenName = $ARGV[5];
    `blastdbcmd -dbtype nucl -db $ARGV[1] -entry $fields[0] -range 1-$ARGV[4] -out $fields[0].5p.$LenName.fa`;
    `blastdbcmd -dbtype nucl -db $ARGV[1] -entry $fields[0] -range $P3_ST-$fields[1] -out $fields[0].3p.$LenName.fa`;
    `blastdbcmd -dbtype nucl -db $ARGV[1] -entry $fields[0] -range $ARGV[4]-$P3_ST -out $fields[0].InterTelo.fa`;
    `jellyfish count -m $ARGV[3] -s $ARGV[2] -t 4 -o $fields[0].5p.$LenName.k$ARGV[3] -C -L 1 $fields[0].5p.$LenName.fa`;
    `echo "Mer\t$fields[0]_5p" > "sorted.$fields[0].5p.$LenName.k$ARGV[3].dump"`;
    `jellyfish dump -c -t $fields[0].5p.$LenName.k$ARGV[3] > $fields[0].5p.$LenName.k$ARGV[3].dump `;
    `sort -k1,1 $fields[0].5p.$LenName.k$ARGV[3].dump >> sorted.$fields[0].5p.$LenName.k$ARGV[3].dump`;
    `jellyfish count -m $ARGV[3] -s $ARGV[2] -t 4 -o $fields[0].3p.$LenName.k$ARGV[3] -C -L 1 $fields[0].3p.$LenName.fa`;
    `echo "Mer\t$fields[0]_3p" > sorted.$fields[0].3p.$LenName.k$ARGV[3].dump`;
    `jellyfish dump -c -t $fields[0].3p.$LenName.k$ARGV[3] > $fields[0].3p.$LenName.k$ARGV[3].dump`;
    `sort -k1,1 $fields[0].3p.$LenName.k$ARGV[3].dump >> sorted.$fields[0].3p.$LenName.k$ARGV[3].dump`;
    `jellyfish count -m $ARGV[3] -s $ARGV[2] -t 4 -o $fields[0].InterTelo_k$ARGV[3] -C -L 1 $fields[0].InterTelo.fa`;
    `echo "Mer\t$fields[0]_IT" > sorted.$fields[0].InterTelo.$LenName.k$ARGV[3].dump`;
    `jellyfish dump -c -t $fields[0].InterTelo_k$ARGV[3] > $fields[0].InterTelo.$LenName.k$ARGV[3].dump`;
    `sort -k1,1 $fields[0].InterTelo.$LenName.k$ARGV[3].dump >> sorted.$fields[0].InterTelo.$LenName.k$ARGV[3].dump`;
}
close (IN);
