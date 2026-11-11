tail -n +2 $1 | split -l $2 -d -a 3 - split.$3.
for file in split.*
do
    head -n 1 $1 > tmp_file
    cat "$file" >> tmp_file
    mv -f tmp_file "$file"
done
