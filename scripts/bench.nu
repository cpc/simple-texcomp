def bench [] { seq 10 | each {
    ./transcoder ~/data/combined_set/* ../test/combined | lines | where $it =~ 'Average encoding rate'
} | split column ':' -c | get Column2 | str to-decimal | math avg }
