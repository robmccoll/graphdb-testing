cd titantesting; make ARGS="$1 $2" | grep 'RSLT:' | sed -e 's/.*RSLT: //' | python -mjson.tool; cd ../
