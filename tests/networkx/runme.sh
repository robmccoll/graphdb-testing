python test_python.py $1 $2 | grep 'RSLT:' | sed -e 's/.*RSLT: //' | python -mjson.tool
