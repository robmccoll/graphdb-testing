export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib/dexcpp-4.7.1/lib/linux64
./main $1 $2 | grep 'RSLT:' | sed -e 's/.*RSLT: //' | python -mjson.tool
