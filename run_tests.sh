# directories
testdir=tests
graphdir=graphs
resultsdir=results

# files
if [ $# -ge 1 ]; then
  graphs="$1"
else
  graphs="tiny small medium large"
fi

if [ $# -ge 2 ]; then
  tests="$2"
else
  tests=$(ls $testdir)
fi

run=$(date "+%Y.%m.%d.%H.%M.%S")
cwd=$(pwd)

echo "Starting with config:"
echo "  graphs: $graphs"
echo "  frameworks: $tests"
echo "  runID: $run"

for g in $graphs
  do
    for f in $tests
      do
	if [ -e $testdir/$f/runme.sh ]
	  then
	    echo "Running Graph: $g, Framework: $f Start Time: $(date '+%Y/%m/%d %H:%M:%S')..."
	    outfile=$cwd/$resultsdir/$run.$f.$g
	    sh sysinfo.sh > $outfile
	    cd $testdir/$f; sh runme.sh $cwd/$graphdir/$g.g $cwd/$graphdir/$g.a >> $outfile; cd $cwd
	    echo "  done."
	  fi
      done
  done
