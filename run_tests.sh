# directories
testdir=tests
graphdir=graphs
resultsdir=results

# files
graphs="tiny small medium large"
tests=$(ls $testdir)

run=$(date "+%Y.%m.%d.%H.%M.%S")
cwd=$(pwd)

for g in $graphs
  do
    for f in $tests
      do
	if [ -e $testdir/$f/runme.sh ]
	  then
	    echo "Running Graph: $g, Framework: $f..."
	    outfile=$cwd/$resultsdir/$run.$f.$g
	    sh sysinfo.sh > $outfile
	    cd $testdir/$f; sh runme.sh $cwd/$graphdir/$g.g $cwd/$graphdir/$g.a >> $outfile; cd $cwd
	    echo "-->done."
	  fi
      done
  done
