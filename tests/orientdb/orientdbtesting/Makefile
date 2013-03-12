.PHONY: all run

run: target/orientdbtesting-1.0-SNAPSHOT.jar
	mvn exec:java -D exec.mainClass=edu.gt.gtri.orientdbtesting.App -Dexec.args="$(ARGS)"

target/orientdbtesting-1.0-SNAPSHOT.jar:  src/main/java/edu/gt/gtri/orientdbtesting/*.java
	mvn package

all: run
