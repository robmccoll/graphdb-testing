.PHONY: all run

run: target/neo4jtesting-1.0-SNAPSHOT.jar
	mvn exec:java -D exec.mainClass=edu.gt.gtri.neo4jtesting.App -Dexec.args="$(ARGS)"

target/neo4jtesting-1.0-SNAPSHOT.jar:  src/main/java/edu/gt/gtri/neo4jtesting/*.java
	mvn package

all: run
