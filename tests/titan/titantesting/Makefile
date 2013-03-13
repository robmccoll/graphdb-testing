.PHONY: all run

run: target/titantesting-1.0-SNAPSHOT.jar
	mvn exec:java -D exec.mainClass=edu.gt.gtri.titantesting.App -Dexec.args="$(ARGS)"

target/titantesting-1.0-SNAPSHOT.jar:  src/main/java/edu/gt/gtri/titantesting/*.java
	mvn package

all: run
