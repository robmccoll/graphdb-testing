.PHONY: all run

run: target/infogridtesting-1.0-SNAPSHOT.jar
	mvn exec:java -D exec.mainClass=edu.gt.gtri.infogridtesting.App -Dexec.args="$(ARGS)"

target/infogridtesting-1.0-SNAPSHOT.jar:  src/main/java/edu/gt/gtri/infogridtesting/*.java
	mvn package

all: run
