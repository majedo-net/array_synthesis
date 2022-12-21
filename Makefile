

build:
	docker build -t array_synth:latest .

run:
	docker run --privileged -it -v C:\Users\mattj\Documents\CSM\array_synthesis\results:/results array_synth:latest