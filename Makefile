

build:
	docker build -t array_synth:latest .

build-spiral-python:
	docker build -t spiral_python:latest spiral/.

run:
	docker run --privileged -it -v /data/pso_results:/results array_synth:latest

run-spiral-python:
	docker run --it --rm \
	--user 1000 \
	--net=host \
	--volume=C:\Users\mattj\Documents\CSM\array_synthesis\results:/results \

