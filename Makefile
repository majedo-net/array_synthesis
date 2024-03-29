

build:
	docker build -t array_synth:latest .

build-spiral-python:
	docker build -t spiral_python:latest spiral/.

build-tessarray:
	docker build -t tessarray:latest -f Dockerfile.python .

build-array:
	docker build -t fullarray:latest -f Dockerfile.python-array .

run:
	docker run --privileged -it -v /data/pso_results:/results array_synth:latest

run-spiral-python:
	docker run -it --rm \
	--user 1000 \
	--net=host \
	--volume=C:\Users\mattj\Documents\CSM\array_synthesis\results:/results \
	spiral_python:latest

run-tessarray:
	docker run -it --rm \
	--user 1000 \
	--net=host \
	--volume=C:\Users\mattj\Documents\CSM\array_synthesis\results:/results \
	tessarray:latest

run-array:
	docker run -it --rm \
	--user 1000 \
	--net=host \
	--volume=C:\Users\mattj\Documents\CSM\array_synthesis\results:/results \
	fullarray:latest