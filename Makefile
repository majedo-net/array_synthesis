

build:
	docker build -t array_synth:latest .

build-spiral-python:
	docker build -t spiral_python:latest spiral/.

build-tessarray:
	docker build -t tessarray:latest -f Dockerfile.python .

build-array:
	docker build -t fullarray:latest -f Dockerfile.python-array .

build-full:
	docker build -t fullarray:latest -f Dockerfile.full-array .

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
	--volume=C:\Users\mattj\Documents\CSM\array_synthesis\results\aws-tess-array\nr12:/results \
	tessarray:latest

run-full:
	docker run -it --rm \
	--user 1000 \
	--net=host \
	--volume=C:\Users\mattj\Documents\CSM\array_synthesis\results\aws-tess-array\dipole\full:/results \
	fullarray:latest