

build:
	docker build -t array_synth:latest .

build-spiral-python:
	docker build -t spiral_python:latest spiral/.

build-tess:
	docker build -t tessarray:latest -f Dockerfile.tess-array .

build-ct:
	docker build -t ct:latest -f Dockerfile.coupling .

build-single:
	docker build -t single:latest -f Dockerfile.single .

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

run-tess:
	docker run -it --rm \
	--user 1000 \
	--net=host \
	--volume=/home/ubuntu/results/tess/y/NP5:/results \
	tessarray:latest

run-full:
	docker run -it --rm \
	--user 1000 \
	--net=host \
	--volume=/home/ubuntu/results/tess/spiral:/results \
	fullarray:latest

run-ct:
	docker run -it --rm \
		--user 1000 \
		--net=host \
		--volume=/home/ubuntu/results/ct/patch:/results \
		ct:latest

run-single:
	docker run -it --rm \
		--user 1000 \
		--net=host \
		--volume=/home/ubuntu/results/spiral:/results \
		single:latest
