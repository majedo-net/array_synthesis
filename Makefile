

build:
	podman build -t array_synth:latest .

run:
	podman run --privileged -it -v /data/pso_results:/results array_synth:latest
