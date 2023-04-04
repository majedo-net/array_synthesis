

build:
	podman build -t array_synth:plots .

run:
	podman run --privileged -it -v /data/pso_results:/results array_synth:plots
