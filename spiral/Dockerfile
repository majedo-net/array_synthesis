FROM openems:latest
WORKDIR "/"
COPY spiral.m ./
COPY entry.py ./
RUN pwd
RUN pip3 install pyswarms
ENTRYPOINT python3 entry.py 
