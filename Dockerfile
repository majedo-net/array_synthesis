FROM 10.2.10.230:5000/openems:v1.1
WORKDIR "/"
COPY spiral.m ./
COPY entry.py ./
COPY hex_rps.py ./
COPY array_funcs.py ./
RUN pwd
RUN pip3 install pyswarms
ENTRYPOINT python3 entry.py 
