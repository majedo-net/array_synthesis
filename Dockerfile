FROM openems:latest
# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1
WORKDIR /
COPY setup.m /opt/openEMS/share/openEMS/matlab/setup.m
COPY hdf5setup.m .
RUN octave --silent hdf5setup.m 
COPY requirements.txt .
RUN pip3 install -r requirements.txt
COPY entry.py .
COPY spiral/spiral.m .
COPY spiral/spiral_ff.m .
COPY circ_rps.py .
COPY hex_rps.py .
COPY array_funcs.py .
COPY ArrayElementException.py .
ENTRYPOINT python3 entry.py
