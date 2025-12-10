FROM python:3

WORKDIR /assessor

COPY requirements.txt /assessor
RUN pip install --no-cache-dir -r requirements.txt

COPY lib /assessor/lib
COPY bin /assessor/bin

CMD ["bin/assess_mzMLs.py"]
