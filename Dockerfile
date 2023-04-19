# This is all python, so ...
FROM python:3.9.16

COPY . .

RUN pip install .

