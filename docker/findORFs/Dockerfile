FROM python:latest
WORKDIR /docker
ADD ./files /files
RUN python -m pip install --upgrade pip
RUN pip install requests
RUN pip install lxml
RUN pip install pandas
RUN pip install numpy
ENTRYPOINT [ "python" ]
