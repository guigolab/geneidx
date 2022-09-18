FROM python:latest
WORKDIR /docker
ADD ./files /files
ADD ./Parameter_files.taxid /data/Parameter_files.taxid
RUN python -m pip install --upgrade pip
RUN pip install requests
RUN pip install lxml
RUN pip install pandas
RUN pip install numpy
ENTRYPOINT [ "python" ]
