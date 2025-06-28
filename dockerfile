
# ใช้ base image ที่รองรับ conda และเหมาะกับ bioinformatics
FROM continuumio/miniconda3

# ตั้ง working directory
WORKDIR /recsalmo

# คัดลอกไฟล์ทั้งหมดของโปรเจกต์เข้า container
COPY . /recsalmo

# ติดตั้ง Java สำหรับโปรแกรมที่ต้องใช้
RUN apt-get update && \
 apt-get install -y default-jdk default-jre && \
 apt-get clean

# สร้าง conda environment และติดตั้ง dependencies
RUN conda create -n recsalmo_env python=3.8 -y && \
 echo "conda activate recsalmo_env" >> ~/.bashrc

# เปิดใช้งาน environment และติดตั้งแพ็กเกจ
RUN /bin/bash -c "source ~/.bashrc && \
 conda activate recsalmo_env && \
 conda install -c bioconda fastmlst sistr_cmd ncbi-amrfinderplus parsnp -y && \
 conda install -c conda-forge openpyxl seaborn -y && \
 fastmlst --update-mlst -t 1 && \
 amrfinder -u"

# ตั้งค่า entrypoint สำหรับการเรียกใช้งาน
ENTRYPOINT ["python", "recsalmo.py"]
