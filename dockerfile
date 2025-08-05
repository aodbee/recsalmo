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
    conda run -n recsalmo_env conda install -c bioconda fastmlst sistr_cmd ncbi-amrfinderplus parsnp -y && \
    conda run -n recsalmo_env conda install -c conda-forge openpyxl seaborn biopython -y && \
    conda run -n recsalmo_env fastmlst --update-mlst -t 1 && \
    conda run -n recsalmo_env amrfinder -u

# ใช้ shell แบบ conda run เพื่อเรียกใช้งาน script
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "recsalmo_env", "python", "recsalmo.py"]
