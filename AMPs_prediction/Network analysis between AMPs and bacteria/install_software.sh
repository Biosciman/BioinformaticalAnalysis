#!/bin/bash
wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2/download -O samtools-1.7.tar.bz2
wget https://sourceforge.net/projects/samtools/files/samtools/1.7/htslib-1.7.tar.bz2/download -O htslib-1.7.tar.bz2
tar jxvf htslib-1.7.tar.bz2
tar jxvf samtools-1.7.tar.bz2

# Cenos, root
yum -y install bzip2-devel
yum -y install ncurses-libs
yum -y install ncurses-devel
yum -y install xz-devel.x86_64
yum -y install zlib-devel
yum -y install curl-devel
yum -y install git

# Ubuntu, root
apt-get install libbz2-dev
apt-get install zlib1g-dev
apt-get install liblzma-dev
apt-get install libncurses5-dev
apt-get install libcurl4-openssl-dev

git clone https://github.com/twestbrookunh/paladin.git
cd paladin/
make
PATH=$PATH:$(pwd)

cd htslib-1.7
./configure
make
make install

cd samtools-1.7
make clean
make
make prefix=/opt/samtools install
echo 'export PATH=$PATH:/opt/samtools/bin' >> /etc/profile
source /etc/profile
```
- Abundence caculation
```bash
# Create index
paladin index -r3 selected_c_AMPs.fa

# Compare to generate bam files. Takes the most time.
paladin align -t 4 selected_c_AMPs.fa test.fa | samtools view -Sb - > t1.bam

# Apply samtools to get the comparison file, re_ab.csv is the abundance of c_AMPs
samtools sort t1.bam -o s1.bam
samtools index s1.bam
samtools idxstats s1.bam >> re_ab.csv
```
