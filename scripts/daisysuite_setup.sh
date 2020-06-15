# Setup and run daisySuite for HGT inference from reads using the docker container.
# cmdoret, 20190527

DIR=/data/legionella/fq/neff_shotgun/daisy

docker run -e LD_LIBRARY_PATH=/usr/local/lib -ti -d \
           --name daisy -v $DIR:/daisy \
           quay.io/biocontainers/daisysuite:1.3.0--0

# Symlinking broken library for pysam in the container
docker exec daisy ln -s /usr/local/lib/libhts.so.1.9 /usr/local/lib/libhts.so.1

# 